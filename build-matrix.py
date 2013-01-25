import os.path
import re
import itertools
import numpy as np

# -------------------------------------
#   data input
# -------------------------------------

All_samples = {
	'delta-lin41': (
			'Delta-Lin41-BC08_S3_L001_R1_001',
			'Delta-Lin41-BC08_S3_L001_R1_002',
			'Delta-Lin41-BC08_S3_L001_R1_003'),
	'lin41-gran': (
			'Lin41-Gran-BC04_S1_L001_R1_001',
			'Lin41-Gran-BC04_S1_L001_R1_002',
			'Lin41-Gran-BC04_S1_L001_R1_003'),
	'lin41-wc': (
			'Lin41-WC-BC06_S2_L001_R1_001',
			'Lin41-WC-BC06_S2_L001_R1_002',
			'Lin41-WC-BC06_S2_L001_R1_003'),
}

def annotated_sam_file_path(replicate_name):
	return '%s/tophat_out/accepted_hits-annotated.sam' % replicate_name

Genes_ignore = set(('alignment_not_unique', 'no_feature'))

# -------------------------------------
#  general util funcs
# -------------------------------------

def close_files(files):
	for file_ in files.itervalues():
		file_.close()

def rewind_files(files):
	for file_ in files.itervalues():
		file_.seek(0)

# return a sorted list of an iterator of elements
def sorted_list(elems):
	l = list(elems)
	l.sort()
	return l

# given a list of unique elements, return a dict of each element and its index in the list
def indices_dict(list_):
	return dict(zip(list_, itertools.count(0)))

# return the number of elements in an iterator
def iterlen(it):
	return sum(1 for _ in it)

# -------------------------------------
#  sample util funcs
# -------------------------------------

def check_files_in_samples(samples):
	all_files_found = True
	for (sample, replicates) in samples.iteritems():
		for replicate in replicates:
			if not os.path.exists(replicate):
				print 'Could not find: %s' % replicate
				all_files_found = False
	return all_files_found

def each_replicate(samples):
	return itertools.chain(*samples.itervalues())

# -------------------------------------
# working with annotated sam files
# -------------------------------------

def open_annot_sam_files(samples):
	return dict((replicate, open(annotated_sam_file_path(replicate), 'r')) for replicate in each_replicate(samples))

Annotated_sam_re = re.compile(r'XF:Z:(.+)$')
Ambiguous_prefix = 'ambiguous['
Ambiguous_method_ignore = 0
Ambiguous_method_each = 1
Ambiguous_method_call_first = 2

def parse_annotated_sam(input_file, ambiguous_method):
	for line in input_file:
		line = line.strip()
		matched = Annotated_sam_re.search(line)
		if not matched: continue # WARNING: make sure Annotated_sam_re matches *every* line in input!
		name = matched.group(1)
		if not name.startswith(Ambiguous_prefix):
			yield name
		else:
			names = name[len(Ambiguous_prefix):-1].split('+')
			if ambiguous_method == Ambiguous_method_ignore:
				continue
			elif ambiguous_method == Ambiguous_method_each:
				for name in names:
					yield name
			elif ambiguous_method  == Ambiguous_method_call_first:
				yield names[0]


# -------------------------------------
#  gene count matrix
# -------------------------------------

def build_all_genes_set(annot_sam_files, ambiguous_method):
	all_genes = itertools.chain(*(parse_annotated_sam(annot_sam_file, ambiguous_method) for annot_sam_file in annot_sam_files.itervalues()))
	return set(all_genes) - Genes_ignore

def build_gene_indices(annot_sam_files, ambiguous_method):
	all_genes_set = build_all_genes_set(annot_sam_files, ambiguous_method)
	all_genes_list = sorted_list(all_genes_set)
	gene_indices = indices_dict(all_genes_list)
	return gene_indices 


def make_gene_count_matrix(samples, ambiguous_method):
	replicate_indices = indices_dict(each_replicate(samples))

	annot_sam_files = open_annot_sam_files(samples)
	gene_indices = build_gene_indices(annot_sam_files, ambiguous_method)
	rewind_files(annot_sam_files)

	gene_count_matrix = np.zeros((len(gene_indices), len(replicate_indices)), np.uint32)

	for replicate in each_replicate(samples):
		col = replicate_indices[replicate]
		annot_sam_file = annot_sam_files[replicate]
		for gene in parse_annotated_sam(annot_sam_file, ambiguous_method):
			if not gene in gene_indices: continue
			row = gene_indices[gene]
			gene_count_matrix[row,col] += 1

	close_files(annot_sam_files)

	return (gene_count_matrix, gene_indices, replicate_indices)



def main(samples, ambiguous_method):
	if not check_files_in_samples(samples):
		return

	(gene_count_matrix, gene_indices, replicate_indices) = make_gene_count_matrix(samples, ambiguous_method)

	print ' '.join(replicate_indices.iterkeys())
	count = 0
	for (gene, row) in gene_indices.iteritems():
		count += 1
		if count == 1000: break
		for (replicate, col) in replicate_indices.iteritems():
			print gene_count_matrix[row,col],
		print

main(All_samples, Ambiguous_method_each)
