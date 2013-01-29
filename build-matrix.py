import re
import itertools
import numpy as np
import cPickle as pickle
from config import All_samples, annotated_sam_file_path, Genes_ignore, Gene_count_matrix_file_path, Gene_indices_file_path 
from util import close_files, rewind_files, sorted_list, indices_dict, check_files_in_samples, each_replicate

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
	gene_set = build_all_genes_set(annot_sam_files, ambiguous_method)
	gene_list = sorted_list(gene_set)
	gene_indices = indices_dict(gene_list)
	return (gene_list, gene_indices)


def make_gene_count_matrix(samples, ambiguous_method):
	annot_sam_files = open_annot_sam_files(samples)
	(gene_list, gene_indices) = build_gene_indices(annot_sam_files, ambiguous_method)
	rewind_files(annot_sam_files)

	replicate_list = list(each_replicate(samples))
	replicate_indices = indices_dict(replicate_list)

	gene_count_matrix = np.zeros((len(gene_list), len(replicate_list)), np.uint32)
	for (replicate, col) in replicate_indices.iteritems():
		annot_sam_file = annot_sam_files[replicate]
		for gene in parse_annotated_sam(annot_sam_file, ambiguous_method):
			if not gene in gene_indices: continue
			row = gene_indices[gene]
			gene_count_matrix[row,col] += 1

	close_files(annot_sam_files)

	return (gene_count_matrix, gene_list, replicate_list)

def main(samples, ambiguous_method):
	if not check_files_in_samples(samples):
		return
	
	(gene_count_matrix, gene_list, replicate_list) = make_gene_count_matrix(samples, ambiguous_method)

	gene_count_matrix_file = open(Gene_count_matrix_file_path, 'wb')
	np.save(gene_count_matrix_file, gene_count_matrix)
	gene_count_matrix_file.close()

	gene_list_file = open(Gene_indices_file_path , 'wb')
	pickle.dump(gene_list, gene_list_file, pickle.HIGHEST_PROTOCOL)
	gene_list_file.close()

main(All_samples, Ambiguous_method_ignore)
