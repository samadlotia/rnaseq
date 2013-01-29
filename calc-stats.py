import cPickle as pickle
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from util import indices_dict, each_replicate
from config import All_samples, Gene_count_matrix_file_path, Gene_indices_file_path
import csv

def load_gene_count_matrix():
	gene_count_matrix_file = open(Gene_count_matrix_file_path, 'rb')
	gene_count_matrix = np.load(gene_count_matrix_file, )
	gene_count_matrix_file.close()
	return gene_count_matrix 

def load_gene_list():
	gene_list_file = open(Gene_indices_file_path , 'rb')
	gene_list = pickle.load(gene_list_file)
	gene_list_file.close()
	return gene_list 

# given a row of numbers, returns True if any of its elements are zero
def any_zeros(row):
	for i in range(len(row)):
		if row[i] == 0:
			return True
	return False

Ensg_re = re.compile(r'ENSG0+(\d+)')
# simplifies ENSG ids; ex: 'ENSG000000001231' => '1231'
# useful for plotting  ENSG ids
def clean_ensg_name(fullname):
	return Ensg_re.match(fullname).group(1)

# given a matrix, where columns are replicates and rows are genes,
# return a new matrix with values normalized by column
def normalize_by_replicate(sample_matrix):
	_, numcols = sample_matrix.shape
	normcols = list()
	for x in xrange(numcols):
		col = sample_matrix[:,x]
		total = float(np.sum(col))
		normcol = col / total
		normcol = np.atleast_2d(normcol).transpose()
		normcols.append(normcol)
	norm_matrix = np.hstack(normcols)
	return norm_matrix

# given all the gene counts across all conditions and replicates,
# return a dictionary of each condition mapping its normalized
# count values, the means per gene, and the logarithmic mean.
# logarithmic means are increased by 7 to make the range of means
# start at around 0 instead of -7. if the mean is 0, the
# logarithmic mean is mapped to 0.
def calc_stats(samples, replicate_indices, gene_count_matrix):
	sample_stats = dict()
	for sample, replicates in samples.iteritems():
		sample_matrix = np.column_stack((gene_count_matrix[:,replicate_indices[replicate]] for replicate in replicates))
		norm_matrix = normalize_by_replicate(sample_matrix)
		means = np.mean(norm_matrix, axis=1)
		logmeans = np.log10(means)
		logmeans += 7.0
		logmeans[np.isinf(logmeans)] = 0.0
		sample_stats[sample] = (norm_matrix, means, logmeans)
	return sample_stats

# writes out csv files for each condition
def output_tables(samples, sample_stats, gene_list):
	for sample, (sample_matrix, means, logmeans) in sample_stats.iteritems():
		replicates = samples[sample]
		with open('%s.csv' % sample, 'w') as output:
			writer = csv.writer(output)
			header = ['gene']
			header.extend(replicates)
			header.extend(['mean', 'logmean'])
			writer.writerow(header)

			for i in range(len(gene_list)):
				row = [gene_list[i]]
				row.extend(sample_matrix[i,:])
				row.append(means[i])
				row.append(logmeans[i])
				writer.writerow(row)
		
# generates a plot of the logmeans comparing two conditions
# acceptable values for sampleA and sampleB parameters:
# 'delta-lin41', 'lin41-gran', 'lin41-wc'
def output_matplotlib(sample_stats, gene_list, sampleA, sampleB):
	_, _, logmeansA = sample_stats[sampleA]
	_, _, logmeansB = sample_stats[sampleB]
	plt.xlabel(sampleA)
	plt.ylabel(sampleB)
	for i in xrange(len(gene_list)):
		gene = gene_list[i]
		logmeanA = logmeansA[i]
		logmeanB = logmeansB[i]
		plt.plot(logmeanA, logmeanB, marker='.', label=gene)
		#plt.text(logmeanA, logmeanB, gene, fontsize=4, alpha=0.3, va='bottom', ha='left')
		plt.show()

def main(samples):
	replicate_indices = indices_dict(each_replicate(samples))
	gene_count_matrix = load_gene_count_matrix()
	gene_list = map(clean_ensg_name, load_gene_list())
	sample_stats = calc_stats(samples, replicate_indices, gene_count_matrix)
	#output_tables(samples, sample_stats, gene_list)
	output_matplotlib(sample_stats, gene_list, 'delta-lin41', 'lin41-gran')

main(All_samples)
