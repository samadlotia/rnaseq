import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from util import indices_dict, each_replicate
from config import All_samples, Gene_count_matrix_file_path, Gene_indices_file_path, Sample_plots
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

def any_zeros(row):
	for i in range(len(row)):
		if row[i] == 0:
			return True
	return False

def calc_stats(samples, replicate_indices, gene_count_matrix, gene_list):
	sample_stats = dict()
	for sample, replicates in samples.iteritems():
		sample_matrix = np.column_stack((gene_count_matrix[:,replicate_indices[replicate]] for replicate in replicates))
		means = np.mean(sample_matrix, axis=1)
		#stddevs = np.std(sample_matrix, axis=1)
		logmeans = np.log(means)
		sample_stats[sample] = (sample_matrix, means, logmeans)
	return sample_stats


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
		

def output_matplotlib(sample_stats, gene_list, plots):
	for sampleA, sampleB in plots:
		sample_matrixA, _, logmeansA, _ = sample_stats[sampleA]
		sample_matrixB, _, logmeansB, _ = sample_stats[sampleB]
		for i in xrange(len(gene_list)):
			gene = gene_list[i]

def main(samples):
	replicate_indices = indices_dict(each_replicate(samples))
	gene_count_matrix = load_gene_count_matrix()
	gene_list = load_gene_list()
	sample_stats = calc_stats(samples, replicate_indices, gene_count_matrix, gene_list)
	output_tables(samples, sample_stats, gene_list)
	output_matplotlib(sample_stats, gene_list, Sample_plots)

main(All_samples)
