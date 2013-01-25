import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from util import indices_dict, each_replicate
from config import All_samples, Gene_count_matrix_file_path, Gene_indices_file_path, Sample_plots
import csv

def load_gene_count_matrix():
	gene_count_matrix_file = open(Gene_count_matrix_file_path, 'rb')
	gene_count_matrix = np.load(gene_count_matrix_file, )
	gene_count_matrix_file.close()
	return gene_count_matrix 

def load_gene_indices():
	gene_indices_file = open(Gene_indices_file_path , 'rb')
	gene_indices = pickle.load(gene_indices_file)
	gene_indices_file.close()
	return gene_indices 

def calc_stats(samples, replicate_indices, gene_count_matrix, gene_indices):
	sample_stats = dict()
	for sample, replicates in samples.iteritems():
		sample_matrix = np.column_stack((gene_count_matrix[:,replicate_indices[replicate]] for replicate in replicates))
		means = np.mean(sample_matrix, axis=1)
		stddevs = np.std(sample_matrix, axis=1)
		sample_stats[sample] = (means, stddevs)
	return sample_stats

def indices_to_list(indices):
	list_ = [0] * len(indices)
	for (elem, index) in indices.iteritems():
		list_[index] = elem
	return list_

def make_plots(sample_stats, plots, gene_list):
	for sample1, sample2 in plots:
		means1, stddevs1 = sample_stats[sample1]
		means2, stddevs2 = sample_stats[sample2]
		diffmeans = means1 - means2
		diffstddevs = ((means1 ** 2. / len(means1)) + (means2 ** 2. / len(means2))) ** .5
		plot = [(gene_list[i], diffmeans[i], diffstddevs[i]) for i in xrange(len(diffmeans))]
		plot.sort(key=lambda x: x[1])
		with open('%s_vs_%s.csv' % (sample1, sample2), 'w') as output:
			writer = csv.writer(output)
			writer.writerows(plot)
		

def main(samples):
	replicate_indices = indices_dict(each_replicate(samples))
	gene_count_matrix = load_gene_count_matrix()
	gene_indices = load_gene_indices()
	gene_list = indices_to_list(gene_indices)
	sample_stats = calc_stats(samples, replicate_indices, gene_count_matrix, gene_indices)
	make_plots(sample_stats, Sample_plots, gene_list)


main(All_samples)
