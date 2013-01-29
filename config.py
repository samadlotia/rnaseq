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

Gene_count_matrix_file_path = 'gene-count-matrix.npy'
Gene_indices_file_path = 'gene-indices.pickle'
