# Usage: compare_rafs.py <tree_dir>.
# Will write out CSV files to current working directory.

import cPickle
import sys
import os
import numpy as np
import tssb

def process(tree, tree_score):
	root = tree.root['node']
	node = root.children()[0]
	ssms = node.get_data()
	with open('%s.csv' % tree_score, 'w') as outf:
		outf.write(','.join(('ssm_name', 'calculated_vaf', 'expected_vaf')) + '\n')
		for ssm in ssms:
			for ssm_name, calculated_raf, expected_raf in calc_raf(tree_score, ssm, ssm.mu_r, ssm.mu_v, 0):
				calculated_vaf, expected_vaf = 1 - calculated_raf, 1 - expected_raf
				outf.write(','.join((ssm_name, str(calculated_vaf), str(expected_vaf))) + '\n')

def calc_raf(tree_score, ssm, mu_r, mu_v, tp, new_state=0):	
	poss_n_genomes = ssm.compute_n_genomes(tp,new_state)
	poss_n_genomes = [x for x in poss_n_genomes if x[1] > 0]
	print len(poss_n_genomes)
	for (nr,nv) in poss_n_genomes:
		mu = (nr * mu_r + nv*(1-mu_r) ) / (nr+ nv)
		yield (ssm.name, mu, float(ssm.a[0]) / ssm.d[0])

def main():
	fdir = sys.argv[1]
	flist = os.listdir(fdir)
	if sum([flist[i]=='.DS_Store' for i in np.arange(len(flist))]): flist.remove('.DS_Store')
	for idx, fname in enumerate(flist):
		with open(fdir + '/' + fname) as f:
			tree = cPickle.load(f)
			tree_score = float(fname)
			process(tree, tree_score)

main()
