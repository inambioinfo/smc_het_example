import tssb
from numpy import *
import cPickle
import os
import sys

def get_n_clusters(outputdir,min_genes=3):
	flist = os.listdir(outputdir+'/best/')
	ns = len(flist) #number of MCMC samples
	nclusters = zeros(ns)
	for idx,fname in enumerate(flist):
		f=open(outputdir+'/best/'+str(fname))
		tssb = cPickle.load(f)
		f.close()
		wts, nodes = tssb.get_mixture()
		nclusters[idx] = sum( [ 1 for node in nodes if len(node.get_data()) > min_genes ] )

	return nclusters

if __name__ == '__main__':
	fname = sys.argv[1]
	nclusters = get_n_clusters(fname+'.output',min_genes=100)
	savetxt(fname+'.nclusters.txt',nclusters)

