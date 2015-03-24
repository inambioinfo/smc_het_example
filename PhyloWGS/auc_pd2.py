# code to generate partial order plots

import cPickle

from numpy	  import *
from numpy.random import *
from tssb		import *
from util		import *
from util2 import *

#from ete2 import *

import sklearn.metrics as mt

#from sklearn.metrics import precision_recall_curve, roc_curve, auc

from clustering import *


import itertools
import string 
import numpy as np

def compute_auc(fdir,dname):	
	f = open("vafs.csv")
	d = f.readlines()
	f.close()

	d= d[1:]
	d = [x.split(',') for x in d]
	d = [ [string.strip(x[1],'\n"')] + [string.strip(x[4],'\n"')] for x in d]
	c_genes = dict(d)
	glist = loadtxt('./'+dname+'.ssm.txt',dtype='string')[1:,:]
	#inds = [i for i,a in enumerate(glist2) if a[1] in c_genes.keys()]
	#glist = [x[1] for x in glist2[inds]]
	glist = [x[1] for x in glist]
	glist = dict([(name,i) for i,name in enumerate(glist)])
	m = len(glist)
	S = zeros((m,m)) #similarity matrix (similar to pyclone)
	
	id_list = glist.keys()
	d = [i for i,m_id in enumerate(id_list) if c_genes[m_id] == "D"]
	b = [i for i,m_id in enumerate(id_list) if c_genes[m_id] == "B"]
	c = [i for i,m_id in enumerate(id_list) if c_genes[m_id] == "C"]
	a = [i for i,m_id in enumerate(id_list) if c_genes[m_id] == 'NA']
	ytrue = construct_matrix([d,c,b,a])
	np.savetxt('ytrue.txt',ytrue)
	flist = os.listdir(fdir)
	ns = len(flist) #number of MCMC samples
	for idx,fname in enumerate(flist):
		if fname.find('Store')>-1: continue
		f=open('./'+fdir+'/'+str(fname))
		tssb = cPickle.load(f)
		wts, nodes = tssb.get_mixture()
				
		for node in nodes:
			data = node.get_data()
			
			pids = [glist[datum.name.lower()] for datum in data if datum.name.lower() in glist.keys()];npids=len(pids)
			pids=sort(pids)
			#ids=[(pids[i],pids[j]) for i in arange(npids) for j in arange(npids)]
			ids= itertools.combinations(pids, 2)			
			for id in ids: 
				S[id[0],id[1]]+=1 # same cluster
				#S[id[1],id[0]]+=1
		f.close()
			
	
	S[arange(m),arange(m)]=ns

	f = open("vafs.csv")
	d = f.readlines()
	f.close()

	d= d[1:]
	d = [x.split(',') for x in d]
	d = [ [string.strip(x[1],'\n"')] + [string.strip(x[4],'\n"')] for x in d if x[4] != "NA\n"]
	c_genes = dict(d)
	glist2 = loadtxt('./'+dname+'.ssm.txt',dtype='string')[1:,:]
	inds = [i for i,a in enumerate(glist2) if a[1] in c_genes.keys()]
	np.savetxt('inds.txt', inds)

	ypred = S*1./ns
	ypred = ypred[inds,:][:,inds]
	

	ytrue = ytrue[inds,:][:,inds]
	np.savetxt('ytrue2.txt',ytrue)

	ids = triu_indices(len(inds))
	

	ytrue=ytrue[ids]	
	ypred=ypred[ids]	
	
	precision, recall, thresholds = mt.precision_recall_curve(ytrue,ypred)
	aucpr = mt.auc(recall, precision)

	fpr, tpr, thresholds = mt.roc_curve(ytrue, ypred)
	aucroc = mt.auc(fpr, tpr)
	
	return aucpr,aucroc	


if __name__ == "__main__":
	dname = sys.argv[1][:-1]
	fdir = dname+'/best/'
	aucpr,aucroc=compute_auc(fdir,dname)
	print dname,round(aucpr,4),round(aucroc,4)
	f = open(dname+'.aupr.txt','w')
	f.write('%s %f %f' % (dname,round(aucpr,4),round(aucroc,4)))
	f.close()

