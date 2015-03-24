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

def compute_auc(fdir,dname):

	c,r,n = [int(x) for x in dname.split('.')[1:4]]
	
	ytrue = get_true(dname)	
	
	m = n*(c) # number of genes / instances
	glist2 = loadtxt('./'+dname,dtype='string')[1:,:]
	inds = [i for i,a in enumerate(glist2)]
	glist = [x[0] for x in glist2[inds]]
	glist = dict([(name,i) for i,name in enumerate(glist)])
	yt2 = array(ytrue)
	yt2 = reshape(yt2,(m,m))
	yt2 = yt2[inds,:][:,inds]
	ytrue = yt2
	m = len(glist)
	S = zeros((m,m)) #similarity matrix (similar to pyclone)
	flist = os.listdir(fdir)
	ns = len(flist) #number of MCMC samples
	
	for idx,fname in enumerate(flist):
		if fname.find('Store')>-1: continue
		f=open('./'+fdir+'/'+str(fname))
		tssb = cPickle.load(f)
		f.close()
	
		#cluster info
		wts, nodes = tssb.get_mixture()
				
		for node in nodes:
			data = node.get_data()
			
			pids = [glist[datum.name.lower()] for datum in data if datum.name.lower() in glist];npids=len(pids)
			pids=sort(pids)
			#ids=[(pids[i],pids[j]) for i in arange(npids) for j in arange(npids)]
			ids= itertools.combinations(pids, 2)			
			for id in ids: 
				S[id[0],id[1]]+=1 # same cluster
				#S[id[1],id[0]]+=1
	
	S[arange(m),arange(m)]=ns
	ypred = S*1./ns
	ids = triu_indices(m)
	
	ytrue=ytrue[ids]	
	ypred=ypred[ids]	
	
	precision, recall, thresholds = mt.precision_recall_curve(ytrue,ypred)
	aucpr = mt.auc(recall, precision)


	fpr, tpr, thresholds = mt.roc_curve(ytrue, ypred)
	aucroc = mt.auc(fpr, tpr)
	
	return aucpr,aucroc	


def get_true(fname):
	c,n = [int(x) for x in fname.split('.')[1:4:2]]
	c = c
	clustering = []
	for i in range(c):clustering.append([x+(i*n) for x in range(n)])
	return construct_matrix(clustering)


def test():

	nc=[3,4,5,6] # number of cluster
	rd=[20,30,50,70,100,200,300] # read depth
	nssms = [5,10,25,50,100,200,500,1000,2000]
	
	res = []
	ctr=0
	for i in nc:
		for j in rd:
			for k in nssms:			
				fdir='./emptysim.'+str(i)+'.'+str(j)+'.'+str(k)+'.ssm.txt.output/best/'
				dname='emptysim.'+str(i)+'.'+str(j)+'.'+str(k)+'.ssm.txt'
				aucpr,aucroc, = compute_auc(fdir,dname)
				print i,j,k,aucpr,aucroc
				ctr+=1
				res.append([int(i),int(j),int(k),aucpr,aucroc])	
				

def post_process():
	dir = './out/pbs-output/'
	flist = os.listdir(dir)
	res = []
	for fname in flist:	
		if fname.find('DS_Store')>=0: continue;
		#if fname.find('python')<0: continue;
		f= open(dir+fname)
		lili = [line.split() for line in f.readlines()]
		f.close()
		if lili[2][0].isdigit():
			res.append(lili[2][0::])
		else:
			print lili
	res = array(res)
	savetxt('result',res,fmt='%s')
	
if __name__ == "__main__":
	if 1:
		c = sys.argv[1]
		r = sys.argv[2]
		n = sys.argv[3]
		fdir='lin2.'+str(c)+'.'+str(r)+'.'+str(n)+'.ssm.txt.output/best/'
		dname='lin2.'+str(c)+'.'+str(r)+'.'+str(n)+'.ssm.txt'
		aucpr,aucroc=compute_auc(fdir,dname)
		print c,r,n,round(aucpr,4),round(aucroc,4)
		f = open('lin.'+c+'.'+r+'.'+n+'.aupr.txt','w')
		f.write('%s %s %s %f %f' % (c,r,n,round(aucpr,4),round(aucroc,4)))
		f.close()
