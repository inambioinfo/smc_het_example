import sys
from numpy import *
from scipy.spatial import distance
from sklearn.metrics import precision_recall_curve,auc,roc_curve

#Construct true co-clustering matrix

def get_true(fname):
	c,n = [int(x) for x in fname.split('.')[1:4:2]]
	c = c -1
	clustering = []
	for i in range(c):
		clustering.append([x+(i*n) for x in range(n)])
	return construct_matrix(clustering)

def get_true_direct(c,n):
	clustering = []
	for i in range(c):
		clustering.append([x+(i*n) for x in range(n)])
	return construct_matrix(clustering)

def get_obs(fname):
	f = open(fname)
	d = f.read()
	d = d.split('\n\n')
	d = d[0].split('\n')
	d = d[1:]
	d = [x.split(',\t') for x in d]
	d = [x[4] for x in d]
	d = [x.split('; ') for x in d]
	d = [[int(y[1:]) for y in x if y] for x in d]
	return construct_matrix(d)
	
def construct_matrix(clustering):
	n = sum([len(x) for x in clustering])
	m = zeros((n,n))
	cpoints = [get_cluster(clustering,x) for x in range(n)]
	for i in range(n):
		for j in range(n):
			if cpoints[i] == cpoints[j]:
				m[i,j] = 1
	return m

def get_cluster(clustering,point):
	return [point in x for x in clustering].index(True)

def score(m1,m2):
	savetxt("true.txt",m1,delimiter=",")
	savetxt("obs.txt",m2,delimiter=",")
	precision,recall,thresholds =  precision_recall_curve(m1.flatten(),m2.flatten())
	pr = auc(recall,precision)
	fpr,tpr,thresholds = roc_curve(m1.flatten(),m2.flatten())
	roc = auc(fpr,tpr)
	cosd = distance.cosine(m1.flatten(),m2.flatten())
	return pr,roc,cosd

if __name__ == '__main__':
	m1 = get_true(sys.argv[1])
	m2 = get_obs(sys.argv[1])
	c,r,n = [int(x) for x in sys.argv[1].split('.')[1:4]]
	pr,roc,cosd = score(m1,m2)
	print c,r,n,pr,roc,cosd
