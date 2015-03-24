
###### code to sample from the paramater posterior p(\phi | data) ########

import numpy
from numpy import *
from data import Datum

from tssb import *

from util import dirichletpdfln
from numpy.random import dirichlet

def metropolis(tssb,iters=1000,std=0.01,burnin=0,n_ss=0,n_cnvs=0):

	wts, nodes = tssb.get_mixture()
	pi = empty(len(wts))

	sample_cons_params(tssb)
	update_params(tssb)
	pi = getpi(tssb)
		
	ctr=0
	pi_new  = empty(pi.shape)	
	for i in arange(-burnin,iters):
		cf=0
		pi_new = sample_cons_params1(tssb,pi,std)
		a=param_posterior(tssb,1)-param_posterior(tssb,0) + correction_term(pi,pi_new,std) - correction_term(pi_new,pi,std)
		if log(rand(1)) < a:
			update_params(tssb)
			pi = pi_new
			ctr+=1		
	return ctr*1./iters


# tree-structured finite-dimensional stick breaking
def sample_cons_params(tssb):
	def descend(root):
	
		if root.parent() is None:
			root.params1 = 1
			root.pi1 = root.params1*rand(1)[0] # break nu stick
		r = root.params1-root.pi1 #mass assigned to children
		p = rand(len(root.children()));p=r*p*1./sum(p)
		index=0
		for child in root.children():			
			#print child.params1
			child.params1 = p[index]# break psi sticks			
			child.pi1 = child.params1*(rand(1)[0]**(len(child.children())>0)) # break nu stick
			index+=1
		for child in root.children():			
			descend(child)	

	descend(tssb.root['node'])
	
			
# no stick breaking, updates pi with small perturbations
def sample_cons_params1(tssb,pi,mh_std=0.01):	
	std = mh_std
	pi = dirichlet(std*pi+1) # for dirichlet proposal
	pi = pi+0.0001; pi=pi/sum(pi)
	pi=list(pi);pi.reverse()
		
	def descend(root):		
		for child in root.children():			
			descend(child)
		root.pi1=pi.pop();
		if len(root.children())==0:
			root.params1 = root.pi1
		else:			
			root.params1 = root.pi1+sum([child.params1 for child in root.children()])			
	descend(tssb.root['node'])
	return getpi(tssb,1)

	
# MH correction terms for asymmetric proposal distribution
def correction_term(pi1,pi2,std):
	return dirichletpdfln(pi1,std*pi2) # for dirichlet proposal

def param_posterior(tssb, new):
	def getposterior(root,new):
		llh = 0
		data = root.get_data()
		if new: 
			p = root.params1
		else:
			p = root.params
		llh = sum([data[i].__log_likelihood__(p,update_tree=False,new_state=new) for i in arange(len(data))])	

		for child in root.children():					
			llh += getposterior(child,new)
		return llh
	return getposterior(tssb.root['node'],new)	

def update_params(tssb):
	def descend(root):			
		for child in root.children():
			descend(child)	
		root.params = root.params1
		root.pi = root.pi1
	descend(tssb.root['node'])

def getpi(tssb,new=0):
	pi = []
	def descend(root):
		for child in root.children():		
			descend(child)
		if new: 
			p = root.pi1
		else:
			p = root.pi
		pi.append(p)
	descend(tssb.root['node'])
	pi=array(pi);pi.shape=len(pi),	
	return pi