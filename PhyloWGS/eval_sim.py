import sys
import math

def eval_accuracy(tkt_file, true_n, min_genes=0):
	f = open(tkt_file)
	d = f.read()
	d = d.split('\n\n')
	d = d[0].split('\n')
	d = d[1:]
	d = [x.split(',\t') for x in d]
	# Remove small clusters
	d = [x for x in d if int(x[3]) > min_genes]
	#Error in number of clusters
	err_clusters = len(d) - true_n
	return err_clusters

def get_rd(fname):
	return int(fname.split('.')[2])

def get_nssms(fname):
	return int(fname.split('.')[3])

if __name__ == '__main__':
	t_num = int(sys.argv[1])
	rds = [int(x) for x in sys.argv[2].split(',')]
	nssms = [int(x) for x in sys.argv[3].split(',')]
	adapt_min_genes = sys.argv[4]=="adapt"
	files = sys.argv[5:]
	res = []
	try:
		for i in range(len(rds)):
			res.append([rds[i]])
			for j in range(len(nssms)):
				if adapt_min_genes:
					min_genes = math.floor(nssms[j]/10)+1
				else:
					min_genes = 2
				matching_f = [x for x in files if get_rd(x) == rds[i] and get_nssms(x) == nssms[j]][0]
				res[i].append(eval_accuracy(matching_f,t_num,min_genes=min_genes))
				print t_num,rds[i],nssms[j],eval_accuracy(matching_f,t_num,min_genes=min_genes)
	except IndexError:
		pass
		#print rds[i],nssms[j]
	#print '\n'.join(['\t'.join([str(y) for y in x]) for x in [nssms]+res])
