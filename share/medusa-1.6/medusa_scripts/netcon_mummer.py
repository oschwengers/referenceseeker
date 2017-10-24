#from mummer_parser_ import *
from mummer_parser import *
import sys,os
#from IPython import embed
import networkx as nx
from cPickle import dump
from optparse import OptionParser,OptionGroup

######################################

if __name__=='__main__':
	usage=""" %prog [options]
	"""
	############ 'python netcon_mummer.py mapping_dir query_genome gexf_out [wscheme] [testing]'
	parser = OptionParser(usage=usage)
	# mandatory
	group1 = OptionGroup(parser, "Mandatory Arguments")
	group1.add_option("-i", "--input", dest="query_genome",
	help="target genome to be scaffolded", metavar="FILE")
	group1.add_option("-f", "--files", dest="mapping_dir",
	help="DIR where the comparison genomes are stored", metavar="DIR")
	parser.add_option_group(group1)
	group1.add_option("-o", "--output", dest="out",
	help="write graph to FILE", metavar="FILE")
	# optional
	group2 = OptionGroup(parser, "Optional Arguments")
	group2.add_option("-w", "--weightingscheme", dest="scheme", action="store_true",default=False,
	help="use a weighting scheme based on sequence similarity")
	group2.add_option("-t", "--testing", dest="testing", action="store_true",default=False,
	help="outputs a pkl object for testing [OLD]")
	group2.add_option("-d", "--distanceEstimation", dest="gap",default=0,
	help="choose how to compute gap length: 0 fixed distance (100bp) [default], 1 mean, 2 median, 3 most similar")
	parser.add_option_group(group2)
	(options, args) = parser.parse_args()
	if not options.query_genome or not options.mapping_dir or not options.out:
		parser.print_help()
		parser.error('Mandatory Arguments missing')
	query_genome,mapping_dir,out,scheme,testing,gap=options.query_genome,options.mapping_dir,options.out,options.scheme,options.testing,int(options.gap)
	if not mapping_dir.endswith('/'): mapping_dir+='/'
	
######################################

def sort_(clusters):
	edges=[]
	for cl in clusters.values():
		if len(cl) == 1: continue
		cl.sort(key=lambda x:int(x.rstart))
		for i in range(len(cl)-1):
			edges.append((cl[i],cl[i+1]))
	return edges

def update_edges_(G,Edge):
	""" old function """
	u,v,distance,orientation,weight,seqSimilarity=Edge.name1,Edge.name2,Edge.distance,Edge.orientation,Edge.weight,Edge.seqSimilarity
	w=G.get_edge_data(u,v,{'weight':0})['weight'] + weight
	d=G.get_edge_data(u,v,{'distance':0})['distance'] + distance
	o=G.get_edge_data(u,v,{'orientation':[]})['orientation']
	s=G.get_edge_data(u,v,{'seqSim':[]})['seqSim'] + [seqSimilarity]
	o.append(orientation)
	G.add_edge(u,v,weight=w,distance=d,orientation=o)
	
def update_edges(G,Edge):
	u,v,distance,orientation,weight,seqSimilarity=Edge.name1,Edge.name2,Edge.distance,Edge.orientation,Edge.weight,Edge.seqSimilarity
	if G.has_edge(u,v):
		G[u][v]['weight'] += weight
		G[u][v]['distance'] += [distance]
		G[u][v]['orientation'] += [orientation]
		G[u][v]['seqSim'] + [seqSimilarity]
	else:
		G.add_edge(u,v,weight=weight,distance=[distance],orientation=[orientation],seqSim=[seqSimilarity])
	
def compute_distances(G,method=0,outlier=1):
	""" Estimate the distance for each edge """
	import numpy as np
	skip=method==0
	for u,v in G.edges():
		if not skip:
			distances=np.array(G[u][v]['distance'])
			seqSims=np.array(G[u][v]['seqSim'])
			if method==1: dist=distanceEstimation_mean(distances,outlier)
			elif method==2: dist=distanceEstimation_median(distances,outlier)
			elif method==3: dist=distanceEstimation_MSH(distances,seqSims,outlier)
			else: dist=distanceEstimation_mean(distances,outlier)
		else: dist=0
		G[u][v]['distance']=dist
		G[u][v]['seqSim']=''

def distanceEstimation_mean(dist_list,method=1):
	""" Estimate the distance between two contigs, by computing the
		average of all distances found for these in the reference 
		genomes. A method for outlier detection is used in order to
		obtain reliable mean."""
	#print "using mean..."
	import numpy as np
	v=np.array(dist_list)
	if method == 1: mask=madBasedOutlier(v)
	else: mask=percentileBasedOutlier(v)
	distance=v[-mask].mean()
	return int(distance)

def distanceEstimation_median(dist_list,method=1):
	""" As the mother function, except that it uses the median instead
		of the mean."""
	#print "using median..."
	import numpy as np
	v=np.array(dist_list)
	return int(np.median(v))
	# try to remove outlier detection?
	#if method == 1: mask=madBasedOutlier(v)
	#else: mask=percentileBasedOutlier(v)
	#distance=np.median(v[-mask])
	#return distance

def distanceEstimation_MSH(dist_list,seqSim_list,method=1):
	""" As the mother function, except that it uses the most similar
		hit's distance."""
	#print "using most similar hit..."
	import numpy as np
	similarities,v=np.array(seqSim_list),np.array(dist_list)
	# should remove outlier detection here
	if method == 1: mask=madBasedOutlier(v)
	else: mask=percentileBasedOutlier(v)
	distance=v[np.where(max(similarities))][0]
	return int(distance)

def madBasedOutlier(points, thresh=3.5):
	""" return outliers based on median-absolute-deviation (MAD) 
		This performs better on small data samples (<20).
		Require numpy object."""
	import numpy as np
	np.seterr(divide='ignore', invalid='ignore')
	#from IPython import embed
	if len(points.shape) == 1: points = points[:,None]
	median = np.median(points, axis=0)
	diff = np.sum((points - median)**2, axis=-1)
	diff = np.sqrt(diff)
	med_abs_deviation = np.median(diff)
	modified_z_score = 0.6745 * diff / med_abs_deviation
	#if med_abs_deviation == 0: embed()
	return modified_z_score > thresh
	
def percentileBasedOutlier(data, threshold=95):
	""" return outliers based on percentiles 
		This performs better on big data samples (>20).
		Require numpy object"""
	import numpy as np
	diff = (100 - threshold) / 2.0
	minval, maxval = np.percentile(data, [diff, 100 - diff])
	return (data < minval) | (data > maxval)


def initialize_graph(genome):
	import networkx as nx
	from Bio.SeqIO import parse
	G=nx.Graph()
	contigs=parse(genome,'fasta')
	for c in contigs:
		id_,length=c.id,len(c.seq)
		G.add_node(id_,length=length)
	return G

def adjust_orientations(G):
	id_=0
	for e in G.edges():
		n1,n2=e
		G[n1][n2]['orientation']=convert_orientations(e,G[n1][n2]['orientation'])
		max_count=G[n1][n2]['orientation'].count(max(G[n1][n2]['orientation'],
					     key=lambda x:G[n1][n2]['orientation'].count(x)))
		#embed()
		G[n1][n2]['orientation_max']=list({tuple(i) for i in G[n1][n2]['orientation'] if G[n1][n2]['orientation'].count(i)==max_count})
		G[n1][n2]['orientation_max']='==='.join(['=='.join(i) for i in G[n1][n2]['orientation_max']])
		l=G[n1][n2]['orientation']
		counts={'_'.join(i):l.count(i)/float(len(l)) for i in l}
		#G[n1][n2]['orientation']='__'.join(['%s&%s' %(k,v) for k,v in counts.items()])
		G[n1][n2]['orientation']=''
		G[n1][n2]['id']=id_
		id_+=1
	#edges=sorted(G.edges(), key=lambda x: G[x[0]][x[1]]['orientation_max'])
	#for e in edges:
		#n1,n2=e[:2]
		#G[n1][n2]['id']=id_
		#id_+=1
		#embed()
	
def format_orientation_string(hit1,hit2):
	orientations=['%s:%s' %(hit1.name,hit1.orientation), '%s:%s' %(hit2.name,hit2.orientation)]
	#orientations.sort()
	return orientations

def convert_orientations(e,ori_list):
	''' convert elements of a list i.e. a:1,b:-1 to b:1,a:-1 '''
	#embed()
	n1_,n2_=e
	ori_new=[]
	for l in ori_list:
		if l[0].split(':')[0] != n1_:
			#print 'inverting!',l[0].split(':')[0],n1_, len(ori_list) 
			n1,v1,n2,v2=[i for j in l for i in j.split(':')]
			v1,v2=int(v1)*-1,int(v2)*-1
			l_=['%s:%s' %(n2,v2),'%s:%s' %(n1,v1)]
		else: l_=l
		ori_new.append(l_)
	#embed()	
	return ori_new


######################################

class Edge(object):
	def __init__(self,hit1,hit2,wscheme=0):
		self.name1,self.name2=hit1.query,hit2.query
		self.distance=hit1.distance_from(hit2)
		self.seqSimilarity=hit1.weight2+hit2.weight2
		doMapWithin(hit1,hit2)
		self.orientation=format_orientation_string(hit1,hit2)
		if wscheme==0:self.weight=1
		else: self.weight=self.seqSimilarity

#####################################

def testForPrune():
	""" function to profile the script """
	mapping_dir='./'
	out='prova.gexf'
	query_genome='test/Rhodobacter_target.fna'
	inputs=[f for f in os.listdir(mapping_dir) if f.endswith('.coords')]
	G=initialize_graph(query_genome)
	for coord in inputs:
		print('using input',coord)
		clusters=parse_mummer(mapping_dir + coord)
		edges=sort_(clusters)
		for e in edges:
			if not testing: update_edges(G,Edge(*e,wscheme=scheme))
			else: update_edges(G,Edge(*e,wscheme=scheme))
	#embed()
	print('adjusting orientations')
	adjust_orientations(G)
	compute_distances(G,method=gap)
	if not testing: nx.write_gexf(G,out)
	else: dump(G,open(out,'w'))

######################################
import sys
if __name__ == '__main__':

	inputs=[f for f in os.listdir(mapping_dir) if f.endswith('.coords')]
	G=initialize_graph(query_genome)
	for coord in inputs:
		print('using input',coord)
		clusters=parse_mummer(mapping_dir + coord)
		edges=sort_(clusters)
		for e in edges:
			if not testing: update_edges(G,Edge(*e,wscheme=scheme))
			else: update_edges(G,Edge(*e,wscheme=scheme))
		#sys.exit()
	#embed()
	print('adjusting orientations')
	adjust_orientations(G)
	compute_distances(G,method=gap)
	if not testing: nx.write_gexf(G,out)
	else: dump(G,open(out,'w'))
