
#######################

if __name__=='__main__':
	usage="""
	python mummer_parser.py coord_file out_file """
	args=sys.argv
	coords=args[1]
	out=args[2]

#######################

def parse(coords_file):
	forbid={0,1,2,3,4}
	for i,l in enumerate(open(coords_file)):
		if i in forbid: continue
		yield mummer_hit(l.strip())

def get_bestHits(hits,attr='covq'):
	query_contigs=set([h.query for h in hits])
	best_hits=[]
	for c in query_contigs:
		best_hit=max([h for h in hits if h.query == c],key=lambda x: float(x.__getattribute__(attr)))
		best_hits.append(best_hit)
	return best_hits

#########################
# start of experimental
#########################

def compare(x,y,attr):
	return float(x.__getattribute__(attr)) > float(y.__getattribute__(attr))


def getBestHits(coords,attr='covq'):
	""" experimental version of get_bestHits with better performances.
		Improvements:
			- use of generators
			- use of a dictionary for best hits:
				. key = hit.query
				. value = hit
			- use of a dictionary for hit clusters:
				. key = hit.reference
				. value = [hit1,hit2,...,hitn]
				. hits are added in order to compare the
				"""
	best_hits={}
	for h in parse(coords):
		if h.query not in best_hits:
			best_hits[h.query]=h
			continue
		if compare(h,best_hits[h.query],attr): best_hits[h.query]=h
	return best_hits.values()

#########################
# end of experimental
#########################


def get_Clusters(best_hits):
	clusters={}
	for h in best_hits:
		clusters[h.reference]=clusters.get(h.reference,[])
		clusters[h.reference].append(h)
	return clusters

def write_Clusters(clusters,out):
	out=open(out,'w')
	for cl in clusters.values():
		cl.sort(key=lambda x:int(x.rstart))
		out.write('\n'.join([c.query for c in cl])+'\n\n')

def parse_mummer(coords):
	hits=parse(coords)
	best_hits=getBestHits(coords)
	clusters=get_Clusters(best_hits)
	return clusters

def parse_mummer2(coords):
	hits=parse(coords)
	best_hits=get_bestHits2(hits)
	clusters=get_Clusters(best_hits)
	return clusters

def do_overlap(a,b):
	'''return true if a (a1,a2) overlap with b (b1,b2)'''
	sol=	((max(a) > min(b)) and (max(a) < max(b))) or \
		((min(a) > min(b)) and (min(a) < max(b)))
	return sol

def doMapWithin(hit1,hit2):
	'''return true if a (a1,a2) maps within b (b1,b2) or viceversa'''
	a,b=[int(hit1.rstart),int(hit1.rend)],[int(hit2.rstart),int(hit2.rend)]
	if ((max(a) > max(b)) and (min(a) < min(b))): print('%s maps within %s !!!' %(hit2.name, hit1.name))
	elif ((max(b) > max(a)) and (min(b) < min(a))): print('%s maps within %s !!!' %(hit1.name, hit2.name))
	return


#######################


class mummer_hit(object):
	def __init__(self,line):
		self.qstart,self.qend,self.rstart,self.rend,self.len1,self.len2,self.percidy,\
		self.lenr,self.lenq,self.covq,self.covr,self.query,self.reference=[i for l in line.split(' | ') for i in l.split()]
		self.name=self.query
		self.weight1=float(self.percidy)
		self.weight2=float(self.percidy)*float(self.covq)
		self.weightNaif=1
		if int(self.rstart)>int(self.rend): self.orientation=-1
		else: self.orientation=1
	def distance_from(self,hit):
		a1,a2,b1,b2=int(self.rstart),int(self.rend),int(hit.rstart),int(hit.rend)
		if do_overlap([a1,a2],[b1,b2]): return 0
		distances=[a1-b1,a1-b2,a2-b1,a2-b2]
		distance=min([abs(i) for i in distances])
		return distance

#######################

if __name__=='__main__':
	clusters=parse_mummer(coords)
