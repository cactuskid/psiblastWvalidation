
from sklearn.cluster import MiniBatchKMeans as MBKM
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
from csb.bio.io.hhpred import HHOutputParser
from colour import Color
import matplotlib.patches as mpatches

import pandas as pd
import multiprocessing as MP
import numpy as np
from csb.bio.io.hhpred import HHOutputParser
import networkx as nx
import math
import glob

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering

def metaclust_filter(ID):
	if 'Representative=' in ID:
		clean = ID.split('Representative=')[-1]

		return clean
	else:
		return None
def return_scaffolds(Representative):
	if Representative is not None:
		if 'scaffold' in Representative:
			return Representative.split('scaffold')[-1]
	else:
		return None

def return_UPI(Representative):
	if Representative is not None:
		if 'UPI' in Representative:
			return Representative.split('Representative=')[-1]
	else:
		return None

def return_gi(Representative):
	if Representative is not None:
		if 'gi' in Representative:
			return Representative.split('|')[1]
	else:
		return None




def cluster_distmat(DM, protlabels, clustersize = 20 , visualize= False , cluster=None):
	print(len(protlabels) / clustersize)
	print(len(protlabels))
	#labels = MBKM(n_clusters= int(len(protlabels) / clustersize ), verbose = True ).fit_predict(DM)
	DM = StandardScaler().fit_transform(DM)
	np.random.seed(0)
	labels = MBKM(n_clusters= int(len(protlabels) / clustersize ) ).fit_predict(DM)


	#labels = DBSCAN( metric = 'precomputed' , n_jobs=-1).fit_predict(DM )

	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	print('Estimated number of clusters: %d' % n_clusters_)
	clusters = {}

	reverse = {}
	if visualize == True:
		plt.imshow(DM)
		plt.show()
	if cluster is not None:
			appendstr = cluster+ '.'
	else:
		appendstr = ''
	for i, label in enumerate(labels):
		if appendstr + str(label) not in clusters:
			clusters[appendstr + str(label)] = []
		clusters[appendstr + str(label) ].append(protlabels[i])
		reverse[protlabels[i]] = appendstr + str(label)
	index = [ protlabels.index(prot)  for prot in reverse]
	return clusters , reverse , labels ,index

def itercluster(DM, protlabels, clustersize = 20 ):
	clusters , reverse , labels, index = cluster_distmat(DM, protlabels , clustersize )
	size = [ len( clusters[cluster]) for cluster in clusters ]
	iterations = 1

	print(clusters.keys())
	print(size)
	print(sum(size))
	while iterations < 4:
		updates = []
		remove = []
		for cluster in clusters:
			if len(clusters[cluster]) > clustersize:
				print('breaking up '+cluster)
				index = [ protlabels.index(prot) for prot in clusters[cluster] ]
				subDM = DM[index , :]
				subDM = subDM[:,index]
				newclusters , reverse , labels ,index = cluster_distmat(subDM, clusters[cluster] , clustersize/2, cluster= cluster  )
				updates.append(newclusters)
				remove.append(cluster)
		for update in updates:
			clusters.update(update)
		for rem in remove:
			del clusters[rem]
		size = [ len(clusters[cluster]) for cluster in clusters ]
		print(size)
		print(sum(size))
		iterations +=1


	reverse ={}
	for cluster in clusters:
		for prot in clusters[cluster]:
			if prot not in reverse:
				reverse[prot]= cluster
			if len(cluster)> len(reverse[prot]):
				reverse[prot]=cluster
	index = [ protlabels.index(prot)  for prot in reverse]
	labels = reverse.values()

	return clusters , reverse, labels, index



def ClusterPlot(DM,color_labels):
	#DM = StandardScaler().fit_transform(DM)
	mapper=MDS(n_components=3, metric=True, n_init=4, max_iter=300, verbose=1, eps=0.001, n_jobs=-1, random_state=0, dissimilarity='precomputed')
	X =mapper.fit_transform(DM)
	unique_labels = list(set(color_labels))
	c1 = Color('red')
	c2 = Color('blue')
	colorvec = list(c1.range_to(c2, len(unique_labels)))
	colordict = { unique_labels[i] : colorvec[i] for i in range(len(colorvec))}

	fig = plt.figure()
	ax = fig.add_subplot(111,  projection='3d')
	ptcolors = [  colordict[k].rgb for k in color_labels ]
	ax.scatter(X[:, 0], X[:, 1] , X[:, 2] ,  color=ptcolors  )

	handles=[]
	for label in unique_labels:
		patch = mpatches.Patch(color=colordict[label].rgb, label=label)
		handles.append(patch)

	plt.legend(handles=handles)
	plt.show()
	return colordict

	"""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	mapper=MDS(n_components=2, metric=True, n_init=4, max_iter=300, verbose=1, eps=0.001, n_jobs=-1, random_state=0, dissimilarity='precomputed')
	X =mapper.fit_transform(DM)
	for filename in unique_files:
		index = np.where( filename == np.asarray(files) )[0]
		marker = (3+unique_files.index(filename), 0, 0 )
		ptcolors = [  colorvec[unique_labels.index(k)].rgb for k in labels[index] ]
		ax.scatter(X[index, 0], X[index, 1] , marker =marker, color=ptcolors  )
	handles=[]
	for label in unique_labels:
		patch = mpatches.Patch(color=colorvec[unique_labels.index(label)].rgb, label=label, alpha = .2)
		handles.append(patch)
	plt.legend(handles=handles)
	plt.show(	)"""


def DMPlot(DM,labels,files):
	#use a different marker for each input file
	#use a different color for each cluster_distmat
	pass

def blastscore(args):
	y = args
	return 1/( 1-math.log(y))

def return_index(args):
	x,y,z = args
	try:
		return (z.index(x), z.index(y))
	except:
		pass

def Blast_parseTo_DM(blastTab,  outdir , E_VALUE_THRESH=.01):
	#make this a sparse mat above a certain number of seqs?
	print('reading blast allvall output, creating distmat and index')
	scores = {}
	df = blastToDF(blastTab)
	df = df[df.evalue < E_VALUE_THRESH]
	subdf = df[df.evalue >= 10**-200 ]
	identitical = df[ df.evalue <= 10**-200 ]
	pool = MP.Pool(30)
	subdf['score'] = pool.map_async( blastscore, subdf['evalue'] , 100).get()
	subdf['score'] = subdf['score']/subdf['length']

	identitical['score'] = 0
	df = pd.concat( [subdf , identitical])

	labels = df.qseqid.unique().tolist()

	tuples = zip( df['sseqid']  , df['qseqid'] , [labels]*df['qseqid'].size )

	df['Coordtuple'] = pool.map_async(return_index , tuples).get()
	df = df[pd.notnull( df['Coordtuple']) ]

	print('blast scores done calculating')
	print(df)

	DM = np.ones((len(labels),len(labels)))*1000
	DM[tuple(zip(*df['Coordtuple']) )] = df['score']
	#symmetrical even though blast usually isn't...
	DM += DM.T
	pool.close()
	print(DM)
	print('DONE')
	return DM, df , labels

def blastToDF(blastTab , columnsstr = 'qseqid qlen slen qcovhsp sseqid staxids bitscore score evalue pident qstart qend sstart send'):
	
	columns = columnsstr.split()
	df = pd.read_csv(blastTab , names = columns , header= None)
	df['evalue']=df['evalue'].convert_objects(convert_numeric=True)
	return df

def distmat_to_txt( labels , distmat, outpath):
	print(labels)
	print(distmat.shape)
	outstr = str(len(labels)) + '\n'

	for i,pdb in enumerate(labels):
		namestr = pdb[0:20]
		outstr += namestr+ ' ' + np.array2string( distmat[i,:], formatter={'float_kind':lambda x: "%.2f" % x} , precision = 8 ).replace('[', '').replace(']', '').replace('\n', '' )  + '\n'

	handle = open(outpath, 'w')
	handle.write(outstr)
	handle.close()
	return  outpath , outstr

def HHSearch_parseTo_DMandNX(hhrs , labels=None ):
	clusternames = []
	for i,hhr in enumerate(hhrs):
		try:
			profile = HHOutputParser(alignments=False).parse_file(hhr)
			if profile.query_name not in clusternames or labels != None:
				if labels == None:
					clusternames.append(profile.query_name)
				else :
					clusternames.append(labels[i])
		except:
			print(hhr)
			pass

	print(clusternames)
	evalDM = np.ones( (len(clusternames),len(clusternames) ))
	pvalDM = np.ones( (len(clusternames),len(clusternames) ))
	scoreDM = np.zeros( (len(clusternames),len(clusternames) ))
	SSDM = np.zeros( (len(clusternames),len(clusternames) ))
	probaDM = np.zeros( (len(clusternames),len(clusternames) ))
	lenDM =  np.ones( (len(clusternames),len(clusternames) ))
	NX = nx.Graph()
	for i,hhr in enumerate(hhrs):
		protlist = []
		profile = HHOutputParser(alignments=False).parse_file(hhr)
		for hit in profile:
			DMscore = float(hit.evalue)
			proba = hit.probability

			if 'anchor' not in hit.id and 'anchor' not in profile.query_name:
				i = clusternames.index(hit.id.strip() )
				j = clusternames.index(profile.query_name.strip())

				if hit.evalue < evalDM[i,j]:
					evalDM[i,j] = hit.evalue
					evalDM[j,i] = evalDM[i,j]

				if hit.pvalue < pvalDM[i,j]:
					pvalDM[i,j] = hit.pvalue
					pvalDM[j,i] = pvalDM[i,j]

				if scoreDM[i,j] < hit.score:
					scoreDM[i,j] = hit.score
					scoreDM[j,i] = scoreDM[i,j]

				if SSDM[i,j] < hit.ss_score:
					SSDM[i,j] = hit.ss_score
					SSDM[j,i] = SSDM[i,j]


				if probaDM[i,j] < hit.probability:
					probaDM[i,j] = hit.probability
					probaDM[j,i] = probaDM[i,j]

				#use smallest of the two prots
				if lenDM[i,j] == 1 or lenDM[i,j] > hit.qlength:
					lenDM[i,j] = hit.qlength
					lenDM[j,i] = lenDM[i,j]

			if hit.id != profile.query_name :
				NX.add_edge( hit.id , profile.query_name )
				NX[hit.id][profile.query_name]['score']= hit.score
	return probaDM, evalDM ,pvalDM,  lenDM , scoreDM, SSDM, NX , clusternames
