
import csv
import pandas as pd
from Bio import SeqIO
from tempfile import NamedTemporaryFile as namedTemp
import time
import pickle
import glob
import functions
import blastcluster
import multiprocessing as mp
import numpy as np
import config
from Bio import Entrez
from matplotlib import pyplot as plt
from random import shuffle

Entrez.email = "dmoi@unil.ch"


fastadict={}
taxadict ={}
genbankdict={}
dfdict ={}



files = glob.glob(config.inputdir +'*.fasta')
for file in files:
	for i, seq_record in enumerate(SeqIO.parse(file, 'fasta')):
		dfdict[seq_record.id] = { 'file':file , 'full':str(seq_record) , 'seq': str(seq_record.seq) , 'id': seq_record.id , 'fasta': '>'+seq_record.id +'\n' + str(seq_record.seq) +'\n'}
	else:
		print(i)

#make seqDF from all results
seqDF = pd.DataFrame.from_dict(dfdict , orient='index')
fastaout = ''.join(seqDF['fasta'].tolist())
handle= open(config.blastdir + 'ALLresults.fasta', 'w')
handle.write(fastaout)
handle.close()

print(seqDF)

print( 'sequences loaded')

if config.blastall == True:

	print( 'blast all v all')

	#blast all v all
	p1, allvall = functions.runBlastALL( config.fastadir+'ALLresults.fasta','blastp', 'makeblastdb', './blastdir/' , mp.cpu_count() )

if config.distmat == True:
	print( 'make blast DM')

	DM, df , protlabels = blastcluster.Blast_parseTo_DM(config.fastadir + 'allVall.tab' ,  './' , E_VALUE_THRESH=.01)
	print(len(protlabels))
	print(len(df.qseqid.unique()))
	#dump blast distance calc
	with open('BlastallOBjects.pkl' , 'wb')as matplickle:
		pickle.dump([DM ,df, protlabels],matplickle,-1)

if config.clusterDistmat == True:
	print( 'clustering with blast DM')

	with open('BlastallOBjects.pkl' , 'rb')as matplickle:
		DM ,df, protlabels = pickle.load(matplickle)
	clusters ,reverse, labels, index = blastcluster.cluster_distmat(DM, protlabels, clustersize = 200)
	seqDF['cluster'] = seqDF['id'].map(reverse)
	#df['compositionScore'] = df['seq'].map(compositionScore)
	#seqDF['compositionScore'] = seqDF['id'].map(lambda x : if 'Representative' in x x.split('Representative=')[-1] else None)
	seqDF['Representative']= seqDF['id'].map(blastcluster.metaclust_filter)
	seqDF['scaffold']= seqDF['Representative'].map(blastcluster.return_scaffolds)
	seqDF['gi']= seqDF['Representative'].map(blastcluster.return_gi)
	seqDF['upi']= seqDF['Representative'].map(blastcluster.return_UPI)
	UPIS = seqDF.upi.unique()
	Scaffolds = list(seqDF.scaffold.unique())[1:]
	print(Scaffolds)
	filenames = [ seqDF['file'][id] for id in labels ]
	with open('clusterOBjects.pkl' , 'wb')as matplickle:
		pickle.dump([clusters ,reverse, protlabels , filenames, index, DM, seqDF ],matplickle,-1)
	seqDF.to_csv('clusters.csv')


if config.makemodels == True:
#make alignments and models
	print( 'aln clusters')

	with open( 'clusterOBjects.pkl', 'rb') as clusterobjects:
	    c = pickle.load(clusterobjects)
	clusters ,reverse, protlabels, files, index, subDM, seqDF   = c
	clusters = seqDF['cluster'].unique()
	for clust in clusters:
		print(len(seqDF['fasta'][seqDF.cluster == clust]))
		#up to 100 random seqs to decrease aln time
		fasta = seqDF['fasta'][seqDF.cluster == clust].tolist()
		shuffle(fasta)
		#subsample big fasta
		if len(fasta) > 300:
			outstr =''.join(fasta[0:300])
		else:
			outstr =''.join(fasta)
		fastafile = config.alndir + str(clust) + '.fasta'
		with open( fastafile , 'w') as fastaout:
			fastaout.write(outstr)
		alnfile = config.alndir + str(clust) + 'aln.fasta'
		with open( alnfile , 'w') as alnout:
			alnout.write(functions.runclustalo(fastafile, verbose= True).replace('*', '-') )

if config.HMMall == True:
	print( 'prepare hmms for each cluster')

	hmms = []
	models = []
	alns = glob.glob(alndir + '*aln.fasta')
	print(alns)
	with open( 'clusterOBjects.pkl', 'rb') as clusterobjects:
		c = pickle.load(clusterobjects)
	clusters ,reverse, protlabels, files, index, subDM, seqDF   = c
	okclusters=[]
	for file in config.seeds:
		okclusters += list(seqDF[seqDF.file == file ].cluster.unique())
	print(okclusters)
	for aln in alns:
		#output = functions.runreformat(aln, aln.split('.fasta')[0] + '.a3m', True)
		#output = functions.runaddSS(aln.split('.fasta')[0] + '.a3m', aln.split('.fasta')[0] + '.a3m' , True)
		output = functions.runHHmake(aln,aln.split('.fasta')[0] + '.hhm' ,verbose=True)
		hmms.append(aln.split('.fasta')[0] + '.hhm')

	#rename HMMs to cluster name
	hmms = glob.glob( alndir + '*.hhm')
	for hmm in hmms:
		print(hmm)
		newstr= ''
		with open(hmm, 'r') as instr:
			for line in instr:
				if 'NAME' in line:
					#q denotes query or seed sequences are present
					for cluster in okclusters:
						if str(cluster)+'aln.hhm' == hmm.split('/')[-1]:
							print(line)
							newstr += 'NAME ' + hmm.split('/')[2]+'q'+'\n'
							break
						else:
							newstr += 'NAME ' + hmm.split('/')[2]+'\n'
							break
				else:
					newstr += line

		with open(hmm, 'w') as outstr:
			outstr.write(newstr)
	output = functions.runHHDB(alndir, 'clusters', True)
	results=[]
	## TODO: add precomputed hmms represent gene clustersize
	for model in hmms:
		print(model)
		outfile= model.replace('hhm', 'hhr')
		output = functions.runHHSearch(model , outfile , './clusters_hhm_db' , verbose = True)
		results.append(outfile)

if config.HHDM_compile == True:
	print('compile all v all hmm results')
	results = glob.glob( './aln/*.hhr')
	probaDM, evalDM ,pvalDM,  lenDM , scoreDM, SSDM, NX , clusternames = blastcluster.HHSearch_parseTo_DMandNX(results , None )
	HHDM =  1 / ( (scoreDM +  3*SSDM) )
	np.fill_diagonal(HHDM , 0 )
	HHDM[HHDM == np.inf ] = 10**5
	HHDM /= 1000
	matfile = config.datadir+'/HHDM.mat'
	outpath,outstr = blastcluster.distmat_to_txt(clusternames , HHDM, matfile )
	outfile = matfile+'_tree.txt'
	output = functions.runFastme(matfile, outfile)
	print(output)

if config.PDB70Validate == True:
	print('validate clusters against pdb and uniprot')
	hmms = glob.glob( alndir + '*.hhm')
	for hmm in hmms:
		if unifirst == True:
			#generate a more diverse MSA by using the uniprot
			functions.runHHBlits(hmm , hmm+'uniprot.hhr' , palfile= config.unipath )
			inhmm =  hmm+'uniprot.hhr.hhm'
		else:
			inhmm = hmm
		#search the PDB
		functions.runHHBlits(inhmm, hmm+'pdb70.hhr')