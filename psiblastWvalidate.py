
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


if config.load_seqdf == True:
	files = glob.glob(config.inputdir +'*.fasta')
	print(files)

	for file in files:
		for i, seq_record in enumerate(SeqIO.parse(file, 'fasta')):
			dfdict[seq_record.id] = { 'file':file , 'full':str(seq_record) , 'seq': str(seq_record.seq) , 'id': seq_record.id , 'fasta': '>'+seq_record.id +'\n' + str(seq_record.seq) +'\n'
			, 'len' : len(seq_record.seq) , 'seq_obj': seq_record }
		else:
			print(i)
	#make seqDF from all results

	seqDF = pd.DataFrame.from_dict(dfdict , orient='index')
	print(seqDF)

	fastaout = ''.join(seqDF['fasta'].tolist())
	print(fastaout)
	handle= open(config.blastdir + 'ALLresults.fasta', 'w')
	handle.write(fastaout)
	handle.close()
	print( 'sequences loaded')
else:
	with open( 'clusterOBjects.pkl', 'rb') as clusterobjects:
	    c = pickle.load(clusterobjects)
	clusters ,reverse, protlabels,  index, subDM, seqDF   = c

"""
if config.length_filter == True:

	seqDF =seqDF[ seqDF.len > config.lower_bound]
	seqDF =seqDF[ seqDF.len < config.upper_bound]
	print(seqDF)
"""
if config.blastall == True:
	print( 'blast all v all')
	#blast all v all
	p1, allvall = functions.runBlastALL( config.fastadir+'ALLresults.fasta','blastp', 'makeblastdb', config.blastdir , mp.cpu_count() )

if config.distmat == True:

	print( 'make blast DM')
	DM, df , protlabels = blastcluster.Blast_parseTo_DM( config.fastadir + 'allVall.tab' ,  config.blastdir , E_VALUE_THRESH=10)
	print(len(protlabels))
	print(len(df.qseqid.unique()))
	#dump blast distance calc
	with open(config.datadir + 'BlastallOBjects.pkl' , 'wb')as matplickle:
		pickle.dump([DM ,df, protlabels],matplickle,-1)

if config.clusterDistmat == True:
	print( 'clustering with blast DM')
	with open(config.datadir + 'BlastallOBjects.pkl' , 'rb')as matplickle:
		DM ,df, protlabels = pickle.load(matplickle)
	clusters ,reverse, labels, index = blastcluster.cluster_distmat(DM, protlabels, clustersize = 50)
	seqDF['cluster'] = seqDF['id'].map(reverse)
	print(seqDF.cluster)
	#df['compositionScore'] = df['seq'].map(compositionScore)
	#seqDF['compositionScore'] = seqDF['id'].map(lambda x : if 'Representative' in x x.split('Representative=')[-1] else None)
	seqDF['Representative']= seqDF['id'].map(blastcluster.metaclust_filter)
	seqDF['scaffold']= seqDF['Representative'].map(blastcluster.return_scaffolds)
	seqDF['gi']= seqDF['Representative'].map(blastcluster.return_gi)
	seqDF['upi']= seqDF['Representative'].map(blastcluster.return_UPI)

	with open(config.datadir+ 'clusterOBjects.pkl' , 'wb')as matplickle:
		pickle.dump([clusters ,reverse, protlabels , index, DM, seqDF ],matplickle,-1)
	seqDF.to_csv('clusters.csv')

if config.makemodels == True:
	#make alignments and models
	alignjobs={}
	oneprot=0
	print( 'aln clusters')
	with open(config.datadir+ 'clusterOBjects.pkl', 'rb') as clusterobjects:
	    c = pickle.load(clusterobjects)
	clusters ,reverse, protlabels,  index, subDM, seqDF   = c
	clusters = seqDF['cluster'].unique()
	print(clusters)
	print(len(clusters))

	for clust in clusters:
		print(len(seqDF['fasta'][seqDF.cluster == clust]))
		print(seqDF['fasta'][seqDF.cluster == clust])

		fasta = seqDF['seq_obj'][seqDF.cluster == clust].tolist()
		shuffle(fasta)

		if len(fasta)> 1:
			if len(fasta) > 50:
				print('subsample to 50')
				fasta =fasta[0:50]
			else:
				oneprot+=1

			fastafile = config.alndir + str(clust) + '.fasta'
			with open( fastafile , 'w' ) as fastaout:
				SeqIO.write(fasta, fastaout , "fasta")
			print(fastafile)
			p = functions.runclustalo(fastafile, verbose= True , wait = False)
			#run all alns in parallel...
			alignjobs[clust] = p
		else:
			pass

	print(len(alignjobs))
	print(len(alignjobs)+oneprot)

	#join all jobs and check aln
	for clust in alignjobs:
		output = alignjobs[clust].communicate()
		aln = output[0].decode()

		print(aln)

		alnfile = config.alndir + str(clust) + 'aln.fasta'
		with open( alnfile , 'w') as alnout:
			alnout.write(functions.runclustalo(fastafile, verbose= True).replace('*', '-') )

if config.hmm_compile == True:
	print( 'prepare hmms for each cluster')

	hmms = []
	models = []
	alns = glob.glob(config.alndir + '*aln.fasta')
	print(alns)
	with open(config.datadir+ 'clusterOBjects.pkl', 'rb') as clusterobjects:
		c = pickle.load(clusterobjects)
	clusters ,reverse, protlabels, index, subDM, seqDF   = c
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
	hmms = glob.glob( config.alndir + '*.hhm')
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
	output = functions.runHHDB(config.alndir, config.alndir+'clusters', verbose = True )


if config.hmm_allvall:
	results=[]
	hmms = glob.glob(config.alndir + '*a3m')+glob.glob(config.alndir + '*hhm')
	for model in hmms:
		print(model)
		outfile= model+'allVall.hhr'
		output = functions.runHHSearch(model, outfile, palfile= config.alndir+'clustersdb_hhm_db',  verbose= True)
		#output = functions.runHHBlits(model , outfile , config.alndir +'clusters_a3m_db' , verbose = True)
		results.append(outfile)

if config.HHDM_compile == True:
	print('compile all v all hmm results')
	results = glob.glob( config.alndir + '*.allvall.hhr')
	probaDM, evalDM ,pvalDM,  lenDM , scoreDM, SSDM, NX , clusternames = blastcluster.HHSearch_parseTo_DMandNX(results , None )
	HHDM =  1 / ( (scoreDM +  3*SSDM) )
	np.fill_diagonal(HHDM , 0 )
	HHDM[HHDM == np.inf ] = 10**5
	HHDM /= 1000
	matfile = config.alndir+'/HHDM.mat'
	outpath,outstr = blastcluster.distmat_to_txt(clusternames , HHDM, matfile )
	outfile = matfile+'_tree.txt'
	#output = functions.runFastme(matfile, outfile)
	#print(output)

if config.PDB70Validate == True:
	print('validate clusters against pdb and uniprot')
	hmms = glob.glob( config.alndir + '*.a3m')+glob.glob(config.alndir + '*hhm')
	for hmm in hmms:
		if config.unifirst == True:
			print(hmm)
			#generate a more diverse MSA by using the uniprot
			functions.runHHBlits(hmm , hmm.replace('.a3m', 'uniprot.hhr' ) , palfile= config.unipath ,verbose=True )
			inhmm =  hmm+'uniprot.hhr.hhm'
		else:
			inhmm = hmm
		#search the PDB
		functions.runHHBlits(inhmm, hmm+'pdb70.hhr', palfile= config.pdb70path , iter = 1)
