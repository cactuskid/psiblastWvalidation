
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
	print(fastaout[0:500])

	handle= open(config.blastdir + 'ALLresults.fasta', 'w')
	handle.write(fastaout)
	handle.close()
	print( 'sequences loaded')

elif config.blast_seqDF == True:

	#for domain specific hits... only grab relevant parts from psiblast out

	files = glob.glob(config.inputdir +'*.csv')
	dfs= []
	print(files)
	columnsstr = 'qseqid qlen slen qcovhsp sseqid staxids bitscore  evalue score pident qstart qend sstart send sseq sscinames'
	for file in files:
		dfs.append(blastcluster.blastToDF(file , columnsstr = columnsstr))

	final = pd.concat(dfs)
	print(final)
	final = final[final.evalue < .01]

	#final = final[final.length > 0]

	final['fasta'] = final.apply(lambda x : '> ' +str(x.sseqid)+'\n'+str(x.sseq).replace('-','')+'\n' , axis = 1)
	final['length'] = final.apply(lambda x : len(x.sseq.replace('-','')) , axis = 1)

	print(final)
	#final['full'] = final[['sseqid','sseq' ] ].apply(lambda x : x.sseqid+'\n'+x.sseq+'\n' , axis = 1)
	#final.sort_values(['sseqid','length'],ascending=[True,False])
	#take =longest sequence of each hit
	fastaout=''



	dfdict= {}
	with open(config.blastdir + 'ALLresults.fasta', 'w') as fastaout:
		for ID in final.sseqid.unique() :
			subdf =final[final.sseqid==ID].sort_values(['length'] ,ascending = False)
			#print(subdf.fasta.iloc[0]+'\n')
			if len(str(subdf.sseq.iloc[0].replace('-','') ))> 0 :
				fastaout.write(subdf.fasta.iloc[0])
				dfdict[ID] = {  'seq': str(subdf.sseq.iloc[0].replace('-','') ) , 'id': ID , 'fasta': subdf.fasta.iloc[0] , 'len': len(str(subdf.sseq.iloc[0].replace('-','') ))}

		fastas= glob.glob(config.inputdir + '*.fasta')
		print(fastas)
		#for fasta in fastas:
		#	with open( fasta, 'r') as fastain:
		#		fastaout.write(fastain.read())
		#	for i, seq_record in enumerate(SeqIO.parse(fasta , 'fasta')):
		#		dfdict[seq_record.id] = {  'seq': str(seq_record.seq) , 'id': seq_record.id , 'fasta': '>'+seq_record.id +'\n' + str(seq_record.seq),  'len' : len(seq_record.seq) }



	#make seqDF from all results
	seqDF = pd.DataFrame.from_dict(dfdict , orient='index')
	print(seqDF)


else:
	with open(config.datadir+ 'clusterOBjects.pkl', 'rb') as clusterobjects:
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
	alndone = glob.glob(config.alndir + '*aln.fasta')

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
		#print(len(seqDF['fasta'][seqDF.cluster == clust]))
		#print(seqDF['fasta'][seqDF.cluster == clust])
		fasta = seqDF['fasta'][seqDF.cluster == clust].tolist()
		shuffle(fasta)

		if len(fasta)> 1:
			if len(fasta) > 50:
				print('subsample to 50')
				fasta =fasta[0:50]
			else:
				oneprot+=1
			fastafile = config.alndir + str(clust) + '.fasta'
			with open( fastafile , 'w' ) as fastaout:
				fastaout.write(''.join(fasta))
			#print(fastafile)
			alnfile = config.alndir + str(clust) + 'aln.fasta'
			if config.overwrite or alnfile not in alndone:
				print(alnfile)
				p = functions.runclustalo(fastafile, alnfile, verbose= True , wait = False)
				#run all alns in parallel...
				alignjobs[clust] = p
			while len(alignjobs)> 10:
				time.sleep(1)
				keys =list( alignjobs.keys() )
				for c in keys:
					if c in alignjobs:
						poll =alignjobs[c].poll()
						if poll is not None:
							print(poll)
							print(alignjobs[c].communicate())
							del alignjobs[c]
							break
	while len(alignjobs)> 0:
		time.sleep(1)
		keys =list( alignjobs.keys() )
		for c in keys:
			if c in alignjobs:
				poll =alignjobs[c].poll()
				if poll is not None:
					print(poll)
					print(alignjobs[c].communicate())
					del alignjobs[c]

if config.hmm_compile == True:
	print( 'prepare hmms for each cluster')
	hmms = []
	models = []
	alns = glob.glob(config.alndir + '*aln.fasta')
	print(alns)
	with open(config.datadir+ 'clusterOBjects.pkl', 'rb') as clusterobjects:
		c = pickle.load(clusterobjects)
	clusters ,reverse, protlabels, index, subDM, seqDF   = c

	for aln in alns:
		#output = functions.runreformat(aln, aln.split('.fasta')[0] + '.a3m', True)
		#output = functions.runaddSS(aln.split('.fasta')[0] + '.a3m', aln.split('.fasta')[0] + '.a3m' , True)
		output = functions.runHHmake(aln,aln.split('.fasta')[0] + '.hhm' ,verbose=True)
		#hmms.append(aln.split('.fasta')[0] + '.hhm')
		#rename HMMs to cluster name
	hmms = glob.glob( config.alndir + '*.hhm')
	for hmm in hmms:
		print(hmm)
		newstr= ''
		with open(hmm, 'r') as instr:
			for line in instr:
				if 'NAME' in line:
					#q denotes query or seed sequences are present
					newstr += 'NAME ' + hmm.split('/')[-1].replace('.hhm','') +'\n'

				else:
					newstr += line

		with open(hmm, 'w') as outstr:
			outstr.write(newstr)

	#output = functions.runHHDB(config.alndir, config.alndir+'clusters', verbose = True )


if config.hmm_allvall:
	results=[]
	hmms = glob.glob(config.alndir + '*a3m')+glob.glob(config.alndir + '*aln.hhm')
	done =  glob.glob(config.alndir + '*allVall.hhr')
	for model in hmms:
		print(model)
		outfile= model+'allVall.hhr'
		if config.overwrite or outfile not in done:
			output = functions.runHHSearch(model, outfile, palfile= config.alndir+'clusters_hhm_db',  verbose= True)
		#output = functions.runHHBlits(model , outfile , config.alndir +'clusters_a3m_db' , verbose = True)
		results.append(outfile)

if config.HHDM_compile == True:
	print('compile all v all hmm results')
	results = glob.glob( config.alndir + '*allVall.hhr')
	print(results)
	names = [ n.split('/')[-1] for n in results]
	probaDM, evalDM ,pvalDM,  lenDM , scoreDM, SSDM, covDM, NX , clusternames = blastcluster.HHSearch_parseTo_DMandNX(results , None )
	with open(config.alndir + 'hmmallvall.pkl', mode='wb') as networkout:
		networkout.write(pickle.dumps([probaDM, evalDM ,pvalDM,  lenDM , scoreDM, SSDM, covDM, NX , clusternames] , 2))

	HHDM =  (1 / np.multiply(scoreDM ,lenDM	))
	HHDM = HHDM/ np.amax(HHDM)


	np.fill_diagonal(HHDM , 0 )
	#HHDM[HHDM == np.inf ] = 10**5
	HHDM *= 1000
	print(HHDM)
	matfile = config.alndir+'/HHDM.mat'
	outpath,outstr = blastcluster.distmat_to_txt(clusternames , HHDM, matfile )
	outfile = matfile+'_tree.txt'


if config.PDB70Validate == True:
	print('validate clusters against pdb and uniprot')
	hmms = glob.glob(config.alndir + '*aln.hhm')
	done = glob.glob(config.alndir + '*uniprot.hhr')
	for hmm in hmms:

		if hmm.replace('.hhm', 'uniprot.hhr' ) not in done or config.overwrite:
			print(hmm)
			functions.runHHBlits(hmm , hmm.replace('.hhm', 'uniprot.hhr' ) , palfile= config.unipath ,verbose=True )
			functions.runHHBlits( hmm.replace('.hhm', 'uniprot.hhm' ), hmm+'pdb70.hhr', palfile= config.pdb70path , iter = 1)
