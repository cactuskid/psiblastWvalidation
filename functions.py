

import subprocess
import shlex
from tempfile import NamedTemporaryFile as namedTemp
import tempfile as tmp
import os
from Bio import Entrez
import multiprocessing as mp
import config


def openprocess(args , inputstr =None , verbose = False , wait = True):
	args = shlex.split(args)
	p = subprocess.Popen(args,  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr= subprocess.PIPE )
	if verbose == True and inputstr is not None:
		print(inputstr.decode())

	if inputstr != None:
		p.stdin.write(inputstr.encode())
	if wait == True:
		output = p.communicate()
		if verbose == True:
			print(output)
		p.wait()
		return output[0].decode()
	else:
		return p
def runBlast(query,db, verbose):
    cmdstr = 'blast -query '+query+ ' -db '+ db
    return (openprocess(cmdstr, verbose=verbose))

def runHHDB(hmmfold,name, verbose):
    cmdstr =  '/scratch/cluster/monthly/dmoi/hhsuite-2.0.16/scripts/hhblitsdb.pl -ihhm '+hmmfold  + ' -o '+ hmmfold+ name
    return (openprocess(cmdstr, verbose=verbose))

def runBlastDB(inputfile,verbose):
    cmdstr = 'makeblastdb -in '+ inputfile + '-dbtype protein '
    return (openprocess(cmdstr , verbose=verbose))

def dlGI(Ginum):
	handle = Entrez.efetch(db="protein", id=Ginum, rettype="gb", retmode="gb")
	return handle.read()

def runFastme(clusterfile ,outfile ):
	cmdstr =  'fastme -i ' + clusterfile + ' -o ' + clusterfile+'_tree.txt'
	return (openprocess(cmdstr) )

def runmafft(fasta, verbose=False):
	cmdstr = 'mafft --localpair --maxiterate 1000  --thread -1 '+fasta
	return (openprocess(cmdstr, verbose = verbose) )

def runclustalo(fasta, inputstr= None,verbose=False , wait = True):
	cmdstr = './clustalo --auto  -i '+fasta + ' --threads='+ str(int(mp.cpu_count()/2))
	return (openprocess(cmdstr, inputstr=inputstr, verbose = verbose , wait=wait) )

def runHHmake(aln,outfile ,verbose=False):
	cmdstr = 'hhmake -add_cons -M 50 -i '+aln + ' -o ' +outfile
	return (openprocess(cmdstr,verbose=verbose))

def runaddSS(aln,outfile ,verbose=False):
	# addss.pl <ali_file> [<outfile>] [-fas|-a3m|-clu|-sto]
	cmdstr = config.scriptdir+'addss.pl '+aln +' '+ outfile +  ' -fas '
	return (openprocess(cmdstr,verbose=verbose))
def runreformat(aln,outfile ,verbose=False):
	#reformat.pl [informat] [outformat] infile outfile [options]
	cmdstr = config.scriptdir+'reformat.pl fas a3m '+aln +' '+ outfile
	return (openprocess(cmdstr,verbose=verbose))

def runHHSearch(aln,outfile, palfile  ,verbose= False):
	cmdstr = '/scratch/cluster/monthly/dmoi/hhsuite-2.0.16/bin/hhsearch -cpu '+ str(mp.cpu_count()) + ' -i '+aln +' -d '+ palfile + ' -o ' +outfile +' -e 1 -p 0 -realign'
	return (openprocess(cmdstr, verbose = verbose))

def runHHBlits(aln,outfile, palfile , iter =3, verbose= False):
	#output the model by default
	#TODO add options for model outputfiles
	cmdstr = 'hhblits -n '+ str(iter) + ' -cpu '+ str(mp.cpu_count()) +' -i '+aln +' -d '+ palfile + ' -o ' +outfile + ' -ohhm ' + outfile.replace( '.hhr', '.hhm' )
	return (openprocess(cmdstr, verbose = verbose))

def runBlastALL(fastafile,blastpath, formatpath, outdir , nCPU ):
	#run blast all against all in a sequence file
	allVall= outdir + fastafile.split('/')[-1].replace('.fasta', '_allVall.tab')
	seqdb = outdir + fastafile.split('/')[-1].replace('.fasta', '.db')
	args1 =  formatpath +' -in '+ fastafile + ' -out ' + seqdb + ' -dbtype prot'
	print( args1)
	args1 = shlex.split(args1)
	p1 = subprocess.call(args1)

	#split queries into tempfile for multiprocessing:
	# count entries
	count = 0
	with open(fastafile, 'r') as infile:
		for line in infile:
			if '>' in line:
				count+=1

	print( str(count) + ' sequences in dataset ' )
	seqsperfile = int(count / nCPU)
	tempstr = ''
	handles = []
	processes = []
	outputfiles = []
	filecount = 0
	seqcount = 0

	#split up blast jobs
	with open(fastafile, 'r')as infile:
		for i,line in enumerate(infile):
			if '>' in line:
				seqcount+=1
			if seqcount== seqsperfile:
				#create tempfile
				seqcount =0
				filehandle,filename = tmp.mkstemp( 'w' , dir = outdir)
				with open(filename , 'w') as tempfile:
					tempfile.write(tempstr)
				query = filename
				output = allVall+str(filecount)
				filecount +=1
				args2 = blastpath+' -db '+seqdb+' -query '+query+' -outfmt="10 qseqid sseqid qstart qend sstart send length evalue sseq"   -out ' + output
				print( args2)
				p2 = subprocess.Popen(shlex.split(args2)  ,  stdout=subprocess.PIPE )
				processes.append(p2)
				outputfiles.append(output)
				tempstr = ''
			tempstr += line
		#leftovers
	filehandle,filename = tmp.mkstemp( 'w' , dir = outdir)
	with open(filename , 'w') as tempfile:
		tempfile.write(tempstr)
	query = filename
	output = allVall+str(filecount)
	args2 = blastpath+' -db '+seqdb+' -query '+query+' -outfmt="10 qseqid sseqid qstart qend sstart send length evalue sseq"  -out ' + output
	print( args2)
	p2 = subprocess.Popen(shlex.split(args2) ,  stdout=subprocess.PIPE  )
	processes.append(p2)
	outputfiles.append(output)

	#wait until everything is finished
	for process in processes:
		process.communicate()

	print( 'done blasting')

	#concatenate output and clean up the files
	allVallfinal = outdir+'allVall.tab'
	with open( allVallfinal , 'w') as finalout :
		for output in outputfiles:
			results = open(output, 'r')
			finalout.write(results.read())
			results.close()
			os.remove(output)
	return p1 , [allVallfinal]
