

csvdir = '/scratch/cluster/monthly/dmoi/archaea/metaclustcsv/'
metaclustdir = '/scratch/cluster/monthly/dmoi/metaclust/'
scaffolddir = '/scratch/cluster/monthly/dmoi/archaea/scaffolds/'

inputdir = './input/'
alndir = './aln/'
blastdir = './blastdir/'

subsample_aln = False

unipath = '/db/SOFTWARE/hhsuite/uniclust_2017_10/'
pdb70path = '/db/SOFTWARE/hhsuite/pdfb70/'

scriptdir = './hhsuitscritps/'
#scriptdir = '/usr/lib/hhsuite/scripts/'
seeds = ['./input/archaeaPositives.fasta' , './input/interprohap2.fasta' , './input/clsuter9_13.fasta']
subsample_aln = False

anchordir = './anchormodels/'

#turn on or off parts of the pipeline
blastall = False
distmat = False
clusterDistmat = True
makemodels = True
HMMall = True
compileDB = True
HHDM_compile = True
PDB70Validate = True
