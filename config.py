#datadir = '/scratch/cluster/monthly/dmoi/sas6/'
csvdir = '/scratch/cluster/monthly/dmoi/archaea/identify/'
datadir = '/scratch/cluster/monthly/dmoi/archaea/identify/'

metaclustdir = '/scratch/cluster/monthly/dmoi/metaclust50/'
scaffolddir = '/scratch/cluster/monthly/dmoi/archaea/scaffolds/'

inputdir = '/scratch/cluster/monthly/dmoi/sas6/input/'
alndir = '/scratch/cluster/monthly/dmoi/sas6/aln/'
blastdir = '/scratch/cluster/monthly/dmoi/sas6/blastdir/'

#inputdir = './input/'
#alndir = './aln/'
#blastdir = './blastdir/'


fastadir = blastdir
subsample_aln = False

unipath = '/db/SOFTWARE/hhsuite/uniclust_2017_10/uniclust30_2017_10_hhsuite/uniclust30_2017_10'
pdb70path = '/db/SOFTWARE/hhsuite/pdfb70/pdb70'


scriptdir = '/scratch/cluster/monthly/dmoi/hhsuitscritps/'
#scriptdir = '/usr/lib/hhsuite/scripts/'
seeds = ['./input/archaeaPositives.fasta' , './input/interprohap2.fasta' , './input/clsuter9_13.fasta' ]
subsample_aln = False

anchordir = './anchormodels/'

#turn on or off parts of the pipeline
load_seqdf = False
blast_seqDF = True

#filter input
length_filter = True
upper_bound = 2000
lower_bound = 300

#blast all v all
blastall = True

#create distance matrix from blast scores
distmat = True

#use kmeans to create clusters
clusterDistmat = True

#align clusters
makemodels = True
#turn alignments into hhm format
hmm_compile = True
hmm_allvall = True
#hmm all v all search
HHDM_compile = True

#validate clusters with external dbs

unifirst = True
PDB70Validate = False
