import pandas as pd
import glob
import blastcluster

from Bio import SeqIO


datadir = '/scratch/cluster/monthly/dmoi/fusexins_all/'
csvdir = '/scratch/cluster/monthly/dmoi/fusexins_all/input2/'
metaclustdir = '/scratch/cluster/monthly/dmoi/metaclust_new/'
scaffolddir = '/scratch/cluster/monthly/dmoi/archaea/scaffolds/'

csvs = glob.glob( config.csvdir + '*.csv')

dfs= []
print(csvs)

for file in csvs:
    try:
      df = blastcluster.blastToDF(file, metaclust = True)
      dfs.append(df)
    except:
      print(file)

blastdf = pd.concat(dfs, axis = 0)
blastdf['evalue']=blastdf['evalue'].convert_objects(convert_numeric=True)
blastdf = blastdf[ blastdf['evalue'] < .001]
blastdf = blastdf[ blastdf['qlen'] > 300 ]

blastdf['grab']= blastdf['sseqid'].map( lambda x: x.split('|')[0])

print(blastdf)

grablist = blastdf.sseqid.unique().tolist()
print(len(grablist))
records = SeqIO.parse(config.metaclustdir + "metaclust50_2017_01.fasta", "fasta")
print(next(records).id)
count =0
with open(config.datadir + 'filterpsiblast.fasta', 'w') as metaout:
    for record in records:
        if record.id in grablist:
            metaout.write('>'+str(record.id) +'\n'+ str(record.seq)+'\n' )
            count+=1
            print(count)
            if count == len(grablist):
                break
