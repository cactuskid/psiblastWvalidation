import pandas as pd
import glob
import blastcluster

from Bio import SeqIO
import config
import functools
import multiprocessing as mp
import dill as pickle
#### filter huge fastas with mp






if __name__ == "__main__":
    csvdir = '/scratch/cluster/monthly/dmoi/archaea/psiblast_wProfileValidation/input/'
    csvs = glob.glob( csvdir + '*.csv')
    dfs= []
    print(csvs)
    for file in csvs:
      df = blastcluster.blastToDF(file, columnsstr =' qseqid qlen slen qcovhsp sseqid staxids bitscore score evalue pident qstart qend sstart send staxids sscinames' )
      dfs.append(df)

    blastdf = pd.concat(dfs, axis = 0)
    print(blastdf.qseqid.unique().tolist())
    blastdf['evalue']=blastdf['evalue'].infer_objects()
    subdf = blastdf[ blastdf.evalue < .00001]

    #print(subdf)

    grablist = subdf.sseqid.unique().tolist()
    print(len(grablist))
    print(grablist[0:100])

    count =0
    record = False
    linecount = 0
    filename = config.metaclustdir + "metaclust_50.fasta"
    output = './filterpsiblast_fusexinsall.fasta'


    def seqgenerator( filename , chunksize=1000):
        with open(filename,  'r')as metain:
            linecount=0
            retstr=''
            retcount = 0
            for i,line in enumerate(metain):
                if '>' in line:
                    linecount+=1
                if linecount >chunksize:
                    linecount=0
                    retcount+=1
                    yield retstr
                    retstr=''
                retstr+=line+'\n'
    def grab(lines, grablist):
            record = False
            thisstr =''
            retarray = []

            for line in lines.split('\n'):
                if '>' in line:
                    if record == True:
                        record = False
                        retarray.append(thisstr)
                        thisstr =''
                    if '#' in line:
                        if line.replace('>', '').split('#')[0].strip() in grablist:
                            record = True


                    elif '|' in line:
                        if line.replace('>', '').split('|') in grablist:
                            record = True

                if record == True:
                    thisstr += line
            return retarray

    seqgen = seqgenerator(filename, chunksize=10000)
    grabentries = functools.partial(grab, grablist= grablist)

    N = 1000
    result= []
    i=0
    saveinterval = 1

    import itertools
    import time
    outstr = ''
    pool = mp.Pool( int(mp.cpu_count()))
    while True:
        g2 = pool.map(grabentries, itertools.islice(seqgen, N))

        if g2:
            result+= [ str for str in [res for res in g2 if len(res) > 0] ]
            print(len(result))
            i+=1
            if i %saveinterval == 0:
                compile =  ''.join( [ res[0]+'\n' for res in result ])

                with open(output, 'w') as metaout:
                    metaout.write(compile)
        else:
            break
