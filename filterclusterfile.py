import pandas as pd
import pickle

clusterout = {}

hits = '/scratch/cluster/monthly/dmoi/fusexins_all/pablosclusters.csv'
clusterfile = '/scratch/cluster/monthly/dmoi/metaclust50/metaclust_50_cluster.fasta'


csvout = './fusexinsout_out.csv'
csvscaffolds = './fusexin_scaffolds.csv'

#output scaffold files here
scaffolds_path = './'

grabhits = False
grabscaffolds = False

formatDf = True
outputscaffolds = True



hitsdf = pd.DataFrame.from_csv(hits)
print(hitsdf)

hitlist = hitsdf.id.unique().tolist()
hitlist = [ entry.replace('>','').strip() for entry in hitlist ]
hitlist = set(hitlist)
print(hitlist)


foundcout = 0
clusters = 0
filenum = 0



record = False
startnum = 0


print(hitlist)
print('pulling hits')
#first pass. grab all fusexin clusters
if grabhits:
    with open('globalout_rev2.fasta' , 'w') as globalout:
        with open( clusterfile, 'r') as clusterin:
            for i,line in enumerate(clusterin):
                if i > startnum:
                    if i % 10000000 == 0:
                        print(i)
                        print('sequences found:'+str(len(clusterout)))
                        print('clusters scanned:'+str(clusters))

                    if '>' in line and '#' not in line:
                        if i < 100:
                            print(line.replace('>','') )
                        clusters+=1
                        if  line.replace('>','').strip()  in hitlist:
                            record = True
                            representative = line.replace('>','')
                            foundcout+=1
                            entry_name = representative
                            clusterout[entry_name] = {}
                            clusterout[entry_name]['representative'] = representative
                            clusterout[entry_name]['sequence'] = ''

                            print(line)
                            globalout.write(line)

                            if len(clusterout)>0 and len(clusterout)% 100 == 0 :
                                print('saving sequences')
                                seqDf = pd.DataFrame.from_dict(clusterout, orient='index')
                                seqDf.to_csv(csvout)
                                print(seqDf)
                        else:
                            record = False
                    elif record==True and '>' in line and '#' in line:
                        #cluster content
                        print(line)
                        globalout.write(line)

                        entry_name = line.replace('>','')
                        clusterout[entry_name] = {}
                        clusterout[entry_name]['representative'] = representative
                        clusterout[entry_name]['sequence'] = ''

                    elif record==True and '>' not in line:
                        #sequence data
                        clusterout[entry_name]['sequence'] += line
                        globalout.write(line)
            else:
                seqDf = pd.DataFrame.from_dict(clusterout, orient='index')
                print(seqDf)
                seqDf.to_csv(csvout)

#second pass grab all prots from all scaffolds found in first pass
print('done pulling hits')
print('pulling scaffolds')
if grabscaffolds:
    seqDf = pd.DataFrame.from_csv(csvout)
    hitlist = seqDf.index.unique().tolist()
    hitlist = [ entry.split('#')[0] if '#' in entry else entry  for entry in hitlist]
    hitlist = set([ ''.join(entry.replace('>','').strip().split('_')[:-1]) for entry in hitlist ])
    clusterout ={}
    print(hitlist)
    with open('scaffoldout.fasta' , 'w') as globalout:
        with open( clusterfile, 'r') as clusterin:
            for i,line in enumerate(clusterin):
                if i > startnum:
                    if i % 10000000 == 0:
                        print(i)
                        print('sequences found:'+str(len(clusterout)))
                        print('clusters scanned:'+str(clusters))

                    if '>' in line and '#' not in line:
                        clusters+=1
                        representative = line

                    if '>' in line and '#' in line:
                        record = False

                        instr = line.split('#')[0]
                        if i < 100:
                            print(''.join(instr.replace('>','').strip().split('_')[:-1] ) )

                        if ''.join(instr.replace('>','').strip().split('_')[:-1]) in hitlist:
                            record = True

                            foundcout+=1
                            entry_name = line.replace('>','')
                            clusterout[entry_name] = {}
                            clusterout[entry_name]['representative'] = representative
                            clusterout[entry_name]['sequence'] = ''
                            print(line)
                            globalout.write(line)

                            if len(clusterout)>0 and len(clusterout)% 100 == 0 :
                                print('saving sequences')
                                seqDf = pd.DataFrame.from_dict(clusterout, orient='index')
                                seqDf.to_csv(csvscaffolds)
                                print(seqDf)
                    
                    elif record==True and '>' not in line:
                        #sequence data
                        clusterout[entry_name]['sequence'] += line
                        globalout.write(line)

            else:
                seqDf = pd.DataFrame.from_dict(clusterout, orient='index')
                seqDf.to_csv(csvscaffolds)
                print(seqDf)


                print('done pulling scaffolds')

if formatDf:
    seqDF = pd.DataFrame.from_csv(csvscaffolds)
    seqDF['JGIproject'] = seqDF.index.map( lambda x: x.split('_')[0])
    seqDF['JGIprojectNumber'] = seqDF.index.map( lambda x: x.split('_')[1].split('.')[0])
    seqDF['scaffoldID'] = seqDF.index.map( lambda x: x.split('ID=')[1].split(';')[0] if 'ID' in x else None )
    seqDF['ORFNR'] = seqDF.scaffoldID.map( lambda x: x.split('_')[1] if x  else None )
    seqDF['scaffoldNR'] = seqDF.scaffoldID.map( lambda x: x.split('_')[0] if x  else None )
    seqDF['fullID'] = seqDF.index.map( lambda x: '>'+x  )
    seqDF['fasta'] = seqDF[['fullID','sequence']].apply( lambda x: str(x[0])+'\n'+str(x[1])+'\n' , axis = 1)
    seqDF.to_csv(csvscaffolds +'reformat.csv')


if outputscaffolds:
    for JGIproj in seqDF.JGIproject.unique():
        subdf = seqDF[ seqDF.JGIproject == JGIproj]
        for scaffold in subdf.scaffoldNR.unique():
            scaffolddf = subdf[ subdf.scaffoldNR == scaffold]
            if len(scaffolddf.index)>1:
                name = JGIproj+'scaffold_'+scaffold
                scaffolddf.to_csv(scaffolds_path+name + '.csv')
                fastastr = ''.join(seqDF.fasta.tolist())
                print(scaffolddf)
                with open(scaffolds_path+name+'.fasta', 'w')as fastout:
                    fastout.write(fastastr)

print('DONE!')
