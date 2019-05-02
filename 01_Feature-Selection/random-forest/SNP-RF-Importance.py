###====== LIBRARIES =========================================================###
import numpy as np
import argparse, re
from sklearn.ensemble import RandomForestClassifier
try:
    from snp import alphabet,animinfo,snpinfo,snp
except ImportError:
    print("Please put the 'snp.py'-Library in the same folder!")
    exit()






###====== FUNCTIONS =========================================================###
def read_mapped(InMap, InPed):
    MapList = []
    PedList = []
    MapLineCount = 0
    PedLineCount = 0

    InMapFile = open(InMap,"r")
    print("\r[  ] Processing MAP-File",end="")
    for idx,Line in enumerate(InMapFile):
        print("\r[  ] Processing MAP-File: line {}".format(idx),end="")
        Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
        if len(Entries) != 4:
            print("!!WARNING!! Line {} in MAP-File has to many entries!\n".format(idx)+Line)
        MapList.append(snpinfo(Entries[0],Entries[1],Entries[2],Entries[3]))
        MapLineCount += 1
    print("\r{}\r[OK] Processing MAP-File".format(" "*80))
    InMapFile.close()

    SnpList = [snp(SnpInfo) for SnpInfo in MapList]
    del MapList

    InPedFile = open(InPed,"r")
    print("\r[  ] Processing PED-File",end="")
    for adx,Line in enumerate(InPedFile):
        print("\r[  ] Processing PED-File: line {}".format(adx),end="")
        Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
        if (len(Entries)-6)/2 != MapLineCount:
            print("!!WARNING!! Line {} has a different amount({}) of snps than required({})!".format(adx,(len(Entries)-6)/2,MapLineCount))
            continue
        PedList.append(animinfo(Entries[0],Entries[1],Entries[2],Entries[3],Entries[4],Entries[5]))
        for idx,jdx in enumerate(range(6,len(Entries),2)):
            SnpList[idx].append(alphabet.encode("{}{}".format(Entries[jdx],Entries[jdx+1])))
        PedLineCount += 1 
    print("\r{}\r[OK] Processing PED-File".format(" "*80))
    InPedFile.close()
    return SnpList
def BuildSnpUnsupRfData(SnpData):
    DataLength = len(SnpData[0])
    SynthData = np.append( np.array(SnpData).T, np.array([ np.random.choice(Row, size=len(Row),replace=True) for Row in np.array(SnpData) ]).T,axis=0)
    SynthLabel = np.repeat([0,1], [DataLength,DataLength],axis=0)
    return SynthData, SynthLabel




###====== VARIABLES =========================================================###
OutMapFileNameCopy = "data.map"
OutPedFileNameCopy = "data.ped"
OutAttrFileNameCopy = "data.attr"

NTree = 500
NThread = -1 #all!




###====== MAIN ==============================================================###
Parser = argparse.ArgumentParser(description='Finds for each SNP the importance via variable importance of unsupervised random forest; Calculates per chromosome.')
Parser.add_argument('--map', metavar='{MAP-File-Name}', type=str, nargs=1, help='Input MAP file, where a subset is taken from', required = True)
Parser.add_argument('--ped', metavar='{PED-File-Name}', type=str, nargs=1, help='Input PED file, where a subset is taken from', required = True)
Parser.add_argument('-n', metavar='<Sample-Set-Sizes>', type=float, nargs=1, help='Return the best n*100 \%  SNPs', required = True)
Parser.add_argument('--tree', metavar='<Tree-Number>', type=int, nargs=1, help='Number of Trees used for random forest, default tree = {}'.format(NTree), required = False)
Parser.add_argument('-t', metavar='<Thread-Number>', type=int, nargs=1, help='Number threads to be used for the random forest, default tree = {} (all)'.format(NTree), required = False)
Parser.add_argument('--omp', metavar='{MAP-Out-Name}', type=str, nargs=1, help='Output-MAP file name', required = False)
Parser.add_argument('--opd', metavar='{PED-Out-Name}', type=str, nargs=1, help='Output-PED file name', required = False)
Parser.add_argument('--out', metavar='{Attribut-List}', type=str, nargs=1, help='Attribute importance list file name', required = False)
Parser.add_argument('--genome', action="store_true", help='Calculate attribute importance on entire genome, not by chromosome', required = False)
Args = Parser.parse_args()
if Args.tree:
    NTree = Args.tree[0]
if Args.t:
    NThread = Args.t[0]
if Args.omp:
    OutMapFileNameCopy = Args.omp[0]
if Args.opd:
    OutPedFileNameCopy = Args.opd[0]
if Args.out:
    OutAttrFileNameCopy = Args.out[0]
OnGenome = Args.genome
SnpPercent = Args.n




#--- Read SNP -----------------------------------------------------------------#
SnpList = read_mapped(Args.map[0],Args.ped[0])

if OnGenome:
    SnpIDs = [repr(Snp) for Snp in SnpList]

    #--- Make Synthetic Data as mentioned by Leo Breiman for unsup. RF ------------#
    print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
    SnpSynthRfData , SnpSynthRfLabel = BuildSnpUnsupRfData(SnpList)
    print("\r[OK] Creating synthetic Data using SNPs for random forest")


    #--- Perform unsupervised random forest ---------------------------------------#
    print("\r[  ] Performing unsupervised random forest",end="")
    SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread)
    SnpRfSynthClf.fit(SnpSynthRfData, SnpSynthRfLabel)
    print("\r[OK] Performing unsupervised random forest")



    #--- Calculate Variable importance and store ----------------------------------#
    print("\r[  ] Writing important attributes to file, sorted",end="")
    OutAttrFile = open(OutAttrFileNameCopy, "w")
    Importance = sorted( [ (ID,importance) for ID, importance in zip(SnpIDs, SnpRfSynthClf.feature_importances_)], key=lambda x:x[1], reverse = True)
    for Tup in Importance:
        OutAttrFile.write("{}\t{}\n".format(Tup[0],Tup[1]))
    OutAttrFile.close()
    print("\r[OK] Writing important attributes to file, sorted")


else:
    SnpChrom = [ Snp.info().chrom_id() for Snp in SnpList ]
    ChromList = sorted(list(set(SnpChrom)))
    ChromListLen = len(ChromList)
    Importance = []
    for Chrom in ChromList:
        SubSnpList = [SnpList[Idx] for Idx,Chr in enumerate(SnpChrom) if Chr == Chrom]
        SubSnpIDs = [repr(Snp) for Snp in SubSnpList]
        
        #--- Make Synthetic Data as mentioned by Leo Breiman for unsup. RF ------------#
        print("\r[  ] Chromosome {}:{} [1/3]: Creating synthetic Data using SNPs for random forest".format(Chrom,ChromListLen),end="")
        SnpSynthRfData , SnpSynthRfLabel = BuildSnpUnsupRfData(SubSnpList)

        #--- Perform unsupervised random forest ---------------------------------------#
        print("\r{}\r".format(" "*81),end="")
        print("\r[  ] Chromosome {}:{} [2/3]: Performing unsupervised random forest".format(Chrom,ChromListLen),end="")
        SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread)
        SnpRfSynthClf.fit(SnpSynthRfData, SnpSynthRfLabel)

        print("\r{}\r".format(" "*81),end="")
        print("\r[  ] Chromosome {}:{} [3/3] Calculating and storing attribute importance".format(Chrom,ChromListLen),end="")
        SubImportance = sorted( [ (ID,Chrom,Importance) for ID, Importance in zip(SubSnpIDs, SnpRfSynthClf.feature_importances_)], key=lambda x:x[1], reverse = True)
        Importance.extend(sorted(SubImportance))
        print("\r{}\r".format(" "*81),end="")
        print("\r[OK] Chromosome {}:{} unsupervised RF and attribute importance".format(Chrom,ChromListLen))
    print("\r[  ] Writing important attributes to file, sorted",end="")
    OutAttrFile = open(OutAttrFileNameCopy, "w")
    for Tup in Importance:
        OutAttrFile.write("{}\t{}\t{}\n".format(Tup[0],Tup[1],Tup[2]))
    OutAttrFile.close()
    print("\r[OK] Writing important attributes to file, sorted")