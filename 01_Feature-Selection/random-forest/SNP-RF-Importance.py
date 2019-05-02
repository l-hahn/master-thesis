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
Parser = argparse.ArgumentParser(description='Finds for each SNP the importance via variable importance of unsupervised random forest.')
Parser.add_argument('--map', metavar='{MAP-File-Name}', type=str, nargs=1, help='Input MAP file, where a subset is taken from', required = True)
Parser.add_argument('--ped', metavar='{PED-File-Name}', type=str, nargs=1, help='Input PED file, where a subset is taken from', required = True)
Parser.add_argument('-n', metavar='<Sample-Set-Sizes>', type=float, nargs=1, help='Return the best n*100 \%  SNPs', required = True)
Parser.add_argument('--tree', metavar='<Tree-Number>', type=int, nargs=1, help='Number of Trees used for random forest, default tree = {}'.format(NTree), required = False)
Parser.add_argument('-t', metavar='<Thread-Number>', type=int, nargs=1, help='Number threads to be used for the random forest, default tree = {} (all)'.format(NTree), required = False)
Parser.add_argument('--omp', metavar='{MAP-Out-Name}', type=str, nargs=1, help='Output-MAP file name', required = False)
Parser.add_argument('--opd', metavar='{PED-Out-Name}', type=str, nargs=1, help='Output-PED file name', required = False)
Parser.add_argument('--out', metavar='{Attribut-List}', type=str, nargs=1, help='Attribute importance list file name', required = False)
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

SnpPercent = Args.n






SnpList = read_mapped(Args.map[0],Args.ped[0])
SnpNames = [repr(Snp) for Snp in SnpList]

print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
SnpSynthRfData , SnpSynthRfLabel = BuildSnpUnsupRfData(SnpList)
print("\r[OK] Creating synthetic Data using SNPs for random forest")

print("\r[  ] Performing unsupervised random forest",end="")
SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread)
SnpRfSynthClf.fit(SnpSynthRfData, SnpSynthRfLabel)
print("\r[OK] Performing unsupervised random forest")

print("\r[  ] Writing important attributes to file, sorted",end="")
OutAttrFile = open(OutAttrFileNameCopy, "w")
Importance = sorted( [ (name,importance) for name, importance in zip(SnpNames, SnpRfSynthClf.feature_importances_)], key=lambda x:x[1], reverse = True)
for Tup in Importance:
    OutAttrFile.write("{}\t{}\n".format(Tup[0],Tup[1]))
OutAttrFile.close()
print("\r[OK] Writing important attributes to file, sorted")