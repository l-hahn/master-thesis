###====== LIBRARIES =========================================================###
import numpy as np
import argparse
from sklearn.ensemble import RandomForestClassifier
try:
    from snp import alphabet,animinfo,snpinfo,snp
except ImportError:
    print("Please put the 'snp.py'-Library in the same folder!")
    exit()






###====== FUNCTIONS =========================================================###
def BuildSnpUnsupRfData(SnpData):
    DataLength = len(SnpData[0])
    SynthData = np.append( np.array(SnpData).T, np.array([ np.random.choice(Row, size=len(Row),replace=True) for Row in np.array(SnpData) ]).T,axis=0)
    SynthLabel = np.repeat([0,1], [DataLength,DataLength],axis=0)
    return SynthData, SynthLabel





###====== VARIABLES =========================================================###
OutFilePrefix = "data"

NTree = 500
NThread = -1 #all!
NFeature = "auto"
Tped = False




###====== MAIN ==============================================================###
Parser = argparse.ArgumentParser(description='Finds for each SNP the importance via variable importance of unsupervised random forest; Calculates per chromosome.')
Parser.add_argument('--file', metavar='{File-Name-Prefix}', type=str, nargs=1, help='Input file prefix for <file>.map and <file>.ped, where a subset is taken from', required = True)
Parser.add_argument('-n', metavar='<Sample-Set-Sizes>', type=int, nargs=1, help='Return the best n SNPs (per chromosome/genome)', required = True)
Parser.add_argument('--tree', metavar='<Tree-Number>', type=int, nargs=1, help='Number of Trees used for random forest, default tree = {}'.format(NTree), required = False)
Parser.add_argument('--features', metavar='<Feature-Number>', type=int, nargs=1, help='Number of feautures used for random forest, default features = {}'.format(NFeature), required = False)
Parser.add_argument('-t', metavar='<Thread-Number>', type=int, nargs=1, help='Number threads to be used for the random forest, default tree = {} (all)'.format(NThread), required = False)
Parser.add_argument('--omp', metavar='{MAP-Out-Name}', type=str, nargs=1, help='Output file name for <omp>.map and <omp>.ped for the n most important SNPs, default out = data.*', required = False)
Parser.add_argument('--out', metavar='{Attribut-List}', type=str, nargs=1, help='Attribute importance list file name, default out = <file>.attr', required = False)
Parser.add_argument('--genome', action="store_true", help='Calculate attribute importance on entire genome, not by chromosome', required = False)
Parser.add_argument('--supervised', action="store_true", help='Consider in addition phenotyp for variable importance during random forest', required = False)
Parser.add_argument('--enc', action="store_true", help='Each SNP will be additive encoded (0/1/2).', required = False)
Parser.add_argument('--tped', action="store_true", help='Use existing <file>.tped and <file>.tfam file.', required = False)
Args = Parser.parse_args()
OutAttrFile = Args.file[0]+".attr"
if Args.tree:
    NTree = Args.tree[0]
if Args.t:
    NThread = Args.t[0]
if Args.omp:
    OutFilePrefix = Args.omp[0]
if Args.out:
    OutAttrFile = Args.out[0]
if Args.features:
    NFeature = Args.features[0]
Tped = Args.tped
OnGenome = Args.genome
EncodeAdditive = Args.enc
Supervised = Args.supervised




#--- Read SNP -----------------------------------------------------------------#
SnpList, PedList = snp.read_mapped(Args.file[0], EncodeAdd = EncodeAdditive, Verbose=True, TPed = Tped)
TopSnps = min(max(Args.n[0],1),len(SnpList))
TopImpSnps = []


if OnGenome:
    SnpIDs = [repr(Snp) for Snp in SnpList]
    if NFeature != "auto":
        NFeature = min(NFeature, len(SnpIDs))

    #--- Make Synthetic Data as mentioned by Leo Breiman for unsup. RF ------------#
    if not Supervised:
        print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
        SnpSynthRfData , SnpSynthRfLabel = BuildSnpUnsupRfData(SnpList)
        print("\r[OK] Creating synthetic Data using SNPs for random forest")
    else:
        SnpSynthRfData , SnpSynthRfLabel = np.array(SnpList).T, np.array([ Animal.phenotype() for Animal in PedList])

    #--- Perform unsupervised random forest ---------------------------------------#
    print("\r[  ] Performing unsupervised random forest",end="")
    SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread, max_features = NFeature)
    SnpRfSynthClf.fit(SnpSynthRfData, SnpSynthRfLabel)
    print("\r[OK] Performing unsupervised random forest")



    #--- Calculate Variable importance and store ----------------------------------#
    print("\r[  ] Writing important attributes to file, sorted",end="")
    OutAttrFile = open(OutAttrFile, "w")
    OutAttrFile.write("ID\tVarImp\tConserv\tMAF\tEntropy\n")
    Importance = sorted( [ (ID,importance,SnpList[idx].conservation(),SnpList[idx].maf(),SnpList[idx].entropy()) for idx,(ID,importance) in enumerate(zip(SnpIDs, SnpRfSynthClf.feature_importances_))], key=lambda x:x[1], reverse = True)
    TopImpSnps = [ SnpList[ SnpIDs.index(SnpID[0]) ] for SnpID in Importance[0:TopSnps]]
    for Tup in Importance:
        OutAttrFile.write("{}\t{}\t{}\t{}\t{}\n".format(Tup[0],Tup[1],Tup[2],Tup[3],Tup[4]))
    OutAttrFile.close()
    print("\r[OK] Writing important attributes to file, sorted")


else:
    SnpChrom = [ Snp.info().chrom_id() for Snp in SnpList ]
    ChromList = sorted(list(set(SnpChrom)))
    ChromListLen = len(ChromList)
    Importance = []
    if Supervised:
        SnpSynthRfLabel = np.array([ Animal.phenotype() for Animal in PedList])
    for Chrom in ChromList:
        SubSnpList = [SnpList[Idx] for Idx,Chr in enumerate(SnpChrom) if Chr == Chrom]
        SubSnpIDs = [repr(Snp) for Snp in SubSnpList]
        SubTopSnps = min(TopSnps,len(SubSnpList))
        if NFeature != "auto":
            SubNFeature = min(NFeature, len(SubSnpIDs))
        else:
            SubNFeature = "auto"
        #--- Make Synthetic Data as mentioned by Leo Breiman for unsup. RF ------------#
        if not Supervised:
            print("\r[  ] Chromosome {}:{} [1/3]: Creating synthetic Data using SNPs for random forest".format(Chrom,ChromListLen),end="")
            SnpSynthRfData , SnpSynthRfLabel = BuildSnpUnsupRfData(SubSnpList), np.array([ Animal.phenotype() for Animal in PedList])
        else:
            SnpSynthRfData = np.array(SubSnpList).T
        #--- Perform unsupervised random forest ---------------------------------------#
        print("\r{}\r".format(" "*81),end="")
        print("\r[  ] Chromosome {}:{} [2/3]: Performing unsupervised random forest".format(Chrom,ChromListLen),end="")
        SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread, max_features = SubNFeature, criterion="gini")
        SnpRfSynthClf.fit(SnpSynthRfData, SnpSynthRfLabel)

        print("\r{}\r".format(" "*81),end="")
        print("\r[  ] Chromosome {}:{} [3/3] Calculating and storing attribute importance".format(Chrom,ChromListLen),end="")
        SubImportance = sorted( [ (ID,Chrom,Importance,SubSnpList[idx].conservation(),SubSnpList[idx].maf(),SubSnpList[idx].entropy()) for idx, (ID, Importance) in enumerate(zip(SubSnpIDs, SnpRfSynthClf.feature_importances_))], key=lambda x:x[2], reverse = True)
        TopImpSnps.extend([ SubSnpList[ SubSnpIDs.index(SnpID[0]) ] for SnpID in SubImportance[0:SubTopSnps]])
        Importance.extend(SubImportance)
        print("\r{}\r".format(" "*81),end="")
        print("\r[OK] Chromosome {}:{} unsupervised RF and attribute importance".format(Chrom,ChromListLen))
    print("\r[  ] Writing important attributes to file, sorted",end="")
    OutAttrFile = open(OutAttrFile, "w")
    OutAttrFile.write("ID\tChrom\tVarImp\tConserv\tMAF\tEntropy\n")
    for Tup in Importance:
        OutAttrFile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(Tup[0],Tup[1],Tup[2],Tup[3],Tup[4],Tup[5]))
    OutAttrFile.close()
    print("\r[OK] Writing important attributes to file, sorted")

TopImpSnps = sorted(TopImpSnps, key=lambda x:(x.info().chrom_id(),x.info().position()))
snp.write_mapped(TopImpSnps, PedList, OutFilePrefix, Verbose = True)
