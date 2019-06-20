import numpy as np
import argparse, operator, subprocess, re
from random import sample, choices, shuffle
from collections import Counter, OrderedDict
from operator import itemgetter 
from sklearn.ensemble import RandomForestClassifier
try:
    from snp import alphabet,animinfo,snpinfo,snp
except ImportError:
    print("Please put the 'snp.py'-Library in the same folder!")
    exit()

def BuildSnpUnsupRfData(SnpData):
    DataLength = len(SnpData[0])
    SynthData = np.append( np.array(SnpData).T, np.array([ np.random.choice(Row, size=len(Row),replace=True) for Row in np.array(SnpData) ]).T,axis=0)
    SynthLabel = np.repeat([0,1], [DataLength,DataLength],axis=0)
    return SynthData, SynthLabel

def BuildSnpUnsupRfDataOld(SnpData):
    DataLength = len(SnpData[0])
    SynthData = np.append( np.array(SnpData).T, np.array([ np.random.randint(min(Row),max(Row)+1,size=len(Row)) for Row in np.array(SnpData) ]).T,axis=0)
    SynthLabel = np.repeat([0,1], [DataLength,DataLength],axis=0)
    return SynthData, SynthLabel

def BuildSnpUnsupRfDataMehmet(SnpData):
    TmpData = []
    DataLength = len(SnpData[0])
    for Row in SnpList:
        Counts = Row.genotypecount()
        Sum = sum(Counts.values())
        Genotypes = OrderedDict(sorted({ Geno:Counts[Geno]/Sum for Geno in Counts if Counts[Geno] > 0}.items(), key = itemgetter(1), reverse = True))
        Genos = list(Genotypes.keys())
        if len(Genotypes) == 3:
            MAFs = list(map(lambda v: v * Row.maf(), [0.5,-0.1,-0.4]))
            #MAFs = list(map(lambda v: v * Row.maf(), [-0.35,-0.15,0.5]))
            Intervals = np.cumsum(np.array(list(Genotypes.values())) + np.array(MAFs))
            print(np.array(list(Genotypes.values())),MAFs,Intervals)
        else:
            MAFs = list(map(lambda v: v * Row.maf(), [-.3,.3]))
            Intervals = np.cumsum(np.array(list(Genotypes.values())) + np.array(MAFs))
            print(np.array(list(Genotypes.values())),MAFs,Intervals)
        NewData = []
        for round in range(DataLength):
            Item = np.random.uniform(0,1)
            for idx,Val in enumerate(Intervals):
                if Val > Item:
                    NewData.append(Row._GenotypesEncoder[Genos[idx]])
                    break
        TmpData.append(NewData)
    SynthData = np.append(np.array(SnpData).T,np.array(TmpData).T,axis=0)
    SynthLabel = np.repeat([0,1], [DataLength,DataLength],axis=0)
    return SynthData,SynthLabel

def BuildSnpUnsupRfShuffle(SnpData,Prob = [0.5,0.5]):
    DataLength = len(SnpData[0])
    SynthData = np.array(SnpData).T
    SynthLabel = np.random.choice([0,1],size=DataLength,replace=True,p=Prob)
    return SynthData, SynthLabel


Data = "data_out_fine"

SnpList, PedList = snp.read_mapped(Data, EncodeAdd = True, Verbose=True,Plink=False)
SnpIDs = [repr(Snp) for Snp in SnpList]

NTree = 500
NThread = 70


# print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
# #SnpSynthRfUnsupData , SnpSynthRfUnsupLabel = BuildSnpUnsupRfData(SnpList)
# SnpSynthRfUnsupData , SnpSynthRfUnsupLabel = BuildSnpUnsupRfShuffle(SnpList)
# SnpSynthRfSupData , SnpSynthRfSupLabel = np.array(SnpList).T, np.array([ Animal.phenotype() for Animal in PedList])
# print("\r[OK] Creating synthetic Data using SNPs for random forest")

# print("\r[  ] Running PLINK Analysis",end="")
# subprocess.Popen(["plink", "--file" , Data,"--assoc","--allow-no-sex", "--out","overlap","--cow"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
# PlinkFile = open("overlap.assoc","r")
# PlinkFile.readline()
# PlinkData =  [ (Entry[0],float(Entry[1])) for Entry in [ operator.itemgetter(2,9)(re.split('[\t\n ]+',Line.rstrip("\n"))) for Line in PlinkFile] ]
# """
# subprocess.Popen(["plink", "--file" , Data,"--model","--allow-no-sex", "--out","overlap","--cow"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
# PlinkFile = open("overlap.model")
# PlinkFile.readline()
# PlinkData =  [ (Entry[0],float(Entry[2])) for Entry in [ operator.itemgetter(2,5,10)(re.split('[\t\n ]+',Line.rstrip("\n"))) for Line in PlinkFile] if Entry[1] == "GENO" and Entry[2] != "NA" ]
# """

# PlinkData.sort(key=lambda k: (k[1]))
# PlinkFile.close()
# print("\r[OK] Running PLINK Analysis")


# Outfile = open("OverlapRfPlink.dat","w")
# Outfile.write("MTry\tOverlapUnsupRfPlink\tOverlapSupRfPlink\n")
# TopSnps = int(0.1*len(SnpList))
# TopPlinkSnps = [Item[0] for Item in PlinkData[0:TopSnps]]
# Start = round(0.1*len(SnpList))
# for NFeature in range(Start,len(SnpList)+1):
#     """
#     #print("NFeature = {}/{}".format(NFeature,len(SnpList)))
#     #print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
#     SnpSynthRfUnsupData , SnpSynthRfUnsupLabel = BuildSnpUnsupRfData(SnpList)
#     SnpSynthRfSupData , SnpSynthRfSupLabel = np.array(SnpList).T, np.array([ Animal.phenotype() for Animal in PedList])
#     #print("\r[OK] Creating synthetic Data using SNPs for random forest")

#     #print("\r[  ] Running PLINK Analysis",end="")
#     subprocess.Popen(["plink", "--file" , Data,"--assoc","--allow-no-sex", "--out","overlap"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
#     PlinkFile = open("overlap.assoc","r")
#     PlinkFile.readline()
#     PlinkData =  [ (Entry[0],float(Entry[1])) for Entry in [ operator.itemgetter(2,9)(re.split('[\t\n ]+',Line.rstrip("\n"))) for Line in PlinkFile] ]
#     PlinkData.sort(key=lambda k: (k[1]))
#     TopPlinkSnps = [Item[0] for Item in PlinkData[0:TopSnps]]
#     PlinkFile.close()
#     #print("\r[OK] Running PLINK Analysis")
#     """

#     #print("\r[  ] Performing unsupervised random forest",end="")
#     SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread, max_features = NFeature)
#     SnpRfSynthClf.fit(SnpSynthRfUnsupData, SnpSynthRfUnsupLabel)
#     #print("\r[OK] Performing unsupervised random forest")
#     Importance = sorted( [ (ID,importance,SnpList[idx].conservation(),SnpList[idx].maf(),SnpList[idx].entropy()) for idx,(ID,importance) in enumerate(zip(SnpIDs, SnpRfSynthClf.feature_importances_))], key=lambda x:x[1], reverse = True)
#     TopImpRfUnsupSnps = [ SnpID[0] for SnpID in Importance[0:TopSnps]] 

#     #print("\r[  ] Performing supervised random forest",end="")
    
#     SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread, max_features = NFeature)
#     SnpRfSynthClf.fit(SnpSynthRfSupData, SnpSynthRfSupLabel)
#     #print("\r[OK] Performing supervised random forest")
#     Importance = sorted( [ (ID,importance,SnpList[idx].conservation(),SnpList[idx].maf(),SnpList[idx].entropy()) for idx,(ID,importance) in enumerate(zip(SnpIDs, SnpRfSynthClf.feature_importances_))], key=lambda x:x[1], reverse = True)
#     TopImpRfSupSnps = [ SnpID[0] for SnpID in Importance[0:TopSnps]]

    

#     OverlapUnsup = len({Item for Item in TopPlinkSnps if Item in TopImpRfUnsupSnps})
#     OverlapSup = len({Item for Item in TopPlinkSnps if Item in TopImpRfSupSnps})
#     OverlapBoth = len({Item for Item in TopImpRfSupSnps if Item in TopImpRfUnsupSnps})


#     Outfile.write("{}\t{}\t{}\n".format(NFeature,OverlapUnsup/len(TopPlinkSnps),OverlapSup/len(TopPlinkSnps)))
#     print("{}\t{}\t{}\t\t{}".format(NFeature,OverlapUnsup,OverlapSup,OverlapBoth))
# Outfile.close()

Outfile = open("OverlapRfPlink_09-1FeautureCompare.dat","w")
Outfile.write("MTry\tOverlapUnsupRfPlink\tOverlapSupRfPlink\tPhenoPercent\n")
#Outfile.write("MTry\tOverlapUnsupRfPlink\tBreiman\tPhenoPercent\n")
for pheno in np.arange(0.5,1,0.05):
    PhenoOffset= round(pheno,2)
    PhenoProb = (PhenoOffset,round(1-PhenoOffset,2))
    
    print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
    #SnpSynthRfUnsupData , SnpSynthRfUnsupLabel = BuildSnpUnsupRfData(SnpList)
    SnpSynthRfUnsupData , SnpSynthRfUnsupLabel = BuildSnpUnsupRfShuffle(SnpList)
    #SnpSynthRfSupData , SnpSynthRfSupLabel = BuildSnpUnsupRfData(SnpList)
    SnpSynthRfSupData , SnpSynthRfSupLabel = np.array(SnpList).T, np.array([ Animal.phenotype() for Animal in PedList])
    print("\r[OK] Creating synthetic Data using SNPs for random forest")

    print("\r[  ] Running PLINK Analysis",end="")
    subprocess.Popen(["plink", "--file" , Data,"--assoc","--allow-no-sex", "--out","overlap","--cow"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
    PlinkFile = open("overlap.assoc","r")
    PlinkFile.readline()
    PlinkData =  [ (Entry[0],float(Entry[1])) for Entry in [ operator.itemgetter(2,9)(re.split('[\t\n ]+',Line.rstrip("\n"))) for Line in PlinkFile] ]
    """
    subprocess.Popen(["plink", "--file" , Data,"--model","--allow-no-sex", "--out","overlap","--cow"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
    PlinkFile = open("overlap.model")
    PlinkFile.readline()
    PlinkData =  [ (Entry[0],float(Entry[2])) for Entry in [ operator.itemgetter(2,5,10)(re.split('[\t\n ]+',Line.rstrip("\n"))) for Line in PlinkFile] if Entry[1] == "GENO" and Entry[2] != "NA" ]
    """

    PlinkData.sort(key=lambda k: (k[1]))
    PlinkFile.close()
    print("\r[OK] Running PLINK Analysis")


    
    TopSnps = int(0.1*len(SnpList))
    TopPlinkSnps = [Item[0] for Item in PlinkData[0:TopSnps]]
    Start = 1
    Start = round(0.9*len(SnpList))
    for NFeature in range(Start,len(SnpList)+1):
    #for NFeature in range(100):
        """
        #print("NFeature = {}/{}".format(NFeature,len(SnpList)))
        #print("\r[  ] Creating synthetic Data using SNPs for random forest",end="")
        SnpSynthRfUnsupData , SnpSynthRfUnsupLabel = BuildSnpUnsupRfData(SnpList)
        SnpSynthRfSupData , SnpSynthRfSupLabel = np.array(SnpList).T, np.array([ Animal.phenotype() for Animal in PedList])
        #print("\r[OK] Creating synthetic Data using SNPs for random forest")

        #print("\r[  ] Running PLINK Analysis",end="")
        subprocess.Popen(["plink", "--file" , Data,"--assoc","--allow-no-sex", "--out","overlap"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
        PlinkFile = open("overlap.assoc","r")
        PlinkFile.readline()
        PlinkData =  [ (Entry[0],float(Entry[1])) for Entry in [ operator.itemgetter(2,9)(re.split('[\t\n ]+',Line.rstrip("\n"))) for Line in PlinkFile] ]
        PlinkData.sort(key=lambda k: (k[1]))
        TopPlinkSnps = [Item[0] for Item in PlinkData[0:TopSnps]]
        PlinkFile.close()
        #print("\r[OK] Running PLINK Analysis")
        """
        #print("\r[  ] Performing unsupervised random forest",end="")
        TempImp = [1 for Idx in range(len(SnpList))]
        for idx in range(10):
            SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread, max_features = NFeature)
            SnpRfSynthClf.fit(SnpSynthRfUnsupData, SnpSynthRfUnsupLabel)
            #print("\r[OK] Performing unsupervised random forest")
            TempImp = [min(Last,Now) for (Last,Now) in (zip(TempImp, SnpRfSynthClf.feature_importances_))]
            
        Importance = sorted( [ (ID,importance,SnpList[idx].conservation(),SnpList[idx].maf(),SnpList[idx].entropy()) for idx,(ID,importance) in enumerate(zip(SnpIDs, TempImp))], key=lambda x:x[1], reverse = True)
        TopImpRfUnsupSnps = [ SnpID[0] for SnpID in Importance[0:TopSnps]] 


        #print("\r[  ] Performing supervised random forest",end="")
        TempImp = [1 for Idx in range(len(SnpList))]
        for idx in range(10):            
            SnpRfSynthClf = RandomForestClassifier(n_estimators = NTree, n_jobs = NThread, max_features = NFeature)
            SnpRfSynthClf.fit(SnpSynthRfSupData, SnpSynthRfSupLabel)
            #print("\r[OK] Performing supervised random forest")
        TempImp = [min(Last,Now) for (Last,Now) in (zip(TempImp, SnpRfSynthClf.feature_importances_))]
        Importance = sorted( [ (ID,importance,SnpList[idx].conservation(),SnpList[idx].maf(),SnpList[idx].entropy()) for idx,(ID,importance) in enumerate(zip(SnpIDs, TempImp))], key=lambda x:x[1], reverse = True)
        TopImpRfSupSnps = [ SnpID[0] for SnpID in Importance[0:TopSnps]]

        

        OverlapUnsup = len({Item for Item in TopPlinkSnps if Item in TopImpRfUnsupSnps})
        OverlapSup = len({Item for Item in TopPlinkSnps if Item in TopImpRfSupSnps})
        OverlapBoth = len({Item for Item in TopImpRfSupSnps if Item in TopImpRfUnsupSnps})


        Outfile.write("{}\t{}\t{}\t{}\n".format(NFeature,OverlapUnsup/len(TopPlinkSnps),OverlapSup/len(TopPlinkSnps),PhenoOffset))
        print("{}\t{}\t{}\t\t{}\t\t{}".format(NFeature,OverlapUnsup,OverlapSup,OverlapBoth,PhenoOffset))
Outfile.close()