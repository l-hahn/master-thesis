#!/usr/bin/python3
import argparse, re
try:
    from snp import alphabet,animinfo,snpinfo,snp
except ImportError:
    print("Please put the 'snp.py'-Library in the same folder!")
    exit()
from random import choice

OutMapFileNameCopy = "data.map"
OutPedFileNameCopy = "data.ped"

RequireSnp = False

MapList = []
PedList = []
ChosenIndexListCopy = set()

MapLineCount = 0
PedLineCount = 0


Parser = argparse.ArgumentParser(description='Sample-Creator for a MAP and PED file from given data. Provide additional file if some SNPs have to occure in the final output.')
Parser.add_argument('--map', metavar='{MAP-File-Name}', type=str, nargs=1, help='Input MAP file, where a subset is taken from', required = True)
Parser.add_argument('--ped', metavar='{PED-File-Name}', type=str, nargs=1, help='Input PED file, where a subset is taken from', required = True)
Parser.add_argument('-n', metavar='<Sample-Set-Sizes>', type=int, nargs="+", help='Final sample-set-size; the amount of SNPs to be choosen.', required = True)
Parser.add_argument('--snp', metavar='{SNP-File-Name}', type=str, nargs=1, help='SNP file, where in the first column SNP-IDS are listed that have to be added to the sample set', required = False)
Parser.add_argument('--omp', metavar='{MAP-Out-Name}', type=str, nargs=1, help='Output-MAP file name', required = False)
Parser.add_argument('--opd', metavar='{PED-Out-Name}', type=str, nargs=1, help='Output-PED file name', required = False)
Parser.add_argument('-b', action="store_true", help='Print integer instead if letters (A,C,G,T)', required = False)
Parser.add_argument('-g', action="store_true", help='Print genotypes instead of alleles', required = False)
Args = Parser.parse_args()
if Args.omp:
    OutMapFileNameCopy = Args.omp[0]
if Args.opd:
    OutPedFileNameCopy = Args.opd[0]
if Args.snp:
    RequireSnp = True
    SnpFileName = Args.snp
SnpSizes = Args.n
Bin = Args.b
Gen = Args.g


InMapFile = open(Args.map[0],"r")
print("\r[  ] Processing MAP-File",end="")
for idx,Line in enumerate(InMapFile):
    #print("\r[  ] Processing MAP-File: line {}".format(idx),end="")
    Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
    if len(Entries) != 4:
        print("!!WARNING!! Line {} in MAP-File has to many entries!\n".format(idx)+Line)
    MapList.append(snpinfo(Entries[0],Entries[1],Entries[2],Entries[3]))
    MapLineCount += 1
print("\r[OK] Processing MAP-File")
InMapFile.close()

SnpList = [snp(SnpInfo) for SnpInfo in MapList]
del MapList

InPedFile = open(Args.ped[0],"r")
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


if RequireSnp:
    SnpsNeeded = []
    InSnpFile = open(Args.snp[0],"r")
    for Line in InSnpFile:
        Entries = Line.rstrip("\n").split("\t")
        SnpsNeeded.append(Entries[0])
    ChosenIndexListCopy = {(SnpList.index(snp(snpinfo(0,MapInfo,0,0))),MapInfo) for MapInfo in SnpsNeeded if any(Snp.info().snp_id() == MapInfo for Snp in SnpList)}
    for SetItem in ChosenIndexListCopy:
        if SetItem[1] not in SnpsNeeded:
            print("!!WARNING!! Snp {} from required file cannot be found in MAP-File!".format(SetItem[1]))
    ChosenIndexListCopy = {SetItem[0] for SetItem in ChosenIndexListCopy}


for SnpSize in SnpSizes:
    ChosenIndexList = set(ChosenIndexListCopy)
    OutMapFileName = OutMapFileNameCopy
    OutPedFileName = OutPedFileNameCopy

    if len(SnpSizes) > 1:
        print("")

    if len(ChosenIndexList) > SnpSize:
        SnpSize = len(ChosenIndexList)
        print("!!WARNING!! Required Snps from file {} exceeds the number of sample snps! Resetting n = {}".format(Args.snp[0],SnpSize))

    print("\r[  ] Selecting Snps for sample set (size: {})".format(SnpSize),end="")
    while len(ChosenIndexList) < SnpSize:
        ChosenIndexList.add(choice(range(0,len(SnpList))))
    NewSnpList = [SnpList[idx] for idx in sorted(ChosenIndexList)]
    print("\r[OK] Selecting Snps for sample set (size: {})".format(SnpSize))

    if len(SnpSizes) > 1:
        if OutMapFileName.rfind("/") <= OutMapFileName.rfind("."):
            OutMapFileName = OutMapFileName[:OutMapFileName.rfind(".")]+"_{}-snps.map".format(len(ChosenIndexList))
        else:
            OutMapFileName = OutMapFileName+"_{}-snps.map".format(len(ChosenIndexList))
        if OutPedFileName.rfind("/") <= OutPedFileName.rfind("."):
            OutPedFileName = OutPedFileName[:OutPedFileName.rfind(".")]+"_{}-snps.ped".format(len(ChosenIndexList))
        else:
            OutPedFileName = OutPedFileName+"_{}-snps.ped".format(len(ChosenIndexList))

    print("\r[  ] Writing Output-MAP-File (target: {})".format(OutMapFileName),end="")
    OutMapFile = open(OutMapFileName,"w")
    for Snp in NewSnpList:
        OutMapFile.write(str(Snp.info())+"\n")
    OutMapFile.close()
    print("\r[OK] Writing Output-MAP-File (target: {})".format(OutMapFileName))

    print("\r[  ] Writing Output-PED-File (target: {})".format(OutPedFileName),end="")
    OutPedFile = open(OutPedFileName,"w")
    for idx,AnimInfo in enumerate(PedList):
        OutPedFile.write(str(AnimInfo))
        for Snp in NewSnpList:
            if not Bin and not Gen:
                DiNu = alphabet.decode(Snp[idx])
                OutPedFile.write(" {} {}".format(DiNu[0],DiNu[1]))
            elif Bin and not Gen:
                DiNu = alphabet.decode(Snp[idx])
                OutPedFile.write(" {} {}".format(alphabet.enc(DiNu[0]),alphabet.enc(DiNu[1])))
            elif not Bin and Gen:
                DiNu = alphabet.decode(Snp[idx])
                OutPedFile.write(" {}".format(DiNu))
            else:
                OutPedFile.write(" {}".format(Snp[idx]))
        OutPedFile.write("\n")
    OutPedFile.close()
    print("\r[OK] Writing Output-PED-file (target: {})".format(OutPedFileName))