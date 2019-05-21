import argparse, re

try:
    from snp import alphabet,animinfo,snpinfo,snp
except ImportError:
    print("Please put the 'snp.py'-Library in the same folder!")
    exit()
from random import choice

OutFileName = "data"
RequireSnp = False
ChosenIndexListCopy = set()


Parser = argparse.ArgumentParser(description='Sample-Creator for a MAP and PED file from given data. Provide additional file if some SNPs have to occure in the final output.')
Parser.add_argument('--file', metavar='{File-Name-Prefix}', type=str, nargs=1, help='File prefix for input MAP and PED file, where a subset is taken from', required = True)
Parser.add_argument('-n', metavar='<Sample-Set-Sizes>', type=int, nargs="+", help='Final sample-set-size; the amount of SNPs to be choosen.', required = True)
Parser.add_argument('--snp', metavar='{SNP-File-Name}', type=str, nargs=1, help='SNP file, where in the first column SNP-IDS are listed that have to be added to the sample set', required = False)
Parser.add_argument('--out', metavar='{MAP-Out-Name}', type=str, nargs=1, help='File prefix for output MAP and PED file', required = False)
Parser.add_argument('-b', action="store_true", help='Print integer instead if letters (A,C,G,T)', required = False)
Parser.add_argument('-g', action="store_true", help='Print genotypes instead of alleles', required = False)
Args = Parser.parse_args()
if Args.out:
    OutFileName = Args.out[0]
if Args.snp:
    RequireSnp = True
    SnpFileName = Args.snp
SnpSizes = Args.n
Bin = Args.b
Gen = Args.g

SnpList, FamList = snp.read_mapped(Args.file[0],Verbose=True)



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
    OutFileSizeName = OutFileName

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
        OutFileSizeName = OutFileSizeName+"_{}-snps".format(len(ChosenIndexList))
    snp.write_mapped(NewSnpList,FamList,OutFileSizeName,Verbose=True)
    