import math, re
from collections import Counter
from itertools import combinations_with_replacement as combinations

class slice():
    def __init__(self, ProPortion=0, Item="-"):
        self._ProPortion = ProPortion
        self._Item = Item
    def __str__(self):
        return "["+str(self._ProPortion)+":"+str(self._Item)+"]"
    def __repr__(self):
        return repr(self._Item)

    def proportion(self):
        return self._ProPortion
    def item(self):
        return self._Item
    def set_proportion(self, ProPortion):
        self._ProPortion = ProPortion
    def set_item(self, Item):
        self._Item = Item


class alphabet():
    _AlphaSize = 25;
    _AlphaSymSize = 15;
    _AlphaSym = True;
    _EncodeNu = {"A":0,"C":1,"G":2,"T":3,"0":4}
    _DecodeNu = {0:"A",1:"C",2:"G",3:"T",4:"0"}
    _EncodeDiNu = {"AA":0, "AC":1, "AG":2, "AT":3, "A0":4, 
                    "CA":5, "CC":6, "CG":7, "CT":8, "C0":9,
                    "GA":10, "GC":11, "GG":12, "GT":13, "G0":14,
                    "TA":15, "TC":16, "TG":17, "TT":18, "T0":19,
                    "0A":20, "0C":21, "0G":22, "0T":23, "00":24}
    _EncodeDiNuSym = {"AA":0, "AC":1, "AG":2, "AT":3, "A0":4, 
                        "CA":1, "CC":5, "CG":6, "CT":7, "C0":8,
                        "GA":2, "GC":6, "GG":9, "GT":10, "G0":11,
                        "TA":3, "TC":7, "TG":10, "TT":12, "T0":13,
                        "0A":4, "0C":8, "0G":11, "0T":13, "00":14}
    
    _DecodeDiNu = {0:"AA", 1:"AC", 2:"AG", 3:"AT", 4:"A0", 
                    5:"CA", 6:"CC", 7:"CG", 8:"CT", 9:"C0",
                    10:"GA", 11:"GC", 12:"GG", 13:"GT", 14:"G0",
                    15:"TA", 16:"TC", 17:"TG", 18:"TT", 19:"T0",
                    20:"0A", 21:"0C", 22:"0G", 23:"0T", 24:"00"}
    _DecodeDiNuSym = {0:"AA", 1:"AC", 2:"AG", 3:"AT", 4:"A0", 
                    5:"CC", 6:"CG", 7:"CT", 8:"C0", 9:"GG",
                    10:"GT", 11:"G0", 12:"TT", 13:"T0",14:"00"}
    def symmetric():
        alphabet._AlphaSym = True;
    def nosymmetric():
        alphabet._AlphaSym = False;
    def size():
        if alphabet._AlphaSym:
            return alphabet._AlphaSymSize
        return alphabet._AlphaSize
    def encode(DiNu):
        if alphabet._AlphaSym:
            return alphabet._EncodeDiNuSym[DiNu]
        return alphabet._EncodeDiNu[DiNu]
    def decode(DiNu):
        if alphabet._AlphaSym:
            return alphabet._DecodeDiNuSym[DiNu]
        return alphabet._DecodeDiNu[DiNu]
    def enc(Nu):
        return alphabet._EncodeNu[Nu]
    def dec(Nu):
        return alphabet._DecodeNu[Nu]


class animinfo():
    def __init__(self, FamID, SampID, PatID, MatID, Sex, Pheno):
        self._Family = FamID
        self._Sample = SampID
        self._Father = PatID
        self._Mother = MatID
        self._Sexual = Sex
        self._Affect = Pheno
    def __str__(self,Delimiter=" "):
        if type(Delimiter) != str:
            Delimiter = str(Delimiter)
        return str(self._Family) + Delimiter + str(self._Sample) + Delimiter \
                + str(self._Father) + Delimiter + str(self._Mother) + Delimiter \
                + str(self._Sexual) + Delimiter + str(self._Affect)
    def __repr__(self):
        return repr(self._Sample)

    def family(self):
        return self._Family
    def sample(self):
        return self._Sample
    def father(self):
        return self._Father
    def mother(self):
        return self._Mother
    def sex(self):
        return self._Sexual
    def phenotype(self):
        return self._Affect

    def set_family(self, FamID):
        self._Family = FamID
    def set_sample(self, SampID):
        self._Sample = SampID
    def set_father(self, PatID):
        self._Father = PatID
    def set_mother(self, MatID):
        self._Mother = MatID
    def set_sex(self, Sex):
        self._Sexual = Sex
    def set_phenotype(self, Pheno):
        self._Affect = Pheno


class snpinfo():
    def __init__(self, ChrmID, SnpID, GenDist, AbsPos):
        self._ChromosomeID = ChrmID
        self._SNPID = SnpID
        self._Distance = GenDist
        self._Position = AbsPos
    def __str__(self, Delimiter=" "):
        if type(Delimiter) != str:
            Delimiter = str(Delimiter)
        return str(self._ChromosomeID) + Delimiter + str(self._SNPID) \
                + Delimiter + str(self._Distance) + str(Delimiter) \
                + str(self._Position)
    def __repr__(self):
        return repr(self._SNPID)

    def __eq__(self, other):
        return self._SNPID == other._SNPID
    def __lt__(self, other):
        return self._Position < other._Position
    def __le__(self, other):
        return self._Position <= other._Position
    def __ge__(self, other):
        return self._Position > other._Position
    def __gt__(self, other):
        return self._Position >= other._Position

    def chrom_id(self):
        return self._ChromosomeID
    def snp_id(self):
        return self._SNPID
    def distance(self):
        return self._Distance
    def position(self):
        return self._Position

    def set_chrom_id(self, ChrmID):
        self._ChromosomeID = ChrmID
    def set_snp_id(self, SnpID):
        self._SNPID = SnpID
    def set_distance(self, GenDist):
        self._Distance = GenDist
    def set_position(self, AbsPos):
        self._Position = AbsPos


class snp():
    def __init__(self, SnpInfo):
        if type(SnpInfo) == snpinfo:
            self._SnpInfo = SnpInfo
        else: 
            self._SnpInfo = snpinfo(0,SnpInfo,0,0)
        self._DiNuList = []
        self._EncodeAdd = False
    def __str__(self):
        return str(self._SnpInfo)
    def __repr__(self):
        return repr(self._SnpInfo)
    def __len__(self):
        return len(self._DiNuList)
    def __getitem__(self,Idx):
        return self._DiNuList[Idx]
    def __setitem__(self,Idx,DiNu):
        self._DiNuList[Idx] = DiNu

    def __eq__(self, other):
        return self._SnpInfo == other._SnpInfo
    def __lt__(self, other):
        return self._SnpInfo < other._SnpInfo
    def __le__(self, other):
        return self._SnpInfo <= other._SnpInfo
    def __ge__(self, other):
        return self._SnpInfo > other._SnpInfo
    def __gt__(self, other):
        return self._SnpInfo >= other._SnpInfo

    def info(self):
        return self._SnpInfo

    def append(self, DiNu):
        self._DiNuList.append(DiNu)
    def extend(self, DiNuList):
        self._DiNuList.extend(DiNuList)
    def insert(self, Idx, DiNu):
        self._DiNuList.insert(Idx,DiNu)
    def remove(self, DiNu):
        self._DiNuList.remove(DiNu)
    def pop(self, DiNuList=[]):
        self._DiNuList.pop(DiNuList)
    def clear(self):
        self._DiNuList.clear()
    def count(self, DiNu):
        self._DiNuList.count(DiNu)
    def conservation(self):
        return max(Counter(self._DiNuList).values())

    def encode_additive(self):
        if not self._EncodeAdd:
            self._EncodeAdd = True
            Alleles = {'A':0,'C':0,'G':0,'T':0, '0':0}
            GenotypesCount = Counter(self._DiNuList)
            for Genotyp in GenotypesCount:
                GenoStr = alphabet.decode(Genotyp)
                Alleles[GenoStr[0]] += GenotypesCount[Genotyp]
                Alleles[GenoStr[1]] += GenotypesCount[Genotyp]
            self._Alleles = [ Counts[0] for Counts in sorted(Alleles.items(), key=lambda kv:kv[1], reverse=True) if Counts[1] > 0 ]
            Genotypes = ["".join(Genotype) for Genotype in combinations(self._Alleles,2)]
            GenotypesEncoder = {alphabet.encode(Genotype):idx for idx,Genotype in enumerate(Genotypes)}
            GenotypesEncoder = {**GenotypesEncoder, **{ alphabet.encode(Genotype[::-1]):GenotypesEncoder[alphabet.encode(Genotype)] for Genotype in Genotypes }}
            self._DiNuList = [ GenotypesEncoder[Genotype] for Genotype in self._DiNuList ]

    def read_mapped(InMap, InPed, EncodeAdd = False , Verbose = False):
        MapList = []
        PedList = []
        MapLineCount = 0
        PedLineCount = 0

        InMapFile = open(InMap,"r")
        if Verbose: print("\r[  ] Processing MAP-File",end="")
        for idx,Line in enumerate(InMapFile):
            if Verbose: print("\r[  ] Processing MAP-File: line {}".format(idx),end="")
            Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
            if len(Entries) != 4 and Verbose:
                print("!!WARNING!! Line {} in MAP-File has to many entries!\n".format(idx)+Line)
            MapList.append(snpinfo(int(Entries[0]),Entries[1],int(Entries[2]),int(Entries[3])))
            MapLineCount += 1
        if Verbose: print("\r{}\r[OK] Processing MAP-File".format(" "*80))
        InMapFile.close()

        SnpList = [snp(SnpInfo) for SnpInfo in MapList]
        del MapList

        InPedFile = open(InPed,"r")
        if Verbose: print("\r[  ] Processing PED-File",end="")
        for adx,Line in enumerate(InPedFile):
            if Verbose: print("\r[  ] Processing PED-File: line {}".format(adx),end="")
            Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
            if (len(Entries)-6)/2 != MapLineCount and Verbose:
                print("!!WARNING!! Line {} has a different amount({}) of snps than required({})!".format(adx,(len(Entries)-6)/2,MapLineCount))
                continue
            PedList.append(animinfo(Entries[0],Entries[1],Entries[2],Entries[3],Entries[4],Entries[5]))
            for idx,jdx in enumerate(range(6,len(Entries),2)):
                SnpList[idx].append(alphabet.encode("{}{}".format(Entries[jdx],Entries[jdx+1])))
            PedLineCount += 1 
        if Verbose: print("\r{}\r[OK] Processing PED-File".format(" "*80))
        InPedFile.close()
        if EncodeAdd:
            for idx in range(len(SnpList)):
                SnpList[idx].encode_additive()
        return SnpList, PedList

    def write_mapped(SnpList, PedList, OutMap, OutPed, Verbose = False):
        if Verbose: print("\r[  ] Writing Output-MAP-File (target: {})".format(OutMap),end="")
        OutMapFile = open(OutMap,"w")
        for Snp in SnpList:
            OutMapFile.write(str(Snp.info())+"\n")
        OutMapFile.close()
        if Verbose: print("\r[OK] Writing Output-MAP-File (target: {})".format(OutMap))

        if Verbose: print("\r[  ] Writing Output-PED-File (target: {})".format(OutPed),end="")
        OutPedFile = open(OutPed,"w")
        for idx,AnimInfo in enumerate(PedList):
            OutPedFile.write(str(AnimInfo))
            for Snp in SnpList:
                DiNu = alphabet.decode(Snp[idx])
                OutPedFile.write(" {} {}".format(DiNu[0],DiNu[1]))
            OutPedFile.write("\n")
        OutPedFile.close()
        if Verbose: print("\r[OK] Writing Output-PED-file (target: {})".format(OutPed))