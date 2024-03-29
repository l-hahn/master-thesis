import math, re, os, subprocess
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
        return str(self._SNPID)

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
        self._AlleleCounted = False
        self._GenotypesCounted = False
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
    def is_additive_encoded(self):
        return self._EncodeAdd

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
        if not self._GenotypesCounted:
            self.count_genotype()
        return max(self._GenotypesCount.values())/len(self)
    def maf(self):
        if not self._AlleleCounted:
            self.count_allel()
        return min([self._AlleleCount[Allele] for Allele in self._Alleles])/(sum(self._AlleleCount.values()))
    def entropy(self):
        return sum([ -Geno/len(self)*math.log2(Geno/len(self)) for Geno in self._GenotypesCount.values() ])
    def count_allel(self):
        ReDo = False
        if self._EncodeAdd:
            self.decode_additive()
            ReDo = True
        if not self._GenotypesCounted:
            self.count_genotype()
        Alleles = {'A':0,'C':0,'G':0,'T':0, '0':0}
        for Genotyp in self._GenotypesCount:
            GenoStr = alphabet.decode(Genotyp)
            Alleles[GenoStr[0]] += self._GenotypesCount[Genotyp]
            Alleles[GenoStr[1]] += self._GenotypesCount[Genotyp]
        self._AlleleCount = Alleles
        self._Alleles = [ Counts[0] for Counts in sorted(Alleles.items(), key=lambda kv:kv[1], reverse=True) if Counts[1] > 0 and Counts[0] != '0']
        self._AlleleCounted = True
        if ReDo:
            self.encode_additive()
        return self._AlleleCount
    def allelcount(self):
        return self._AlleleCount
    def count_genotype(self):
        self._GenotypesCount = Counter(self._DiNuList)
        self._Genotypes = [ Counts[0] for Counts in sorted(self._GenotypesCount.items(), key=lambda kv:kv[1], reverse=True) ]
        self._GenotypesCounted = True
        return self._GenotypesCount
    def genotypecount(self):
        return self._GenotypesCount

    def has_plink():
        try:
            subprocess.Popen(["plink"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True).communicate()
        except OSError as Except:
            if Except.errno == os.errno.ENOENT:
                return False
        return True

    def encode_additive(self):
        if not self._EncodeAdd:
            if not self._AlleleCounted:
                self.count_allel()
            Genotypes = ["".join(Genotype) for Genotype in combinations(self._Alleles,2)]
            GenotypesEncoder = {alphabet.encode(Genotype):(idx) for idx,Genotype in enumerate(Genotypes)}
            if self._AlleleCount['0'] != 0:
                GapEncoder = { alphabet.encode(Genotype):-1 for Genotype in [ Allele+"0" for Allele in self._Alleles+['0']]}
                GenotypesEncoder = {**GenotypesEncoder, **GapEncoder}
            self._GenotypesEncoder = {**GenotypesEncoder, **{ alphabet.encode(Genotype[::-1]):GenotypesEncoder[alphabet.encode(Genotype)] for Genotype in Genotypes }}
            self._DiNuList = [ self._GenotypesEncoder[Genotype] for Genotype in self._DiNuList ]
            self._EncodeAdd = True
    def decode_additive(self):
        if self._EncodeAdd:
            GenotypesDecoder = {Enc:Geno for Geno,Enc in self._GenotypesEncoder.items()}
            self._DiNuList = [ GenotypesDecoder[Genotype] for Genotype in self._DiNuList ]
            self._EncodeAdd = False

    def read_mapped(InFile, EncodeAdd = False , Verbose = False, Plink = True, TPed = False):
        MapList = []
        FamList = []
        MapLineCount = 0
        PedLineCount = 0
        
        if (snp.has_plink() and Plink) or TPed == True:
            if not TPed:
                if Verbose: print("\r[  ] Running PLINK for data transposition",end="")
                InMapFile = open(InFile+".map","r")
                Chrs = max([int(re.split('[\t\n ]+',Line.rstrip("\n"))[0]) for Line in InMapFile])
                InMapFile.close()
                subprocess.Popen(["plink", "--file" , InFile,"--chr-set", str(Chrs), "--recode", "transpose","--allow-no-sex", "--out", InFile], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()            
                if Verbose: print("\r[OK] Running PLINK for data transposition")

            if Verbose: print("\r[  ] Processing TFAM-File",end="")
            InTfamFile = open(InFile+".tfam","r")
            for adx,Line in enumerate(InTfamFile):
                if Verbose: print("\r[  ] Processing TFAM-File: line {}".format(adx),end="")
                Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
                if len(Entries) != 6 and Verbose:
                    print(" !!WARNING!! Line {} in TFAM-File has wrong entry number!\n".format(idx)+Line)
                FamList.append(animinfo(Entries[0],Entries[1],Entries[2],Entries[3],Entries[4],Entries[5]))
                PedLineCount += 1
            InTfamFile.close()
            if Verbose: print("\r{}\r[OK] Processing TFAM-File".format(" "*80))

            if Verbose: print("\r[  ] Processing TPED-File",end="")
            InTpedFile = open(InFile+".tped","r")
            SnpList = []
            for adx,Line in enumerate(InTpedFile):
                #if Verbose: print("\r[  ] Processing TPED-File: line {}".format(adx),end="")
                Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
                if (len(Entries)-4)/2 != PedLineCount and Verbose:
                    print(" !!WARNING!! Line {} has a different amount({}) of genotpyes than required({})!".format(adx,(len(Entries)-4)/2,PedLineCount))
                    continue
                SnpList.append(snp(snpinfo(int(Entries[0]),Entries[1],int(Entries[2]),int(Entries[3]))))
                SnpList[adx].extend( [ alphabet.encode("{}{}".format(Entries[idx],Entries[idx+1])) for idx in range(4,len(Entries),2) ] )
            InTpedFile.close()
            if Verbose: print("\r{}\r[OK] Processing TPED-File".format(" "*80))

            if EncodeAdd:
                for idx in range(len(SnpList)):
                    SnpList[idx].encode_additive()
            return SnpList, FamList

        else:
            if Verbose: print("\r[NO] Running PLINK for data transposition")
            if Verbose: print("\r[  ] Processing MAP-File",end="")
            InMapFile = open(InFile+".map","r")
            for idx,Line in enumerate(InMapFile):
                if Verbose: print("\r[  ] Processing MAP-File: line {}".format(idx),end="")
                Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
                if len(Entries) != 4 and Verbose:
                    print(" !!WARNING!! Line {} in MAP-File has to many entries!\n".format(idx)+Line)
                MapList.append(snpinfo(int(Entries[0]),Entries[1],int(Entries[2]),int(Entries[3])))
                MapLineCount += 1
            InMapFile.close()
            if Verbose: print("\r{}\r[OK] Processing MAP-File".format(" "*80))

            SnpList = [snp(SnpInfo) for SnpInfo in MapList]
            del MapList

            if Verbose: print("\r[  ] Processing PED-File",end="")
            InPedFile = open(InFile+".ped","r")
            for adx,Line in enumerate(InPedFile):
                if Verbose: print("\r[  ] Processing PED-File: line {}".format(adx),end="")
                Entries = re.split('[\t\n ]+',Line.rstrip("\n"))
                if (len(Entries)-6)/2 != MapLineCount and Verbose:
                    print(" !!WARNING!! Line {} has a different amount({}) of snps than required({})!".format(adx,(len(Entries)-6)/2,MapLineCount))
                    continue
                FamList.append(animinfo(Entries[0],Entries[1],Entries[2],Entries[3],Entries[4],Entries[5]))
                for idx,jdx in enumerate(range(6,len(Entries),2)):
                    SnpList[idx].append(alphabet.encode("{}{}".format(Entries[jdx],Entries[jdx+1])))
                PedLineCount += 1 
            if Verbose: print("\r{}\r[OK] Processing PED-File".format(" "*80))
            InPedFile.close()
            if EncodeAdd:
                for idx in range(len(SnpList)):
                    SnpList[idx].encode_additive()
            return SnpList, FamList

    def write_mapped(SnpList, FamList, OutFile, Verbose = False):
        if Verbose: print("\r[  ] Writing Output-MAP-File (target: {})".format(OutFile+".map"),end="")
        OutMapFile = open(OutFile+".map","w")
        for Snp in SnpList:
            OutMapFile.write(str(Snp.info())+"\n")
        OutMapFile.close()
        if Verbose: print("\r[OK] Writing Output-MAP-File (target: {})".format(OutFile+".map"))

        if Verbose: print("\r[  ] Writing Output-PED-File (target: {})".format(OutFile+".ped"),end="")
        OutPedFile = open(OutFile+".ped","w")
        for Snp in SnpList:
            if Snp.is_additive_encoded():
                Snp.decode_additive()
        for idx,AnimInfo in enumerate(FamList):
            OutPedFile.write(str(AnimInfo))
            for Snp in SnpList:
                DiNu = alphabet.decode(Snp[idx])
                OutPedFile.write(" {} {}".format(DiNu[0],DiNu[1]))
            OutPedFile.write("\n")
        OutPedFile.close()
        if Verbose: print("\r[OK] Writing Output-PED-file (target: {})".format(OutFile+".ped"))
        
        """
        For more output:
            for Snp in SnpList:
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
        """