# This program will calculate the melting temperature of primers
# The equation used will be: Tm=(DeltaH/(Deltas+R*ln(C/4))-273.15+16.6log[salt]
# U need a file in text format of in the order: name forwardprimer reverseprimer

# Data of DeltaH, DeltaS and DeltaG was extracted from: Breslauer, K. J., R. Frank, et al.  (1986). 

from math import log, log10
import AUC, PrimerDimer


DeltaG = {'AA':1.9, 'AC':1.3, 'GT':1.3, 'AG':1.6,'CC':3.1,
'TT':1.9, 'CG':3.6, 'TC':1.6,'GG':3.1, 'GC':3.1, 'AT':1.5, 'GA':1.6,
'TG':1.9, 'TA':0.9, 'CA':1.9, 'CT':1.6}

DeltaH = {'AA':9.1, 'AT':8.6, 'AC':6.5, 'AG':7.8, 'TA':6.0, 'TT':9.1,
'TC':5.6, 'TG':5.8, 'CA':5.8, 'CT':7.8, 'CC':11.0, 'CG':11.9, 'GA':5.6,
'GT':6.5, 'GC':11.1, 'GG':11.0}

DeltaS = {'AA':24, 'AT':23.9, 'AC':17.3, 'AG':20.8, 'TA':16.9, 'TT':24.0,
'TC':13.5, 'TG':12.9,'CA':12.9, 'CT':20.8, 'CC':26.6, 'CG':27.8, 'GA':13.5,
'GT':17.3, 'GC':26.7, 'GG':26.6}

RevDNA={'A':'T','T':'A','G':'C','C':'G','N':'N'}

#list of all codon combinations
Bases = ['T', 'C', 'A', 'G']
Codons = [a+b+c for a in Bases for b in Bases for c in Bases]

#amino-acid sequence translation
AAdict = {"TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser",
              "TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp",
              "TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu",
              "CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro",
              "CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg",
              "CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met",
              "ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn",
              "AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg",
              "GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala",
              "GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu",
              "GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}

"""
Atributtes of the class Primer:
	Conc -	Conc is the primer concentration in uM It is automatically  set to 0.2uM when the class is created.
	SaltConc - SaltConc is the salt concentration in mM. It is automatically  set to 50mM.
	Type -	Type is FP for forward primer and RP for reverse primer.
Methods of the class Primer
	SetConc - use this method to change the concetration of the primer. Submit concentration in uM units.
	SetSaltConc -	use this method to change the salt concentration. submit in nm
	Tm -	this methods returns the melting temperature of the primer.
	PerfectMatch - Returns the part of the sequence that binds to the template. At least 10 bases long.
	ReverseMatch -	Returns the part of the sequence that binds when the primer is of type RP.
	RestOfPrimer -	This returns the overhang of the primer. it takes the result of MaxMatch as input.
	MatchAndRest - This method combines both the PerfectMatch and the RestOfPrimer. It returns a tupple (PerfectMatch, RestOfPrimer)
	OligoDeltaG - Calculates the DeltaG of the oligonucleotide (or primer).
	DeltaGWindow - A sliding window that runs through the sequence, for each window it calculates the DeltaG and returns a list of values.

Attributes of the class PcrMachine:
	FP -	the forward primer.
	RP -	the reverse primer.
	Template - the template you will use in the PCR.
	Tm - melting temperature of the primer.
	GC - GC content of the primer (percent)
	Last5- number of G+C in the last 5 nucleotides.
	
Methods of the class PcrMachine:
	RawPCR -	this method returns the PCR product but it allows no mismatches. It is used once you have the parts
			of the primers that you know match perfectly (i.e. no hanging sequence allowed).
	MakePCR -	this is the fully active method. basically its the only one you need. It will check the primers,
			figure out the hangover and return the full PCR product.
"""

# Complementary reverse function. return the complementary sequence of any DNA strand.
# DNA must be uppercase.
def CompRev(Sequence):
	NewSeq=''
	for i in Sequence:
		NewSeq=RevDNA[i]+NewSeq
	return NewSeq

def SeqToAA(seq):
	""""
	Translates gene sequence, including numbers
	"""
	trans_output = ''
	for i in xrange(0, len(seq),3):
		#reference the dicitonary and translate codon
		codon = seq[i: i+3]
		#if codon is incomplete, mark it as such; it will be rendered as stop codon
		if len(codon)!=3: 
			codon = 'incomplete' #replace with incomplete
			a_acid = 'XXX'
		else:
			a_acid = AAdict[codon][2:len(AAdict[codon])]
			#AAcount = trans_output.count(AAdict[codon][2:len(AAdict[codon])])+1
			AAcount = i/3+1
			trans_output += a_acid+str(AAcount)+'_'		
	return(trans_output)
	
	

# return the number of G's and C's in a DNA sequence.
# DNA must be uppercase.
def GCcount(Seq):
	return float(Seq.count('G') + Seq.count('C'))

# return the percentace of GC in a DNA sequence
# DNA must be uppercase.
def GCcontent (Seq):
	GC = Seq.count('G')+Seq.count('C')
	return (float(GC)/len(Seq))*100

def OligoDeltaG(Seq):
	Sum = 0
	for i in range(len(Seq)-1):
		pair = Seq[i:i+2]
		Sum+=DeltaG[pair]
	return Sum
	
def DeltaGWindow(Seq, Window=5):
	Values = []
	for i in range(len(Seq)-Window+1):
		SubSeq = Seq[i:i+Window]
		Values.append(OligoDeltaG(SubSeq))
	return Values

def primerDimer(FP, RP):
	 return PrimerDimer.FindDimers(FP,RP)
	 
	
	
def Tm(Seq, Conc=0.2/1000000.0, SaltConc=50/1000.0):
	PrimerDeltaH = 0.0
	PrimerDeltaS = 0.0
	for i in range(len(Seq)-1):
		couple = Seq[i] + Seq[i+1]
		if couple in DeltaH.keys():
			PrimerDeltaH = PrimerDeltaH + DeltaH[couple]
			PrimerDeltaS = PrimerDeltaS + DeltaS[couple]
		else:
			raise KeyError ("unknown characters in the sequence.")
	PrimerDeltaH = PrimerDeltaH * (-1000)
	PrimerDeltaS = (PrimerDeltaS * (-1))-10.8
	
	Tmelting = (PrimerDeltaH/(PrimerDeltaS + (1.987*log(Conc/4.0)))) - 273.15 + 16.6*log10(SaltConc)
	return Tmelting

def TmProduct(Seq, SaltConc=50/1000.0):
	return 81.5 + 16.6*log10(SaltConc) + 0.41*GCcontent(Seq) - 675.0/len(Seq)

def mean(List):
	sum = 0.0
	for i in List:
		sum += i
	return sum/len(List)
	
class Primer:
	def __init__(self, Seq, Type):
		self.Conc = 0.2/1000000.0
		self.SaltConc = 50/1000.0
		if Type=="FP" or Type=="RP":
			self.Type = Type
		else:
			raise TypeError ('the value of Type can only be "FP" or "RP"')
		self.Seq = Seq.upper()
		self.Analyze()

	
	def Analyze (self):
		# length of the primer:
		self.Length = len(self.Seq)
		# GC content in percent
		self.GC = GCcontent(self.Seq)
		# calculate the melting temperature for the primer
		self.Tm = Tm(self.Seq, self.Conc, self.SaltConc)
		# calculate the internal stability
		Stability = DeltaGWindow(self.Seq)
		self.IntStab3 = mean(Stability[-5:])
		self.IntStab5 = mean(Stability[:5])
		


	def SetSaltConc (self, NewSaltConc):
		self.SaltConc = NewSaltConc/1000.0
	
	def SetPrimerConc (self, NewConc):
		self.Conc = NewConc/1000000.0
			
	#this function returns the 3' part of the primer that binds.
	# at least 10 bases must bind on the 3'.	
	def TemplateMatch(self,Sequence):
		if self.Type == "RP":
			Sequence = CompRev(Sequence)
		
		MustBind = self.Seq
		for i in range(len(self.Seq)-10):
			X = Sequence.find(MustBind)
			if X == -1:
				MustBind = self.Seq[i:]
			else:
				break
		if X == -1:
			self.Bind = "No Match"
			self.Hang = ""
		else:
			self.Bind = MustBind
			if self.Bind == self.Seq:
				self.Hang = ''
			else:
				self.Hang = self.Seq[:i-1]
		return self.Bind, self.Hang
		
	
class PCR_Product:
	def __init__(self, sequence):
		if sequence:
			self.sequence = sequence
			self.Tm = TmProduct(sequence)
			self.GC = GCcontent(sequence)
				
			
		else:
			self.sequence = "No Product"
			self.Tm = ""
			self.GC = ""
			self.TaOpt = ""
			self.AUCGC = ""
			self.RatioGC = ""
			self.NormAUCGC = ""
			self.AUCTM = ""
			self.RatioTM = ""
			self.NormAUCTM = ""
			
	def CalcTaOpt(self, FP_Tm, RP_Tm):
		self.TaOpt = 0.3 * min(FP_Tm, RP_Tm) + 0.7 * self.Tm - 14.9
		return self.TaOpt
		
	def AUCGC(self):
		self.AUCGC, self.RatioGC = AUC.GetAUC(self.sequence, GCcontent,21, 61)
		self.NormAUCGC = self.AUCGC * self.RatioGC
		return self.AUCGC, self.RatioGC, self.NormAUCGC
			
	def AUCTM(self):
		self.AUCTM, self.RatioTM = AUC.GetAUC(self.sequence, Tm, 21, 74)
		self.NormAUCTM = self.AUCTM * self.RatioTM
		return self.AUCTM, self.RatioTM, self.NormAUCTM
	
	
class PcrMachine:
	def __init__  (self, FP, RP, Template):
		self.FP = Primer(FP, "FP")
		self.RP = Primer(RP, "RP")
		self.Template = Template.upper()
		
	 # A function that receives the pars of the primers that must bind and returns the PCR sequence.
	def RawPCR(self, MustBindFP, MustBindRP, Sequence):
		if MustBindFP=="No Match" or MustBindRP=="No Match":
			return ""
		locationFP=Sequence.find(MustBindFP)
		locationRP=Sequence.find(CompRev(MustBindRP))
		if locationFP==-1 or locationRP==-1:
			return ''
		else:
			return Sequence[locationFP:locationRP+len(MustBindRP)]
	
	
	 # that's the action function which combine all the functions to one virtual pcr machine.
	def MakePCR(self, Frame=0):
		self.FP.TemplateMatch(self.Template)
		self.RP.TemplateMatch(self.Template)
		
		if self.FP.Bind != "No Match":
			self.RP.TemplateMatch(self.Template)
			if self.RP.Bind != "No Match":
				Result= self.RawPCR(self.FP.Bind, self.RP.Bind, self.Template)
				Product = self.FP.Hang + Result + CompRev(self.RP.Hang)
				return Product
	
		
	
"""
FPrimer="AAAAAGCAGGCTTGGTTAAGCTGAATGAACGAATA"
RPrimer="AGAAAGCTGGGTATGACTCTGGCAACTGATT"
MyTemplate="ATGGCATTCCGGACAATTTGCGTGTTGGTTGGAGTATTTATTTGTTCTATCTGTGTGAAAGGATCTTCCCAGCCCCAAGCAAGAGTTTATTTAACATTTGATGAACTTCGAGAAACCAAGACCTCTGAATACTTCAGCCTTTCCCACCATCCTTTAGACTACAGGATTTTATTAATGGATGAAGATCAGGACCGGATATATGTGGGAAGCAAAGATCACATTCTTTCCCTGAATATTAACAATATAAGTCAAGAAGCTTTGAGTGTTTTCTGGCCAGCATCTACAATCAAAGTTGAAGAATGCAAAATGGCTGGCAAAGATCCCACACACGGCTGTGGGAACTTTGTCCGTGTAATTCAGACTTTCAATCGCACACATTTGTATGTCTGTGGGAGTGGCGCTTTCAGTCCTGTCTGTACTTACTTGAACAGAGGGAGGAGATCAGAGGACCAAGTTTTCATGATTGACTCCAAGTGTGAATCTGGAAAAGGACGCTGCTCTTTCAACCCCAACGTGAACACGGTGTCTGTTATGATCAATGAGGAGCTTTTCTCTGGAATGTATATAGATTTCATGGGGACAGATGCTGCTATTTTTCGAAGTTTAACCAAGAGGAATGCGGTCAGAACTGATCAACATAATTCCAAATGGCTAAGTGAACCTATGTTTGTAGATGCACATGTCATCCCAGATGGTACTGATCCAAATGATGCTAAGGTGTACTTCTTCTTCAAAGAAAAACTGACTGACAATAACAGGAGCACGAAACAGATTCATTCCATGATTGCTCGAATATGTCCTAATGACACTGGTGGACTGCGTAGCCTTGTCAACAAGTGGACCACTTTCTTAAAGGCGAGGCTGGTGTGCTCGGTAACAGATGAAGACGGCCCAGAAACACACTTTGATGAATTAGAGGATGTGTTTCTGCTGGAAACTGATAACCCGAGGACAACACTAGTGTATGGCATTTTTACAACATCAAGCTCAGTTTTCAAAGGATCAGCCGTGTGTGTGTATCATTTATCTGATATACAGACTGTGTTTAATGGGCCTTTTGCCCACAAAGAAGGGCCCAATCATCAGCTGATTTCCTATCAGGGCAGAATTCCATATCCTCGCCCTGGAACTTGTCCAGGAGGAGCATTTACACCCAATATGCGAACCACCAAGGAGTTCCCAGATGATGTTGTCACTTTTATTCGGAACCATCCTCTCATGTACAATTCCATCTACCCAATCCACAAAAGGCCTTTGATTGTTCGTATTGGCACTGACTACAAGTACACAAAGATAGCTGTGGATCGAGTGAACGCTGCTGATGGGAGATACCATGTCCTGTTTCTCGGAACAGATCGGGGTACTGTGCAAAAAGTGGTTGTTCTTCCTACTAACAACTCTGTCAGTGGCGAGCTCATTCTGGAGGAGCTGGAAGTCTTTAAGAATCATGCTCCTATAACAACAATGAAAATTTCATCTAAAAAGCAACAGTTGTATGTGAGTTCCAATGAAGGGGTTTCCCAAGTATCTCTGCACCGCTGCCACATCTATGGTACAGCCTGTGCTGACTGCTGCCTGGCGCGGGACCCTTATTGCGCCTGGGATGGCCATTCCTGTTCCAGATTCTACCCAACTGGGAAACGGAGGAGCCGAAGACAAGATGTGAGACATGGAAACCCACTGACTCAATGCAGAGGATTTAATCTAAAAGCATACAGAAATGCAGCTGAAATTGTGCAGTATGGAGTAAAAAATAACACCACTTTTCTGGAGTGTGCCCCCAAGTCTCCGCAGGCATCTATCAAGTGGCTGTTACAGAAAGACAAAGACAGGAGGAAAGAGGTTAAGCTGAATGAACGAATAATAGCCACTTCACAGGGACTCCTGATCCGCTCTGTTCAGGGTTCTGACCAAGGACTTTATCACTGCATTGCTACAGAAAATAGTTTCAAGCAGACCATAGCCAAGATCAACTTCAAAGTTTTAGATTCAGAAATGGTGGCTGTTGTGACGGACAAATGGTCCCCGTGGACCTGGGCCAGCTCTGTGAGGGCTTTACCCTTCCACCCGAAGGACATCATGGGGGCATTCAGCCACTCAGAAATGCAGATGATTAACCAATACTGCAAAGACACTCGGCAGCAACATCAGCAGGGAGATGAATCACAGAAAATGAGAGGGGACTATGGCAAGTTAAAGGCCCTCATCAATAGTCGGAAAAGTAGAAACAGGAGGAATCAGTTGCCAGAGTCATAA"

a_FPrimer = Primer(FPrimer, "FP")
print a_FPrimer.GC
print a_FPrimer.Length

Reaction = PcrMachine(FPrimer, RPrimer, MyTemplate)
print Reaction.FP.Tm
print Reaction.RP.Tm
X = Reaction.MakePCR()
print Reaction.CheckFrame(X)

#print Reaction.FP.PerfectMatch(Reaction.Template)
#print Reaction.RP.PerfectMatch(Reaction.Template)
"""