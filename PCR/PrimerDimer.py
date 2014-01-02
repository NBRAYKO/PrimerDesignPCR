"""
This script will find primer dimers. Basically it slides 2 sequences one over the other and find the
best way they can bind. The information is stored in the Dimer object.
By calling the function FindDimers (ForwardPrimer, ReversePrimer) you get back a list ofk
Dimer objects.
FP - The full sequence of the forward primer
RP - The full sequence of the reverse primer
FPfrag - the fragment from the forward primer which produced the hit.
RPfrag - the fragment from the reverse primer which produced the hit.
domains - a list of tuples of consequtive matches [(start,end)] 
Total - Total number of matches including all domains.
TotalPercent - number of matches in percent.
BindingSeq - The sequence that matches
MatchString - a 3 line string showing how the FP and RP match. Blast style.
"""


import VirtualPCR

RevDNA={'A':'T','T':'A','G':'C','C':'G','N':'N'}
#compl reverse seq
def CompReverse(Sequence):
	NewSeq=''
	for i in Sequence:
		NewSeq=RevDNA[i]+NewSeq
	return NewSeq	
	
# Reverse the oreder of a sequenc	
def Rev(Seq):
	NewSeq = ""
	for i in Seq:
		NewSeq = i + NewSeq
	return NewSeq	



class Dimer:
	def __init__ (self,HitTuple, FP, RP): #, Inverse=0):
		self.FP = FP
		self.RP = RP
		#self.Inverse = Inverse
		self.FPfrag = HitTuple[1]
		self.RPfrag = HitTuple[2]
		self.domains = HitTuple[0][0]
		self.Total = HitTuple[0][1]
		self.TotalPercent = float(self.Total)/len(self.FPfrag)
		self.MatchString = self.MakeString()
		self.DimerStability()
	
		
	def DimerStability (self):
		OnlyBinding = ""
		for i in self.domains:
			OnlyBinding += self.FPfrag[i[0]:i[1]+1]
		self.BindingSeq = OnlyBinding
		self.Stability = VirtualPCR.OligoDeltaG(self.BindingSeq)
	
	# make the match string
	def MakeString(self):
		StartFP = self.FP[:self.FP.rfind(self.FPfrag)]
		self.Start = len(StartFP)
		self.Length = self.Start + len(self.RP)
		TheString = self.FP+"\n"+ " "*len(StartFP) + "|"
		MatchList = []
		for i in self.domains:
			MatchList += range(i[0],i[1]+1)
		for q in range(1,len(self.FPfrag)-1):
			if q in MatchList:
				TheString+="|"
			else:
				TheString+=" "
		TheString += ("|\n" + " "*len(StartFP) + Rev(self.RP))
		return TheString




# This is the function that you need to call. Seq1 is the Forward primer and Seq2 is the Reverse.
# Edge is the number of consequtive matches you want to have in order to define a domain.	


def FindDimers (Seq1,Seq2,Edge=0):
	Match = []
	Seq2c = CompReverse(Seq2)
	counter = 1
	#for i in range(len(Seq1)):
	while 1:
		try:
			Primer1 = Seq1[-3-counter:]
			Primer2 = Seq2c [:3+counter]
		except IndexError:
			break
		if len(Primer1)==len(Seq1):
			break
		counter += 1
		#print Primer1
		#print Primer2, "\n"
		Test = CompareSeqs(Primer1, Primer2, Edge)
		if len(Test) > 0:
			X = Dimer((Test, Primer1,Primer2),Seq1, Seq2)
			Match.append(X)
	return Match


# This function only computes the number of matches
def CompareSeqs (Seq1,Seq2,Edge):
	if Seq1[0] != Seq2[0] or Seq1[-1] != Seq2[-1]:
		return []
	Identical = 0
	Matches = []
	for i in range(len(Seq1)):
		try:
			if Seq1[i]==Seq2[i]:
				Matches.append(i)
		except IndexError:
			break
	if len(Matches) > 0:
		return Consequtive(Matches, Edge)
	else:
		return []
	

def Consequtive (List, Edge):
# This function gets a list of matches and returns a list of tupples [(start,end)]
	start = List[0]
	next = List[0]+1
	domains = []
	for i in List[1:]:
		if i == next:
			next+=1
		else:
			if next-start-1 >=Edge:
				domains.append((start,next-1))
			start = i
			next = start +1
		if i==List[-1] and (i-start >= 3):
			domains.append((start,i))
	Total = 0
	for i in domains:
		Total += (i[1]-i[0]+1)
	if Total >= 3:
		return domains, Total
	else:
		return []


"""
FP = "AAAAAGCAGGCTTGAGCAATAATGACGATGAGGA"
RP = "AGAAAGCTGGGTATCGGACTTGCTGCTCCAGACT"

#print FP
#print RevRP
a = FindDimers(FP, RP,0)
for i in a:
	print i.FP
	print i.RP
	print i.FPfrag
	print i.RPfrag
	print i.domains
	print i.Total
	print i.TotalPercent
	print i.MatchString
	print i.BindingSeq
	print i.Stability
"""