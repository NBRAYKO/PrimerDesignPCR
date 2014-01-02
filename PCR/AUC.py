from math import log, log10

class HighRegion:
	def __init__ (self, Start, ValueList ,WindowSize,CutOff):
		self.Start = Start+1
		self.ActualStart = self.Start + WindowSize/2 # use this start so that overlaps don't occur.
		self.Length = len(ValueList)
		self.End = self.Start + self.Length
		self.ActualEnd = self.End - WindowSize/2
		self.Values = ValueList
		self.Auc = AucAboveLine(ValueList, CutOff)
	

# This class hold information on all the regions in one DNA sequence.
class AllRegions:
	def __init__(self,Label,RegionList):
		self.Label = Label # this is just the name of the sequence.
		self.NumOfRegions = len(RegionList)
		self.Regions = {}
		
		for i in RegionList:
			self.Regions[i.ActualStart] = i 	
		
		self.StartValues = self.Regions.keys()
		self.StartValues.sort()
		self.TotalAuc, self.Length = self.SumRegions()
	
	# Here I summerize information from all regions. Its better if you want to characterize per sequence.		
	def SumRegions(self):
		AUC = 0.0
		Length = 0
		for i in self.Regions.keys():
			AUC += self.Regions[i].Auc
			Length += self.Regions[i].Length
		return AUC, Length
		

	def MaxRegion(self):
		MaxAUC = 0
		MaxRegion = 0
		for i in self.Regions.keys():
			if self.Regions[i].Auc > MaxAUC:
				MaxAUC = self.Regions[i].Auc
				MaxRegion = self.Regions[i].Auc
		self.MaxRegion = MaxRegion
	


# This function gets a list and makes a region list. Here I also specify criterions for the region.
# It returns a dictionary: key = start position, value = list of Tm values.
def AnalyzeList(ValueList, Edge):
	RegionsList = {}
	Region = []
	for i in range(len(ValueList)):
		if ValueList[i] >= Edge:
			Region.append(ValueList[i])
		else:
			if len(Region) >= 10:
				RegionsList[i-len(Region)] = Region
			Region = []
	if len(Region) >= 10:
		RegionsList[i-len(Region)] = Region
	return RegionsList


def AucAboveLine(List,Edge = 0):
	AUC = 0.0
	for i in range(len(List)-1):
		AUC += (List[i] + List[i+1] - Edge * 2.0)/2.0
	AUC += (List[0] - Edge)/2.0 + (List[-1] - Edge)/2.0
	return AUC
	

# The function that generates the list of values using a sliding window.
def SlidingWindow (Seq, Func, Window):
	List = []
	for i in range(len(Seq)-Window+1):
		SubSeq=Seq[i:i+Window]
		X = Func(SubSeq)
		List.append(X)
	return List

def GetAUC(Seq, Func,Window, CutOff):
	X = SlidingWindow(Seq, Func, Window)
	Xlist = AnalyzeList(X, CutOff)
	Xobjects = []
	for i in Xlist.keys():
		Y = HighRegion(i, Xlist[i], Window, CutOff)
		Xobjects.append(Y)
	
	All = AllRegions("Test", Xobjects)
	return All.TotalAuc, float(All.Length)/len(X) 

"""
# ---------- example--------------
#Seq = "CCTGGACTCAACAGCCACGACCCGCACTCGGACGAGGACACGCCGACGTCGGACGACCTGGAGCAGTTCGCCAAGCAGTTCAAGCAGCGGCGCATCAAGCTGGGCTTCACGCAGGCCGACGTGGGGTTGGCGCTGGGCACACTCTACGGCAACGTGTTCTCGCAGACCACCATCTGCCGCTTCGAGGCCCTGCAGCTGAGCTTCAAGAACATGTGCAAGCTCAAGCCGCTGCTGAACAAGTGGCTGGAGGAGGCGGACTCAAGCACCGGCAGCCCCACAAGCATCGACAAGATCGCGGCGCAGGGCCGCAAGCGCAAGAAGCGGACCTCTATCGAGGTGAGCGTCAAGGGCGCGCTGGAGAGCCACTTCCTCAAGTGCCCCAAGCCCTCCGCGCAGGAGATCACCAACCTGGCCGACAGCCTGCAGCTCGAGAAGGAGGTGGTGCGGGTCTGGTTCTGCAATCGGCGCCAAAAGGAGAAGCGCATGACGCCGCCCGGGATCCAACAGCAGACGCCCGACGACGTCTACTCGCAGGTGGGCACCGTGAGCGCCGACACGCCGCCGCCTCACCACGGGCTGCAGACGAGCGTTCAG"
file = open ("/Users/ybenita/Desktop/AllSeq.tab","r")

for mybuffer in file.xreadlines():
	Label, dummy, Seq = mybuffer.split()
	AUCGC, RatioGC = GetAUC(Seq, GCcontent, 61)
	AUCTM, RatioTM = GetAUC(Seq, Tm, 75)
	
	print "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (Label, AUCGC, RatioGC, (AUCGC*RatioGC), AUCTM, RatioTM, (AUCTM*RatioTM))
	#break
	
file.close()
"""