from PCR import VirtualPCR
from sys import argv

def MakeString(List):
	MyString = ""
	for i in List:
		if type(i) is type(0.0):
			MyString += "%.2f" % i + "\t"
		else:
			MyString += str(i) + "\t"
	return MyString[:-1]

# in this line write the path to your input file
InputFile =	argv[1]

# output style can be either Tab or Human. Tab is a tab delimited table and
# Human is for a small number of sequences.
OutputStyle = argv[2]


file = open (InputFile,"r")

# file header
#print "Label\tForward Primer\tFP Bind\tFP Hang\tFP Tmelting\tFP Int Stab 3\tFP Int Stab 5\tFP GC content\tFP GC Diff\tReverse Primer\tRP Bind\tRP Hang\tRP Tmelting\tRP Int Stab 3\tRP Int Stab 5\tRP GC content\tTemplate\tPCR Product\tProduct Tmelting\tProduct GC\tProduct TaOpt\tProduct AUCGC\tProduct Ratio GC\tProduct NormAUCGC\tProduct AUCTM\tProduct Ratio TM\tProduct NormAUCTM\tPrimer Dimer\tPrimer Tmelting diff\tPrimer GC Diff"
Header = ["Label",
"Forward Primer", "FP Bind", "FP Hang", "FP Tmelting", "FP Int Stab 3", "FP Int Stab 5", "FP GC content", "FP GC Diff",
"Reverse Primer", "RP Bind", "RP Hang", "RP Tmelting", "RP Int Stab 3", "RP Int Stab 5", "RP GC content", "RP GC Diff",
"Template", "PCR Product", "Product Tmelting", "Product GC", "Product TaOpt", "Product AUCGC", "Product Ratio GC", "Product NormAUCGC", "Product AUCTM", "Product Ratio TM", "Product NormAUCTM",
"Primer Dimer", "Primer Tmelting diff"]

print MakeString(Header)

for mybuffer in file.xreadlines():
	Label, FP, RP, Template = mybuffer.split()

	X = VirtualPCR.PcrMachine(FP, RP, Template)
	# make an instance of PCR_Product
	P = VirtualPCR.PCR_Product(X.MakePCR())
	P.CalcTaOpt(X.FP.Tm, X.RP.Tm)
	P.AUCGC()
	P.AUCTM()
	Dimer =	VirtualPCR.primerDimer(X.FP.Seq, X.RP.Seq)
	Dimer.sort(lambda x,y: cmp(-x.Stability, -y.Stability))
	
	if len(Dimer)>=1:
		PrimerDimer = Dimer[0].Stability
	else:
		PrimerDimer = 0.0

	Data = [Label,
			FP, X.FP.Bind, X.FP.Hang, X.FP.Tm, X.FP.IntStab3, X.FP.IntStab5, X.FP.GC,X.FP.GC/P.GC,
			RP, X.RP.Bind, X.RP.Hang, X.RP.Tm, X.RP.IntStab3, X.RP.IntStab5, X.RP.GC,X.RP.GC/P.GC,
			Template, P.sequence, P.Tm, P.GC, P.TaOpt, P.AUCGC, P.RatioGC, P.NormAUCGC, P.AUCTM, P.RatioTM, P.NormAUCTM,
			PrimerDimer,abs(X.FP.Tm-X.RP.Tm)]
	
	if OutputStyle == 'Tab':
		print MakeString(Data)
	
	elif OutputStyle == 'Human':
		print "-"*20 + " %s: %s " % (Header[0], Data[0]) + "-"*20
		for i in range(1,len(Header)):
			if type(Data[i]) is type(0.0):
				print Header[i] + ":\t" + "%.2f" % Data[i]
			else:
				print Header[i] + ":\t" + str(Data[i])
		print "-"* (40+len(Header[0])+len(Data[0])+4) + "\n"
	else:
		raise TypeError("The output style must be set to Tab or Human.")	
	
	
file.close()