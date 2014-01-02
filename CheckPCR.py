"""
Compute parameters to predict the success of PCR (GC content, dimers, reaction temperature etc.). 

The porgram does not design primers, but provides a fast and efficient way of comparing detailed statistics of the PCR reactions for proposed primers

To design primers, use the PickPCR module, which uses a custom implementation of the EMBOSS eprimer3 program. 

INPUT: A file with organism name, strain, gene, sense and antisense primer, target region
RETURN: Text tab delimited file containing: Bind, Outflanks, Tempreature, GC ratio, AA sequence, product sequence etc. 

Source: Existing code from Y. Benita was expanded to download reference sequence using Entrez, loop over a simpler input file... 
in addition, the output formatting was tweaked in several small ways (return the AA sequence of interest, show amplicon size etc.) 
Original Benita code is available http://xavierlab.mgh.harvard.edu/benita/PCRSupInfo/ (Accessed in December 2013)

For testing: 
Put all the file
os.chdir('/Users/nbraykov/Google Drive/SYNC w HOME/Grad Skool/IBS 796/St00pidPCR')
input_filename = 'cip_test.txt'
"""

from Bio import Entrez
import re
import urllib
import PCR
from PCR import VirtualPCR
from PCR import PrimerDimer
from PCR import AUC
import os





def MakeString(List):
	MyString = ""
	for i in List:
		if type(i) is type(0.0):
			MyString += "%.2f" % i + "\t"
		else:
			MyString += str(i) + "\t"
	return MyString[:-1]

def EntrezCompile (SearchTerm):
	Entrez.email = 'nbrayko@emory.edu'
	# accession id works, see here for more 
	handle=Entrez.esearch(db='gene', retmax=10, term=SearchTerm )
	result = Entrez.read(handle, validate=True)
	if int(result["Count"]) !=0:
		GeneRefNum = result["IdList"][0]
		handle = Entrez.esummary(db="gene", id=GeneRefNum)
		temp=Entrez.read(handle, validate=True)
		temp=temp[0]
		GeneRefSpecies = temp['Orgname'].replace(" ", "_").replace(".", "_").replace(",", "_")
		GeneRefName = temp['Name'].replace(" ", "_").replace(".", "_").replace(",", "_")
		print "Found %s result(s) sorted by relevance; \n returning the first one, corresponding to %s, %s (%s)" % (result["Count"], GeneRefSpecies, GeneRefName, temp['Description'])
		GenomicPos=temp["GenomicInfo"][0]
		return(GeneRefNum, GeneRefName, GeneRefSpecies, GenomicPos)
	else:
		print "Note: no results for %s \n Repeat search" % (SearchTerm)	


def GeneLookup(Species, Strain, Gene, Genespec =0):
	"""
	Searches NCBI Gene database for the reference number of the gene sequence of interest, returns the first result; also, prints a list of other available result
	Search term has to be submitted as ("staphylococcus aureus"[Organism]) AND "gyra"[Gene Name] AND "newman"
	
	Genespec set to 0, meaning that if sequence is not found for the strain specified in the input file, a search will be perform w/o that parameter. 
	"""
	SearchTerm = Species+"[Organism] AND "+Strain+" AND "+Gene+"[Gene]"
	#if more than 1 result returned, fetch the ID for teh first one, look it up and print all the dets for verification 
	handle=Entrez.esearch(db='gene', retmax=10, term=SearchTerm )
	result = Entrez.read(handle, validate=True)
	if int(result["Count"]) ==0:
		if Genespec ==1:
			print "Did not find anything for %s \n YOU TRIED, lol" % (SearchTerm)			
		else:
		#attemps search without strain-specific info if Genespec =0
			SearchTerm2 =Species+"[Organism] AND "+Gene+"[Gene]"
			print "Repeating search for %s \n WARNING: not strain-specific" % (SearchTerm2)	
			return(EntrezCompile (SearchTerm2))			
	else:
		if Strain == '':
			print "WARNING: Entrez search not strain-specific"
		return(EntrezCompile(SearchTerm))	

			
		
def GeneTemplate(GeneRefNum, GeneRefName, GeneRefSpecies, GenomicPos, Outflank=0):
	"""
	Queries NCBI using the Gene ID, downloads the reference gene in FASTA format The gene sequence is in caps, the 'outflank' non-coding regions are in small letters
	"""
	# ref.: http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html
	# replace with your real email (optional):
	Entrez.email = 'nbrayko@emory.edu'
	
	outflank = int(Outflank)
	# accession id works, see here for more 
	handle=Entrez.efetch(db='nucleotide', id=GenomicPos['ChrAccVer'], seq_start = GenomicPos['ChrStart']+1-outflank, seq_stop=GenomicPos['ChrStop']+1+outflank, rettype="fasta", retmode="text" )
	temp_seq = handle.read()
	# store locally:
	template_filename = 'Template_Seq'+str(outflank)+'_'+GeneRefName+'_'+GeneRefSpecies +'.txt'
	template_seqf=open(template_filename, 'w')
	template_seqf.write(temp_seq)
	handle.close()
	template_seqf.close()
	return(temp_seq, template_filename)
	
def IntermLineCompile(template_filename, primer_sense, primer_antisense, interest_regions):
	"""
	Merges input of primer sequences with the gene template; returns  FP, RP, Template 	that is to be fed into the VirtualPCR; 
	"""
	template = open(template_filename)
	Outflank = re.findall('\d+', template_filename.split('_')[1])[0]
	outflank = int(Outflank)
	template.readline()
	#read entire FASTA sequence as one string (includes hang
	gene_seq_h=[seqs.replace('\n','') for seqs in template.readlines()]
	gene_seq_h=''.join(gene_seq_h)
	template.close()
	#gene_seq excludes the hang:
	if outflank!=0:
		gene_seq= gene_seq_h[outflank:-outflank]
		#gene_seg_f is formatted, combining the hangs as small caps
		gene_seq_f = gene_seq_h[:outflank].lower()+ gene_seq + gene_seq_h[-outflank:].lower()
	else:
		gene_seq = gene_seq_h
		gene_seq_f = gene_seq
	#translate codon sequence to AA sequence w numbers
	aa_seq = VirtualPCR.SeqToAA(gene_seq)
	#visualize the target relativet o entire aa sequence; basically, scan the aa sequence, draw it as a ---, replace teh dash w the aa target name if it matches.  
	if interest_regions != ['']:
		aa_viz = ['*']*aa_seq.count('_')
		for region in interest_regions:
			for aa in aa_seq.split('_'):
				#if there is a perfect match, return them in uppercase
				if aa ==region.strip():
					aa_viz[aa_seq.split('_').index(aa)] = aa.upper() 
				#if just numbers match, return in lowercase
				elif re.findall('\d+',aa)==re.findall('\d+',region.strip()):
					aa_viz[aa_seq.split('_').index(aa)] = aa
				else:					
					next 
		aa_viz = ''.join(aa_viz)			
	#if no target region is specified, translate entire aa sequence
	else:
		aa_viz = aa_seq	
	#if outflanked, add '-' to the ends
	aa_viz = outflank/3*'-'+aa_viz+outflank/3*'-'
	#find the positions of the sense and anti-sense primers relative to template; if not found, return (No exact match)
	sense_postest = re.search(primer_sense, gene_seq, flags=0)
	if sense_postest == None:
		sense_pos = "Not in range"
	else:
		sense_pos = sense_postest.span()
	#Note: use the reverse comp for the antisense check
	antisense_postest = re.search(VirtualPCR.CompRev(primer_antisense), gene_seq, flags=0)
	if antisense_postest == None:
		antisense_pos = "Not in range"
	else:
		antisense_pos = antisense_postest.span()
	#open the input file, read in the primer sequences
	InputLineData = [gene_seq_f, aa_viz, primer_sense, primer_antisense, sense_pos, antisense_pos]
	InputLineDataStr = MakeString(InputLineData)
	return(InputLineDataStr)

def IntermInputCompile(input_filename, Outflank= 0):
	"""
	Parses the input, use all the functions above. Outflank option is overrided by the Outflank speciefied in input file, if one is specified at all. 
	"""
	outflank = int(Outflank)
	infile = open(input_filename +'.txt', 'rU')
	infile.readline()
	outfile=open('_Interm_Output'+input_filename + '.txt', 'w')
	counter = 1
	gen_count = ''
	for line in infile:
		Species, Strain, Gene, Outflank, primer_sense, primer_antisense, interest_regions  = line.split('\t')
		if Outflank!='': 
			outflank = int(Outflank)
		primer_antisense = primer_antisense.rstrip('\n').rstrip('\t')
		interest_regions = interest_regions.strip('\n' '\t' '"').split(',')
		if gen_count != Gene:
			temp_query = GeneLookup(Species, Strain, Gene)
		temp_data = GeneTemplate(temp_query[0], temp_query[1], temp_query[2], temp_query[3], outflank)
		outline = Gene+'_'+str(counter) +'\t'+IntermLineCompile(temp_data[1], primer_sense, primer_antisense, interest_regions) + '\t'+ temp_data[1]+'\n'
		counter = counter+1
		gen_count = Gene
		outfile.write(outline)
		
	outfile.close()
	infile.close()
	

def runCheckPCR(input_filename, Outflank=0):
	'''
	Runs the program to check in-silico PCR stats for an input file modeled after cip_test; 
	Outflank parameter specifies an area upstream and downstream from gene's coding region that will also be downloaded. 
	The purpose of this feature is if you want primers that bind to non-coding regions to flank the entire gene 
	NOTE: outflank is different from the target_region, which is the mutation of interest	
	Program does not design primers, it just compiles their detailed chracteristics to assist in picking. 
	For primer design, run the PickPCR module, which will return a file that looks like cip_test_pr3
	'''
	#input_filename = cip_test_lit
	IntermInputCompile(input_filename, Outflank)
	interm_input_filename = '_Interm_Output'+input_filename + '.txt'
	infile = open (interm_input_filename,"rU")
	Header = ["Label",
	"Forward Primer", "FP Bind", "FP Hang", "FP Tmelting", "FP Int Stab 3", "FP Int Stab 5", "FP GC content", "FP GC Diff",
	"Reverse Primer", "RP Bind", "RP Hang", "RP Tmelting", "RP Int Stab 3", "RP Int Stab 5", "RP GC content", "RP GC Diff",
	"Template", "Template filename" ,"PCR Product", "Product length", "Product lcoation", "Product Tmelting", "Product GC", "Product TaOpt", "Product AUCGC", "Product Ratio GC", "Product NormAUCGC", "Product AUCTM", "Product Ratio TM", "Product NormAUCTM",
	"Primer Dimer Max Stability", "Primer Dimers", "Primer Tmelting diff", "Template AA sequence", "Product AA sequence"]
	outfile=open('_PCR_result_' +input_filename+'.txt', 'w')
	outfile.write(MakeString(Header)+'\n')
	for line in infile:
		Label, Template, AASequence, FP, RP, PosSense, PosAnitsense, TemplFile = line.split('\t')
		X = VirtualPCR.PcrMachine(FP, RP, Template.upper())
		# make an instance of PCR_Product
		P = VirtualPCR.PCR_Product(X.MakePCR())
		if P.sequence == 'No Product':
			Data = [Label, FP, RP, '-> **** FAIL. No Product. Nice try.****']
		else:
			P.CalcTaOpt(X.FP.Tm, X.RP.Tm)
			P.AUCGC()
			P.AUCTM()
			Dimer =	VirtualPCR.primerDimer(X.FP.Seq, X.RP.Seq)
			Dimer.sort(lambda x,y: cmp(-x.Stability, -y.Stability))	
			if len(Dimer)>=1:
				PrimerDimerStab = Dimer[0].Stability
				PrimerDimerSeq = Dimer[0].MatchString.replace('\n','')
			else:
				PrimerDimer = 0.0
				PrimerDimerSeq='No dimers'	
			#the sequence and positon of amplified region
			product_lookup = re.search(P.sequence, Template.upper(), flags=0)
			if product_lookup == None:
				P_location = "Not found"
				AAproduct = "Not found"
			else:
				P_location = str(product_lookup.span())
				AAproduct = AASequence[int(product_lookup.span()[0]/3)-5:int(product_lookup.span()[1]/3)+5]		
			Data = [Label,FP, X.FP.Bind, X.FP.Hang, X.FP.Tm, X.FP.IntStab3, X.FP.IntStab5, X.FP.GC,X.FP.GC/P.GC, 
				RP, X.RP.Bind, X.RP.Hang, X.RP.Tm, X.RP.IntStab3, X.RP.IntStab5, X.RP.GC,X.RP.GC/P.GC,Template, 
				TemplFile.strip('\n'), P.sequence, len(P.sequence),  P_location,  P.Tm, P.GC, P.TaOpt, P.AUCGC, P.RatioGC, 
				P.NormAUCGC, P.AUCTM, P.RatioTM, P.NormAUCTM, PrimerDimerStab, PrimerDimerSeq, abs(X.FP.Tm-X.RP.Tm), 
				AASequence, AAproduct]	
		outfile.write(MakeString(Data)+'\n')		
	outfile.close()	
	infile.close()
	


