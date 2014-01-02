'''
Set of functions to look up and download fasta files and design primers for particular regions using Primer3. 
The input file for the main function PrimerParser has to be modeled after cip_test_pr3 and contain organism, strain (optional), gene and target mutation site (optional)

Target mutation site has to specify at least a codon number. The exact name of the amino-acid can be in any form (e.g., Glu-84, G84 etc.). 
Multiple codon numbers can be entered if separated by comma. If >1 codon is entered as a target site, program picks the midpoint 
(so, if interested in sites too far apart, make it a separate query) 

Notes on how to set up the dependent scripts for Emboss and Primer3 (this is srsly byzantine):
1. INSTALL PRIMER3 using brew https://github.com/Homebrew/homebrew-science
2. INSTALL EMBOSS (need the eprimer3 interface script in addition to primer3?), use brew install https://raw.github.com/Homebrew/homebrew-science/master/emboss.rb
3. use the biopython Primer3Commandline interface to run things from python
4. Note: to get things running, need to have the primer3_config dir copied in the current working dir
5. Note: also, need to download the old version of primer3 1.xx, and use terminal to set the emboss path to recognize this as teh primer3_core file'\
Type $ export EMBOSS_PRIMER3_CORE="/usr/local/Cellar/primer3/1.1.4/bin/primer3_core"
6. finally, this python script raises exceptions   File "//anaconda/lib/python2.7/site-packages/Bio/Application/__init__.py", line 444, in __call__
Despite the fact the program runs fine in bash; however, this environmental variable is not stored when running python, so have to enter the command os.environ['EMBOSS_PRIMER3_CORE'] before each run

Finally, had to modify line 44 in '/Bio/Application/__init__.py' so it does not raise errors based on Eprimer outpit codes, as it was not interpreting them correctly 
'''
import os
import Bio.Emboss.Applications
Bio.Emboss.Applications.Primer3Commandline()
from Bio.Emboss.Primer3 import *
from Bio import SeqIO
import CheckPCR
import re
from Bio import Entrez

#point to where the primer_3 file is; have to use the older 1.1 primer3 version as emboss can;t handle newer ones. 
os.environ['EMBOSS_PRIMER3_CORE'] = '/usr/local/Cellar/primer3/1.1.4/bin/primer3_core'


def Primer3Picker(seq_file, target_codon=0, out_num = 10,prod_length=200):
	'''
	Takes in fasta sequence file; uses emboss and primer3 commandline tool to 
	design primers. User can specify product length, and target sites for mutation regions
	by entering the codon number. Only single codon number can be entered. 
	
	Default product length is 200; target site is 0, meaning Primer3 picks the "most optimal" 
	Returns a .pr3 file and a status message
	'''
	cline = Bio.Emboss.Applications.Primer3Commandline(sequence=seq_file, task =1)
	in_root = str(seq_file.split('.')[0])
	cline.outfile = in_root+'_site'+str(target_codon)+'.pr3'
	cline.numreturn = out_num
	cline.maxsize = 26
	#cline.otm = 58
	cline.mintm = 52
	cline.mingc = 35
	cline.maxgc = 75
	cline.psizeopt = prod_length
	cline.prange = "100-400"
	region_size = int(cline.psizeopt/2)
	#default behavior is to ask for ~200bp amplicon such  that the target mutation is right in the middle
	if target_codon!=0:
		cline.target = str(3*(target_codon-region_size/3)) + ',' + str(3*(target_codon+region_size/3))     
	stdout, stderr = cline()
	outfile = cline.outfile
	return(outfile, stderr)


#write function to loop over pr3 files and collate the output. adapted from https://github.com/Jingping/BiteTools/blob/d0d34e792f4c952373653220b0d3ac24bbaac33a/primer_pick.py
def PrimerFileComp(pr3_path):
	'''
	Gets the primer3 file from path  (e.g.,'Template_Seq_gyrA0_Staphylococcus_aureus_subsp__aureus_str__Newman_site0.pr3'), 
	parse it and spit out text files of proposed primers for each gene and site
	'''
	gen_id = pr3_path.split('_')[2]+pr3_path.split('_')[-1].strip('.pr3')+'_pr3'
	outfile_path = gen_id+"_output.txt"
	outfile = open(outfile_path, "w")
	outfile.write("rec_id\tforward_seq\treverse_seq\tforward_start\treverse_start\tproduct_length\n")
	pr3_outfile = file(pr3_path, "r")
	primer_record = read(pr3_outfile)
	counter = 1
	for primer in primer_record.primers:	
		rec_id = gen_id + '_'+str(counter)
		product_len = -primer.forward_start+primer.reverse_start+primer.reverse_length
		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % \
                     (rec_id,primer.forward_seq,primer.reverse_seq,primer.forward_start,primer.reverse_start,product_len))
		counter =  counter+1
	outfile.close()
	return(outfile_path)

def PrimerParser(input_filename, out_num=10, prod_length=200):
	'''
	Takes input that looks like cip_test_pr3; 
	
	The input file for the main function PrimerParser has to be modeled after cip_test_pr3 and contain organism, strain (optional), gene and target mutation site (optional)
	Target mutation site has to specify at least a codon number. The exact name of the amino-acid can be in any form (e.g., Glu-84, G84 etc.). 
	Multiple codon numbers can be entered if separated by comma. If >1 codon is entered as a target site, program picks the midpoint 
	(so, if interested in sites too far apart, make it a separate query)if no value is given for target mutation, default primer3 'optimal product' is returned
	'''
	infile = open(input_filename +'.txt', 'rU')
	infile.readline()
	outfile=open('_pr3_out_'+input_filename + '.txt', 'w')
	#create separate text files with each primer3 output for each line
	outfile.write(("Organism\tStrain\tGene\tOutflank bp size\tSense primer\tAntisense primer\tFlanked region (codons; mid-point if >1)\n"))
	#note: all the downloaded genomes have outflank=0
	outflank=0
	for line in infile:
		Species, Strain, Gene, Outflank, interest_regions, Entire  = line.split('\t')
		interest_regions_l = re.findall('\d+', interest_regions)
		interest_regions_l = [int(l) for l in interest_regions_l]
		if interest_regions_l !=[]:
			target_codon = int(sum(interest_regions_l)/len(interest_regions_l))
		else:
			target_codon=0
		temp_query = CheckPCR.GeneLookup(Species, Strain, Gene)
		temp_data = CheckPCR.GeneTemplate(temp_query[0], temp_query[1], temp_query[2], temp_query[3], Outflank=outflank)
		#if Entire!=' \n':
		#	print "Looking up primers to sequence the entirety of %s " % (Gene)
		#	temp_pr3file=Primer3Picker(temp_data[1], target_codon, out_num, temp_query[3]['ChrStop']-temp_query[3]['ChrStart']-75)
		#else:
		#	temp_pr3file=Primer3Picker(temp_data[1], target_codon, out_num, prod_length)
		temp_pr3file=Primer3Picker(temp_data[1], target_codon, out_num, prod_length)
		temp_txtfile=PrimerFileComp(temp_pr3file[0])
		#have subloop to record each line into master outfile
		subfile = open(temp_txtfile, 'rU')
		subfile.readline()
		for line2 in subfile:
			ID, Sense, Antisense, SenseStart, AntisenseStart, ProdLength= line2.split('\t')
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                     (Species,Strain,Gene,outflank, Sense, Antisense,interest_regions.rstrip('\n')))
		subfile.close()
	outfile.close()
	infile.close()
