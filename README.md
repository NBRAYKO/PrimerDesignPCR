PrimerDesignPCR_
================

Python scripts to design and check primer stats for PCR 

'''===============

PickPCR has functions to look up and download fasta files and design primers for particular regions using Primer3. The input file for the main function PrimerParser has to be modeled after cip_test_pr3.txt and contain organism, strain (optional), gene and target mutation site (optional)

Target mutation site has to specify at least a codon number. The exact name of the amino-acid can be in any form (e.g., Glu-84, G84 etc.). Multiple codon numbers can be entered if separated by comma. If >1 codon is entered as a target site, program picks the midpoint (so, if interested in sites too far apart, make it a separate query)

INPUT: A file name with organism name, strain, gene, target region (see cip_test_pr3), product length and the number of primers to be returned for each line) RETURN: Adds primers to the input (see cip_test.txt ), and individual files with more details (see grlAsite82_pr3_output.txt)

Notes on how to set up the dependent scripts for Emboss and Primer3 (this is srsly byzantine): 1. INSTALL PRIMER3 using brew https://github.com/Homebrew/homebrew-science 2. INSTALL EMBOSS (need the eprimer3 interface script in addition to primer3?), use brew install https://raw.github.com/Homebrew/homebrew-science/master/emboss.rb 3. use the biopython Primer3Commandline interface to run things from python 4. Note: to get things running, need to have the primer3_config dir copied in the current working dir 5. Note: also, need to download the old version of primer3 1.xx, and use terminal to set the emboss path to recognize this as teh primer3_core file'\ Type $ export EMBOSS_PRIMER3_CORE="/usr/local/Cellar/primer3/1.1.4/bin/primer3_core" 6. finally, this python script raises exceptions File "//anaconda/lib/python2.7/site-packages/Bio/Application/init.py", line 444, in call Despite the fact the program runs fine in bash; however, this environmental variable is not stored when running python, so have to enter the command os.environ['EMBOSS_PRIMER3_CORE'] before each run

Finally, had to modify line 44 in '/Bio/Application/init.py' so it does not raise errors based on Eprimer outpit codes, as it was not interpreting them correctly '''

"""===============

CheckPCR computes parameters to predict the success of PCR (GC content, dimers, reaction temperature etc.).

The porgram does not design primers, but provides a fast and efficient way of comparing detailed statistics of the PCR reactions for proposed primers

To design primers, use the PickPCR module, which uses a custom implementation of the EMBOSS eprimer3 program.

INPUT: A file with organism name, strain, gene, sense and antisense primer, target region RETURN: Text tab delimited file containing: Bind, Outflanks, Tempreature, GC ratio, AA sequence, product sequence etc. (see PCRresult_cip_test.txt)

Source: Existing code from Y. Benita was expanded to download reference sequence using Entrez, loop over a simpler input file... in addition, the output formatting was tweaked in several small ways (return the AA sequence of interest, show amplicon size etc.) Original Benita code is available http://xavierlab.mgh.harvard.edu/benita/PCRSupInfo/ (Accessed in December 2013)

For testing: Put all the file os.chdir('/Users/.../PrimerDesignPCR') input_filename = 'cip_test.txt' """
