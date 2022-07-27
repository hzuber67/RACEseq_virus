#!/usr/bin/env python

# Copyright (c) 2016-2020 Institute of Plant Molecular Biology, CNRS
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Helene Zuber <helene.zuber@ibmp-cnrs.unistra.fr>

###############################################################################################################
# Script for RACEseq analysis 
# GOAL: - determine the 3 prime extremities of the target RNAs and identify untemplated tail
#				1) reads 2 are mapped to the reference sequence of the corresponding Arabidopsis thaliana mRNA 
#       		2) to identify reads with untemplated tails and map their 3 prime end position, the sequences 
#				   of the unmatched reads 2 are  successively trimmed from their 3 prime end
###############################################################################################################

###required modules: 
import os, sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex, re
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

########################################################################################################################
#function for mapping

def mapping():
	#give file names
	global Counttotal 
	global Countmapped 
	global Countunmapped
	Counttotal = 0
	Countmapped = 0
	Countunmapped = 0
	for title, seq, qual in FastqGeneralIterator(readfn_file):
		Done = 'FALSE'
		Counttotal += 1
		Gene = title.split('\t')[-2]
		Target_seq = DICT[Gene]
		Target_seq_rev = str(Seq(Target_seq).reverse_complement())
		Length = len(Target_seq_rev)
		stringfilter = '(^T)(T){30}'
		if regex.compile(stringfilter).search(seq):
			other_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
			Countunmapped +=1
			Done = 'TRUE'
		else:
			for i in range(0, 31): #from 0 to 30
				x= i 
				objend =seq[x:x+4]
				obj =seq[x+4:x+30]
				if len(seq[x:x+30]) < 30:
					Done = 'FALSE'
					break
				else:
					stringtest= '(%s){e<=0}(%s){i<=2,d<=2,s<=4,e<=6}'% (objend,obj)
					if regex.compile(stringtest).search(Target_seq_rev):
						for k in range(0, 6):
							string= '(%s){e<=0}(%s){i<=2,d<=2,s<=4,e<=%i}'% (objend,obj,k)
							if regex.compile(string).search(Target_seq_rev):
								research = regex.compile(string, regex.BESTMATCH)
								Countmapped += 1
								for match in research.finditer(Target_seq_rev):
									output_file.write("%s\t%s\t%s\t%s\t%i\n" % (title, Length-match.start(), seq[x:], seq[:x], x))
									Done = 'TRUE'
							if Done == 'TRUE':
								break
				if Done == 'TRUE':
						break
					
			if Done == 'FALSE':
				other_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
				Countunmapped +=1
				Done = 'TRUE'
	print '%s\t%i\t%i\t%i\t%i\n' % (readfn.split('/')[-1], x, Counttotal, Countmapped, Countunmapped )





########################################################################################################################
#main program
########################################################################################################################
### header for files

def main() :
	global readfn
	global readfn_file
	global output_file
	global DICT
	global other_file
	readfn = sys.argv[1] #file with trimmed read2
	output_fn = sys.argv[2]
	other = sys.argv[3]
	filedir = sys.argv[4]

	# open file with reference mRNA sequences
	#file with the precursor sequence
	DICT = {}
	Target_sequences = SeqIO.parse(open(filedir),'fasta')
	for fasta in Target_sequences:
		DICT[fasta.id]=str(fasta.seq)
	readfn_file=gzip.open(readfn,'r')
	output_file=open(output_fn,'w')
	print 'File\tPosition\ttotal\tmapped\tunmapped' # write header for counting files
	other_file = gzip.open(other,'w')
	mapping()

if __name__ == "__main__" :
    main()

