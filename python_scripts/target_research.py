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

#################################################################################
# Script for RACEseq analysis 
# GOAL: extract reads that match to mRNA
#################################################################################


###import needed python and biopython modules
import regex
import os, sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
from collections import defaultdict

# the data directory
readfn = sys.argv[1] # deduplicated read1 fastq files
outfn = sys.argv[2] #to adjust directory where output files will be saved

DICT = dict([arg.split('=') for arg in sys.argv[4:]])

print DICT
	
###Main program
def main() :
	print 'Gene\tAll\tTarget'
	COUNT = dict.fromkeys(DICT, 0)
	Infile = gzip.open(readfn,'r')
	outfile = gzip.open(outfn,'w')
	Countall = 0
	global Count
	global Nb_mmR1
	Nb_mmR1 = int(sys.argv[3]) # number of other accepted nucleotides
	Count = 0
	for title, seq, qual in FastqGeneralIterator(Infile):
		Countall += 1
		global key
		for key in DICT:
			if regex.findall('(%s){s<=%i}'%(DICT[key], Nb_mmR1), seq): # 1 mm acccepted
				outfile.write("@%s\t%s\n%s\n+\n%s\n" % (title, key, seq, qual))
				COUNT[key] += 1
	outfile.close()
	for key in COUNT:
		print '%s\t%i\t%i' % (key, Countall, COUNT[key])

		
if __name__ == "__main__" :
    main()


