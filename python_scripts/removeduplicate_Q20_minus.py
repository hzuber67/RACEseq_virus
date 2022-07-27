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
# GOAL: deduplication of sequences with identical nucleotides in the 1st to 15th cycles
#   (randomized bases in 3 prime adapter) and in the 20 and 50th cycles of read 2 
# Filter low quality reads: filter out reads with at least one nucleotide < Q10 in the 1st to 15th cycles
#   or in the 20 and 50th cycles of read 2 
#################################################################################

import re
import gzip
import os, sys
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools


badstring = "[!\"#$%&'()*+,-./012345]" ##remove sequence with low quality nucleotide in the 30 first nucleotides of read2 (Q20)


#definition of the main function
def main() :
	print "Total\tBad_quality\tDeduplicated"
	readfn = sys.argv[1] #read 1 fastq files
	readfn2 = sys.argv[2] #read 1 fastq files
	outname =  sys.argv[3] #directory where deduplicated read 1 files will be saved
	outname2 =  sys.argv[4] #directory where deduplicated read 1 files will be saved
	Unique_seqs=set()
	Infile1 = gzip.open(readfn, "r")
	outfile1 = gzip.open(outname,"w")
	Infile2 = gzip.open(readfn2, "r")
	outfile2 = gzip.open(outname2,"w")
	Countall=0
	Count=0
	Countbad=0
	Unique_seqs=set()
	f_iter = FastqGeneralIterator(Infile1)
	r_iter = FastqGeneralIterator(Infile2)
	for (title, seq, qual), (title2, seq2, qual2) in itertools.izip(f_iter, r_iter):
		assert title.split(None,1)[0] == title2.split(None,1)[0]
		seq2b =seq2[0:8]
		qualc =str(qual2[0:8]+qual2[13:43])
		Countall+=1
		if not re.compile(badstring).search(qualc):
			if str(seq2b) not in Unique_seqs:
				Unique_seqs.add(str(seq2b))
				outfile1.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
				outfile2.write("@%s\n%s\n+\n%s\n" % (title2, seq2, qual2))
				Count+=1
		else:
			Countbad+=1
	print "%i\t%i\t%i" % (Countall, Countbad, Count)	
	outfile1.close()
	outfile2.close()
	
	
if __name__ == "__main__" :
    main()


