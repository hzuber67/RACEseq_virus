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
# GOALS: - Research of potential long A tail (> 10nt) in the read 2 that didn't  map target sequence. 
#		- Research of potential other modification (U or C, G) on these tails  
#		- Reverse complement read 2 
###############################################################################################################

###import needed python and biopython modules
import re, regex
import os, sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
import gzip

####################################################
#function to extract polyA extension
####################################################
def extractpolyA(SearchStringx, Category):
	global Count2
	global Sequence
	global title
	Result=regex.search(SearchStringx, Sequence)
	PolyA=Result.group(1) + Result.group(2) + Result.group(3)
	PolyAsize=len(PolyA)
	Modification = ''
	Modificationsize=0
	Extensionsize= PolyAsize + Modificationsize
	Extension= ''.join([Modification, PolyA])
	PolyA = Seq(PolyA).reverse_complement()
	Extension = Seq(Extension).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (title.split('\t')[0], title.split('\t')[-2], PolyA,PolyAsize, Modification, Modificationsize, Extension, Extensionsize, Category, title.split('\t')[-1]))
	Count2 += 1	
	
####################################################
#function to extract extension
####################################################
def extract(SearchStringx, Category):
	global Count2
	global Sequence
	global title
	Result=regex.search(SearchStringx, Sequence)
	PolyA=Result.group(2) + Result.group(3) + Result.group(4)
	PolyAsize=len(PolyA)
	Modification = Result.group(1)
	Modificationsize=len(Modification)
	Extensionsize= PolyAsize + Modificationsize
	Extension= ''.join([Modification, PolyA])
	PolyA = Seq(PolyA).reverse_complement()
	Extension = Seq(Extension).reverse_complement()
	Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (title.split('\t')[0], title.split('\t')[-2], PolyA,PolyAsize, Modification, Modificationsize, Extension, Extensionsize, Category, title.split('\t')[-1]))
	Count2 += 1	
	
####################################################
#function to extract polyX extension
####################################################
def extractpolyX(SearchStringx, Category):
	global Count2
	global Sequence
	global title
	Result=regex.search(SearchStringx, Sequence)
	PolyA=''
	PolyAsize=0
	Modification = Result.group(1) + Result.group(2) + Result.group(3)
	Modificationsize=len(Modification)
	Extensionsize= PolyAsize + Modificationsize
	Extension= ''.join([Modification, PolyA])
	PolyA = Seq(PolyA).reverse_complement()
	Extension = Seq(Extension).reverse_complement()
	Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (title.split('\t')[0], title.split('\t')[-2], PolyA,PolyAsize, Modification, Modificationsize, Extension, Extensionsize, Category, title.split('\t')[-1]))
	Count2 += 1	
	
		
####################################################
###Main program
####################################################

def main() :
	###define regex regular expression
	global Nb_mm
	Nb_mm = int(sys.argv[4]) # number of other nucleotides accepted in the polyA
	##at least 31 by accepting 1 other nucleotides
	A_SearchStringpolyA = '^([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	B_SearchStringpolyAandU = '(^[A]{1,29})([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	C_SearchStringpolyAandG = '(^[C]{1,29})([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	D_SearchStringpolyAandC = '(^[G]{1,29})([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	
	E_SearchStringpolyU = '(^A{3})([A]{25,}){s<=%i}([A]{3})(.*$)' % (Nb_mm) #polyU only
	F_SearchStringpolyG = '(^C{3})([C]{25,}){s<=%i}([C]{3})(.*$)' % (Nb_mm) #polyG only
	G_SearchStringpolyC = '(^G{3})([G]{25,}){s<=%i}([C]{3})(.*$)' % (Nb_mm) #polyC only
	
	H_SearchStringpolyAandhetero= '(^[TACG]{1,29}?)([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm) 

	
	##at least 8 by accepting 0 other nucleotides
	A_SearchStringpolyA_2 = '^([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	B_SearchStringpolyAandU_2 = '(^[A]{1,29})([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	C_SearchStringpolyAandG_2 = '(^[C]{1,29})([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	D_SearchStringpolyAandC_2 = '(^[G]{1,29})([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	
	E_SearchStringpolyU_2 = '(^A{3})([A]{5}){s<=%i}([A]{2,})(.*$)' % (Nb_mm) #polyU only
	F_SearchStringpolyG_2 = '(^C{3})([C]{5}){s<=%i}([C]{2,})(.*$)' % (Nb_mm) #polyG only
	G_SearchStringpolyC_2 = '(^G{3})([G]{5}){s<=%i}([C]{2,})(.*$)' % (Nb_mm) #polyC only
	
	H_SearchStringpolyAandhetero_2= '(^[TACG]{1,29}?)([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm) #polyA + heteropolymeric
	
	
	global j
	global File
	global Sequence
	global title
	global Count2
	global CountA
	global CountB
	global CountC
	global CountD
	global CountE
	global CountF
	global CountG
	global CountH
	#global CountI
	global CountJ
	global Count3
	Count1= 0
	Count2= 0
	Count3= 0
	CountA= 0
	CountB= 0
	CountC= 0
	CountD= 0
	CountE= 0
	CountF= 0
	CountG= 0
	CountH= 0
	#CountI= 0
	CountJ= 0

	#external argument
	readfn = sys.argv[1]
	outfn = sys.argv[2]
	unfn = sys.argv[3]
	Infile = gzip.open(readfn,'r')
	File = open(outfn,'w')
	Other = open(unfn,'w')
		
	File.write("Read.ID\tTarget\tPolyA\tPolyA.size\tModification\tModification.size\tExtension\tExtension.size\tClassification\tRandom\n") 
	for title, Sequence, qual in FastqGeneralIterator(Infile) :
		Count1+= 1
		if regex.findall(A_SearchStringpolyA, Sequence):#1 extract polyA extension without other modification
			extractpolyA(A_SearchStringpolyA, 'A_polyA')
			CountA+=1
		elif regex.findall(B_SearchStringpolyAandU, Sequence):#2 extract polyA + U 
			extract(B_SearchStringpolyAandU, 'B_polyAandU')	
			CountB+=1
		elif regex.findall(C_SearchStringpolyAandG, Sequence):#3 extract polyA extension with G
			extract(C_SearchStringpolyAandG, 'C_polyAandG')	
			CountC+=1
		elif regex.findall(D_SearchStringpolyAandC, Sequence):#4extract polyA extension with C
			extract(D_SearchStringpolyAandC, 'D_polyAandC')
			CountD+=1
		elif regex.findall(E_SearchStringpolyU, Sequence): #5 extract polyA extension with AU modification
			extractpolyX(E_SearchStringpolyU, 'E_polyU')
			CountE+=1
		elif regex.findall(F_SearchStringpolyG, Sequence):#1 extract polyA extension without other modification
			extractpolyX(F_SearchStringpolyG, 'F_polyG')
			CountF+=1
		elif regex.findall(G_SearchStringpolyC, Sequence):#2 extract polyA + U 
			extractpolyX(G_SearchStringpolyC, 'G_polyC')	
			CountG+=1
		elif regex.findall(H_SearchStringpolyAandhetero, Sequence):#3 extract polyA extension with G
			extract(H_SearchStringpolyAandhetero, 'H_polyAandhetero')	
			CountH+=1	
		#elif regex.findall(I_SearchStringheteropolyA, Sequence):#4extract polyA extension with C
			#extractpolyA(I_SearchStringheteropolyA, 'I_heteropolyA')
			#CountI+=1
		elif regex.findall(A_SearchStringpolyA_2, Sequence):#1 extract polyA extension without other modification
			extractpolyA(A_SearchStringpolyA_2, 'A_polyA')
			CountA+=1
		elif regex.findall(B_SearchStringpolyAandU_2, Sequence):#2 extract polyA + U 
			extract(B_SearchStringpolyAandU_2, 'B_polyAandU')	
			CountB+=1
		elif regex.findall(C_SearchStringpolyAandG_2, Sequence):#3 extract polyA extension with G
			extract(C_SearchStringpolyAandG_2, 'C_polyAandG')	
			CountC+=1
		elif regex.findall(D_SearchStringpolyAandC_2, Sequence):#4extract polyA extension with C
			extract(D_SearchStringpolyAandC_2, 'D_polyAandC')
			CountD+=1
		elif regex.findall(E_SearchStringpolyU_2, Sequence): #5 extract polyA extension with AU modification
			extractpolyX(E_SearchStringpolyU_2, 'E_polyU')
			CountE+=1
		elif regex.findall(F_SearchStringpolyG_2, Sequence):#1 extract polyA extension without other modification
			extractpolyX(F_SearchStringpolyG_2, 'F_polyG')
			CountF+=1
		elif regex.findall(G_SearchStringpolyC_2, Sequence):#2 extract polyA + U 
			extractpolyX(G_SearchStringpolyC_2, 'G_polyC')	
			CountG+=1
		elif regex.findall(H_SearchStringpolyAandhetero_2, Sequence):#3 extract polyA extension with G
			extract(H_SearchStringpolyAandhetero_2, 'H_polyAandhetero')	
			CountH+=1
		#elif regex.findall(I_SearchStringheteropolyA_2, Sequence):#4extract polyA extension with C
			#extractpolyA(I_SearchStringheteropolyA_2, 'I_heteropolyA')
			#CountI+=1

		
		else:
			Count3+=1
			rev_seq = Seq(Sequence).reverse_complement()
			Other.write("@%s\n%s\n+\n%s\n" % (title, rev_seq, qual))
	File.close()
	Other.close()
	print "total\tanalyzed\tundetermined\tA_polyA\tB_polyAandU\tC_polyAandG\tD_polyAandC\tE_polyU\tF_polyG\tG_polyC\tH_polyAandhetero\tJ_no_tail" 
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (Count1, Count2, Count3, CountA, CountB, CountC, CountD, CountE, CountF, CountG, CountH, CountJ)

	
if __name__ == "__main__" :
    main()

	