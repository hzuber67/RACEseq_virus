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
#function to extract polyU extension
####################################################
def extractpolyU(SearchStringx, Category):
	global Count2
	global Sequence
	global title
	Result=regex.search(SearchStringx, Sequence)
	PolyU=Result.group(1) + Result.group(2) + Result.group(3)
	PolyUsize=len(PolyU)
	Modification = ''
	Modificationsize=0
	Extensionsize= PolyUsize + Modificationsize
	Extension= ''.join([Modification, PolyU])
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (title.split('\t')[0], title.split('\t')[-2], PolyU,PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, title.split('\t')[-1]))
	Count2 += 1	
	
####################################################
#function to extract extension
####################################################
def extract(SearchStringx, Category):
	global Count2
	global Sequence
	global title
	Result=regex.search(SearchStringx, Sequence)
	PolyU=Result.group(2) + Result.group(3) + Result.group(4)
	PolyUsize=len(PolyU)
	Modification = Result.group(1)
	Modificationsize=len(Modification)
	Extensionsize= PolyUsize + Modificationsize
	Extension= ''.join([Modification, PolyU])
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (title.split('\t')[0], title.split('\t')[-2], PolyU,PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, title.split('\t')[-1]))
	Count2 += 1	
	
####################################################
#function to extract polyX extension
####################################################
def extractpolyX(SearchStringx, Category):
	global Count2
	global Sequence
	global title
	Result=regex.search(SearchStringx, Sequence)
	PolyU=''
	PolyUsize=0
	Modification = Result.group(1) + Result.group(2) + Result.group(3)
	Modificationsize=len(Modification)
	Extensionsize= PolyUsize + Modificationsize
	Extension= ''.join([Modification, PolyU])
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (title.split('\t')[0], title.split('\t')[-2], PolyU,PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, title.split('\t')[-1]))
	Count2 += 1	
	
		
####################################################
###Main program
####################################################

def main() :
	###define regex regular expression
	global Nb_mm
	Nb_mm = int(sys.argv[4]) # number of other nucleotides accepted in the polyU
	##at least 31 by accepting 1 other nucleotides
	A_SearchStringpolyU = '^([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	B_SearchStringpolyUandA = '(^[A]{1,29})([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	C_SearchStringpolyUandG = '(^[G]{1,29})([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	D_SearchStringpolyUandC = '(^[C]{1,29})([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	
	E_SearchStringpolyA = '(^A{3})([A]{25,}){s<=%i}([A]{3})(.*$)' % (Nb_mm) #polyU only
	F_SearchStringpolyG = '(^G{3})([G]{25,}){s<=%i}([G]{3})(.*$)' % (Nb_mm) #polyG only
	G_SearchStringpolyC = '(^C{3})([C]{25,}){s<=%i}([C]{3})(.*$)' % (Nb_mm) #polyC only
	
	H_SearchStringpolyUandhetero= '(^[TACG]{1,29}?)([T]{3})([T]{25,}){s<=%i}([T]{3})(.*$)' % (Nb_mm) 

	
	##at least 8 by accepting 0 other nucleotides
	A_SearchStringpolyU_2 = '^([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	B_SearchStringpolyUandA_2 = '(^[A]{1,29})([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	C_SearchStringpolyUandG_2 = '(^[G]{1,29})([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	D_SearchStringpolyUandC_2 = '(^[C]{1,29})([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm)
	
	E_SearchStringpolyA_2 = '(^A{3})([A]{5}){s<=%i}([A]{2,})(.*$)' % (Nb_mm) #polyU only
	F_SearchStringpolyG_2 = '(^G{3})([G]{5}){s<=%i}([G]{2,})(.*$)' % (Nb_mm) #polyG only
	G_SearchStringpolyC_2 = '(^C{3})([C]{5}){s<=%i}([C]{2,})(.*$)' % (Nb_mm) #polyC only
	
	H_SearchStringpolyUandhetero_2= '(^[TACG]{1,29}?)([T]{3})([T]{2,}){s<=%i}([T]{3})(.*$)' % (Nb_mm) #polyU + heteropolymeric
	
	
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
		
	File.write("Read.ID\tTarget\tPolyU\tPolyU.size\tModification\tModification.size\tExtension\tExtension.size\tClassification\tRandom\n") 
	for title, Sequence, qual in FastqGeneralIterator(Infile) :
		Count1+= 1
		if regex.findall(A_SearchStringpolyU, Sequence):#1 extract polyU extension without other modification
			extractpolyU(A_SearchStringpolyU, 'A_polyU')
			CountA+=1
		elif regex.findall(B_SearchStringpolyUandA, Sequence):#2 extract polyU + U 
			extract(B_SearchStringpolyUandA, 'B_polyUandA')	
			CountB+=1
		elif regex.findall(C_SearchStringpolyUandG, Sequence):#3 extract polyU extension with G
			extract(C_SearchStringpolyUandG, 'C_polyUandG')	
			CountC+=1
		elif regex.findall(D_SearchStringpolyUandC, Sequence):#4extract polyU extension with C
			extract(D_SearchStringpolyUandC, 'D_polyUandC')
			CountD+=1
		elif regex.findall(E_SearchStringpolyA, Sequence): #5 extract polyU extension with AU modification
			extractpolyX(E_SearchStringpolyA, 'E_polyA')
			CountE+=1
		elif regex.findall(F_SearchStringpolyG, Sequence):#1 extract polyU extension without other modification
			extractpolyX(F_SearchStringpolyG, 'F_polyG')
			CountF+=1
		elif regex.findall(G_SearchStringpolyC, Sequence):#2 extract polyU + U 
			extractpolyX(G_SearchStringpolyC, 'G_polyC')	
			CountG+=1
		elif regex.findall(H_SearchStringpolyUandhetero, Sequence):#3 extract polyU extension with G
			extract(H_SearchStringpolyUandhetero, 'H_polyUandhetero')	
			CountH+=1	
		#elif regex.findall(I_SearchStringheteropolyU, Sequence):#4extract polyU extension with C
			#extractpolyU(I_SearchStringheteropolyU, 'I_heteropolyU')
			#CountI+=1
		elif regex.findall(A_SearchStringpolyU_2, Sequence):#1 extract polyU extension without other modification
			extractpolyU(A_SearchStringpolyU_2, 'A_polyU')
			CountA+=1
		elif regex.findall(B_SearchStringpolyUandA_2, Sequence):#2 extract polyU + U 
			extract(B_SearchStringpolyUandA_2, 'B_polyUandA')	
			CountB+=1
		elif regex.findall(C_SearchStringpolyUandG_2, Sequence):#3 extract polyU extension with G
			extract(C_SearchStringpolyUandG_2, 'C_polyUandG')	
			CountC+=1
		elif regex.findall(D_SearchStringpolyUandC_2, Sequence):#4extract polyU extension with C
			extract(D_SearchStringpolyUandC_2, 'D_polyUandC')
			CountD+=1
		elif regex.findall(E_SearchStringpolyA_2, Sequence): #5 extract polyU extension with AU modification
			extractpolyX(E_SearchStringpolyA_2, 'E_polyA')
			CountE+=1
		elif regex.findall(F_SearchStringpolyG_2, Sequence):#1 extract polyU extension without other modification
			extractpolyX(F_SearchStringpolyG_2, 'F_polyG')
			CountF+=1
		elif regex.findall(G_SearchStringpolyC_2, Sequence):#2 extract polyU + U 
			extractpolyX(G_SearchStringpolyC_2, 'G_polyC')	
			CountG+=1
		elif regex.findall(H_SearchStringpolyUandhetero_2, Sequence):#3 extract polyU extension with G
			extract(H_SearchStringpolyUandhetero_2, 'H_polyUandhetero')	
			CountH+=1
		#elif regex.findall(I_SearchStringheteropolyU_2, Sequence):#4extract polyU extension with C
			#extractpolyU(I_SearchStringheteropolyU_2, 'I_heteropolyU')
			#CountI+=1

		
		else:
			Count3+=1
			#rev_seq = Seq(Sequence).reverse_complement()
			Other.write("@%s\n%s\n+\n%s\n" % (title, Sequence, qual))
	File.close()
	Other.close()
	print "total\tanalyzed\tundetermined\tA_polyU\tB_polyUandU\tC_polyUandG\tD_polyUandC\tE_polyA\tF_polyG\tG_polyC\tH_polyUandhetero\tJ_no_tail" 
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (Count1, Count2, Count3, CountA, CountB, CountC, CountD, CountE, CountF, CountG, CountH, CountJ)

	
if __name__ == "__main__" :
    main()

	