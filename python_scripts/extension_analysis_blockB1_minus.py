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
# GOALS: - Analysis of 3' tail found after read2 mapping. 
#		- Research of potential other modification (U or C, G) on these tails 
#		- Reverse complement read 2 
###############################################################################################################

###import needed python and biopython modules
import re, regex
import os, sys
from Bio.Seq import Seq
from fractions import Fraction
from math import *
import gzip

####################################################
#function in extension length>0 
####################################################

def extension(x):
	global Infile2
	global Count2
	global Count3
	global Extension
	global Done
	global Name
	global Gene
	global Random
	global Position
	global Sequence
	global CountA
	global CountB
	global CountC
	global CountD
	global CountE
	global CountF
	global CountG
	global CountH
	#global CountI
	global Nb_mm
	Nb_mm = int(sys.argv[5]) # number of other nucleotides accepted in the polyA
	Done = 'FALSE'
	if x > 7:
		A_SearchStringpolyU = '(^T)(T)(T+){s<=%i}($)' % (Nb_mm) #polyA only, 1 mm in polyA accepted
		B_SearchStringpolyUandA = '(^[A]{1,29})(T)(T+){s<=%i}($)' % (Nb_mm) #polyA + U, 1 mm in polyA accepted
		C_SearchStringpolyUandG = '(^[G]{1,29})(T)(T+){s<=%i}($)' % (Nb_mm) #polyA + G, 1 mm in polyA accepted => we should have no G modification
		D_SearchStringpolyUandC = '(^[C]{1,29})(T)(T+){s<=%i}($)' % (Nb_mm) #polyA + C, 1 mm in polyA accepted 
		
		E_SearchStringpolyA = '(^A)(A+){s<=%i}($)' % (Nb_mm) #polyU only, 1 mm in polyU accepted
		F_SearchStringpolyG = '(^G)(G+){s<=%i}($)' % (Nb_mm) #polyG only, 1 mm in polyG accepted => we should have no G modification
		G_SearchStringpolyC = '(^C)(C+){s<=%i}($)' % (Nb_mm) #polyC only, 1 mm in polyC accepted
		
		H_SearchStringpolyUandhetero= '(^[TACG]{1,29}?)(T)(T+){s<=%i}($)' % (Nb_mm) #polyA + heteropolymeric,	
		#I_SearchStringheteropolyA = '(^TTT)(T+){s<=%i}($)' % (3) #hetero with up to 3 of other nucleotide
	
		if regex.findall(A_SearchStringpolyU, Extension):
			extractpolyA(A_SearchStringpolyU, 'A_polyU')
			Done = 'TRUE'
			CountA+=1
		elif regex.findall(B_SearchStringpolyUandA, Extension):
			extract(B_SearchStringpolyUandA, 'B_polyUandA')
			Done = 'TRUE'
			CountB+=1
		elif regex.findall(C_SearchStringpolyUandG, Extension):
			extract(C_SearchStringpolyUandG,'C_polyUandG')
			Done = 'TRUE'
			CountC+=1
		elif regex.findall(D_SearchStringpolyUandC, Extension):
			extract(D_SearchStringpolyUandC,'D_polyUandC')
			Done = 'TRUE'
			CountD+=1
		elif regex.findall(E_SearchStringpolyA, Extension):
			extractpolyX(E_SearchStringpolyA,'E_polyU')
			Done = 'TRUE'
			CountE+=1
		elif regex.findall(F_SearchStringpolyG, Extension):
			extractpolyX(F_SearchStringpolyG,'F_polyG')
			Done = 'TRUE'
			CountF+=1
		elif regex.findall(G_SearchStringpolyC, Extension):
			extractpolyX(G_SearchStringpolyC,'G_polyC')
			Done = 'TRUE'
			CountG+=1
		elif regex.findall(H_SearchStringpolyUandhetero, Extension):
			extracthetero(H_SearchStringpolyUandhetero,'H_polyUandhetero')
			Done = 'TRUE'
	elif (8 > x > 1):
		A_SearchStringpolyU = '(^T)(T)(T*){s<=0}($)' #polyA only, 1 mm in polyA accepted
		B_SearchStringpolyUandA = '(^[A]{1,29})(T)(T*){s<=0}($)' #polyA + U, 1 mm in polyA accepted
		C_SearchStringpolyUandG = '(^[G]{1,29})(T)(T*){s<=0}($)' #polyA + G, 1 mm in polyA accepted
		D_SearchStringpolyUandC = '(^[C]{1,29})(T)(T*){s<=0}($)' #polyA + C, 1 mm in polyA accepted
		
		E_SearchStringpolyA = '(^A)(A*){s<=0}($)' #polyU only, 1 mm in polyU accepted
		F_SearchStringpolyG = '(^G)(G*){s<=0}($)' #polyG only, 1 mm in polyG accepted
		G_SearchStringpolyC = '(^C)(C*){s<=0}($)' #polyC only, 1 mm in polyC accepted
		
		H_SearchStringpolyUandhetero= '(^[TACG]{1,29}?)(T)(T*){s<=0}($)' #polyA + heteropolymeric,	
		#I_SearchStringheteropolyA = '(^T)(T+){s<=%i}($)' % (2) #hetero with up to 20% of other nucleotide
		
		if regex.findall(A_SearchStringpolyU, Extension):
			extractpolyA(A_SearchStringpolyU, 'A_polyU')
			Done = 'TRUE'
			CountA+=1
		elif regex.findall(B_SearchStringpolyUandA, Extension):
			extract(B_SearchStringpolyUandA, 'B_polyUandA')
			Done = 'TRUE'
			CountB+=1
		elif regex.findall(C_SearchStringpolyUandG, Extension):
			extract(C_SearchStringpolyUandG,'C_polyUandG')
			Done = 'TRUE'
			CountC+=1
		elif regex.findall(D_SearchStringpolyUandC, Extension):
			extract(D_SearchStringpolyUandC,'D_polyUandC')
			Done = 'TRUE'
			CountD+=1
		elif regex.findall(E_SearchStringpolyA, Extension):
			extractpolyX(E_SearchStringpolyA,'E_polyA')
			Done = 'TRUE'
			CountE+=1
		elif regex.findall(F_SearchStringpolyG, Extension):
			extractpolyX(F_SearchStringpolyG,'F_polyG')
			Done = 'TRUE'
			CountF+=1
		elif regex.findall(G_SearchStringpolyC, Extension):
			extractpolyX(G_SearchStringpolyC,'G_polyC')
			Done = 'TRUE'
			CountG+=1
		elif regex.findall(H_SearchStringpolyUandhetero, Extension):
			extracthetero(H_SearchStringpolyUandhetero,'H_polyUandhetero')
			Done = 'TRUE'
		#elif regex.findall(I_SearchStringheteropolyA, Extension):
			#extractpolyA(I_SearchStringheteropolyA,'I_heteropolyA')
			#Done = 'TRUE'
			#CountI+=1
	elif x ==1:
		A_SearchStringpolyU = '(^T)($)' #polyA only, 1 mm in polyA accepted
		E_SearchStringpolyA = '(^A)($)' #polyU only, 1 mm in polyU accepted
		F_SearchStringpolyG = '(^G)($)' #polyG only, 1 mm in polyG accepted
		G_SearchStringpolyC = '(^C)($)' #polyC only, 1 mm in polyC accepted
		if regex.findall(A_SearchStringpolyU, Extension):
			polyA_short(A_SearchStringpolyU, 'A_polyU')
			Done = 'TRUE'
			CountA+=1
		elif regex.findall(E_SearchStringpolyA, Extension):
			extract_short(E_SearchStringpolyA,'E_polyA')
			Done = 'TRUE'
			CountE+=1
		elif regex.findall(F_SearchStringpolyG, Extension):
			extract_short(F_SearchStringpolyG,'F_polyG')
			Done = 'TRUE'
			CountF+=1
		elif regex.findall(G_SearchStringpolyC, Extension):
			extract_short(G_SearchStringpolyC,'G_polyC')
			Done = 'TRUE'
			CountG+=1
	if Done == 'FALSE':
		title2 = ("%s\t%s\t%s") % (Name, Gene, Random)
		Sequence2 = Extension + Sequence
		Infile2.write("@%s\n%s\n+\n%s\n" % (title2, Sequence2, Sequence2))
		Count3 += 1

####################################################
#function to extract polyA extension
####################################################
def extractpolyA(SearchStringx, Category):
	global Count2
	global Extension
	global Extensionsize 
	global Name
	global Random
	global Gene
	global Position
	Result=regex.search(SearchStringx, Extension)
	PolyU=Result.group(1) + Result.group(2) + Result.group(3)
	PolyUsize=len(PolyU)
	Modification = ''
	Modificationsize=0
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	Count2 += 1				
				
####################################################
#function to extract extension
####################################################
def extract(SearchStringx, Category):
	global Count2
	global Extension
	global Extensionsize
	global Name
	global Random
	global Gene
	global Position
	Result=regex.search(SearchStringx, Extension)
	PolyU=Result.group(2) + Result.group(3)
	PolyUsize=len(PolyU)
	Modification = Result.group(1)
	Modificationsize=len(Modification)
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	Count2 += 1

####################################################
#function to extract polyX
####################################################
def extractpolyX(SearchStringx, Category):
	global Count2
	global Extension
	global Extensionsize
	global Name
	global Random
	global Gene
	global Position
	Result=regex.search(SearchStringx, Extension)
	PolyU=''
	PolyUsize=0
	Modification = Result.group(1) + Result.group(2)
	Modificationsize=len(Modification)
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	Count2 += 1

####################################################
#function to extract extension
####################################################
def polyA_short(SearchStringx, Category):
	global Count2
	global Extension
	global Extensionsize
	global Name
	global Random
	global Gene
	global Position
	Result=regex.search(SearchStringx, Extension)
	PolyU=Result.group(1)
	PolyUsize=len(PolyU)
	Modification = ''
	Modificationsize=0
	Extensionsize = PolyUsize + Modificationsize
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	Count2 += 1
	
####################################################
#function to extract extension
####################################################
def extract_short(SearchStringx, Category):
	global Count2
	global Extension
	global Extensionsize
	global Name
	global Random
	global Gene
	global Position
	Result=regex.search(SearchStringx, Extension)
	PolyU=''
	PolyUsize=0
	Modification = Result.group(1)
	Modificationsize=len(Modification)
	Extensionsize = PolyUsize + Modificationsize
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
	Count2 += 1
	
####################################################
#function to extract hetero extension
####################################################
def extracthetero(SearchStringx, Category):
	global Count2
	global CountH
	global Extension
	global Extensionsize
	global Name
	global Random
	global Gene
	global Position
	global Done
	Result=regex.search(SearchStringx, Extension)
	PolyU=Result.group(2) + Result.group(3)
	PolyUsize=len(PolyU)
	Modification = Result.group(1)
	Modification2 = Result.group(1)
	Modificationsize=len(Modification)
	Extensionsize = PolyUsize + Modificationsize
	#PolyU = Seq(PolyU).reverse_complement()
	#Extension = Seq(Extension).reverse_complement()
	#Modification = Seq(Modification).reverse_complement()
	j= ceil(Fraction(len(Modification2)*2,3))#70%
	Nrich = '(.*A.*){%i}|(.*T.*){%i}|(.*C.*){%i}|(.*G.*){%i}' % (j, j, j, j)
	if regex.findall(Nrich, Modification2):
		CountH+=1
		File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
		File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
		Count2 += 1
	else:
		Done="FALSE"
	
####################################################
###Main program
####################################################

def main() :
	SearchLine = '^(.+)\t1:N:0:\d*\t(.+)\t([ATCGN]{8})\t(\d+)\t([ATCGN]+)\t([ATCG]{0,30})\t(\d{1,2})$'
	Searchstring0 = '(^T+).+$'
	global File
	global File2
	global Count2
	global Count3
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
	
	global Extension
	global Extensionsize
	global Sequence
	global Name
	global Random
	global Gene
	global Position
	global Infile2
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
	readfn = sys.argv[1]
	outfn = sys.argv[2]
	outfn2 = sys.argv[3]
	readfn2 = sys.argv[4]
	Infile = open(readfn,'r')
	Infile2 = gzip.open(readfn2,'a')
	File = open(outfn,'w')
	File2 = open(outfn2,'w')
	File.write("Read.ID\tTarget\tPolyU\tPolyU.size\tModification\tModification.size\tExtension\tExtension.size\tClassification\tRandom\n")
	File2.write("Read.ID\tTarget\tPosition\tPolyU\tPolyU.size\tModification\tModification.size\tExtension\tExtension.size\tClassification\tRandom\n") 
	for Line in Infile:
		Count1 += 1 
		if re.compile(SearchLine).search(Line):
			Result=re.search(SearchLine, Line) #parse line
			Name = Result.group(1)
			Gene= Result.group(2)
			Random = Result.group(3)
			Position= Result.group(4)
			Sequence = Result.group(5)
			Extension = Result.group(6)
			Extensionsize = int(Result.group(7))
			if re.compile(Searchstring0).search(Sequence):
				Result2=re.search(Searchstring0, Sequence)
				Add = Result2.group(1)
				Extension = Extension + Add
				Extensionsize = len(Extension)
			if Extensionsize==0:
				Category = 'J_no_tail'
				PolyU=''
				PolyUsize=0
				Modification = ''
				Modificationsize=0
				Extension = Seq(Extension).reverse_complement()
				File.write("%s\t%s\t%s\t%i\t%s\t%i\t%s\t%i\t%s\t%s\n" % (Name, Gene, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
				File2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"  % (Name, Gene, Position, PolyU, PolyUsize, Modification, Modificationsize, Extension, Extensionsize, Category, Random))
				Count2 += 1 
				CountJ += 1
			elif Extensionsize > 0:
				extension(Extensionsize)
	print "total\tanalyzed\tundetermined\tA_polyU\tB_polyUandA\tC_polyUandG\tD_polyUandC\tE_polyA\tF_polyG\tG_polyC\tH_polyUandhetero\tJ_no_tail" 
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (Count1, Count2, Count3, CountA, CountB, CountC, CountD, CountE, CountF, CountG, CountH, CountJ)

if __name__ == "__main__" :
    main()
