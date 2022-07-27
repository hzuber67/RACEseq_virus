#!/bin/sh
#
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


#This script allows for RACE-seq analyses. This is the version for mRNA analyses. Input files fastq files generated from MiSeq
#Require tools:
#python v2.7
#biopython
#regex

source config.conf


########################################################
###Configuration part that doesn't need to be adjusted:
########################################################
### Data directories:
Outdir_d=$Outdir"dedupl/" #output directory for deduplicated reads
Outdir_dR1=$Outdir_d"R1/" #output directory for deduplicated read 1
Outdir_dR2=$Outdir_d"R2/" #output directory for deduplicated read 2
Outdir_targetR1=$Outdir_d"target_R1/" #output directory for target read 1
Outdir_targetR2=$Outdir_d"target_R2_delim/" #output directory for target read 2 with delimiter
Outdir_trimmedR2=$Outdir_d"trimmed_R2/" #output directory for trimmed read 2 

Outdir_ext=$Outdir"R2_mapping/" #output directory for extension analysis
Outdir_other=$Outdir_ext"unmapped/" #output directory for non mapped read 2
Outdir_results=$Outdir_ext"mapped/" #output directory for extension analysis results

Outdir_final=$Outdir"process_results/" #final output directory
Outdir_mapped=$Outdir_final"mapped_reads/" #results read2 that could map to the target sequence after trimming
Outdir_mapped_PA=$Outdir_final"position/" #results read2 that could map to the target sequence after trimming
Outdir_combined=$Outdir_final"combined/" #final output directory (after concatening step 1 and 2 and reverse complementing)
Outdir_lost=$Outdir_final"filtered_out_results/" #reads lost during this step


if [ $Option == '5RACE' ]; then
	Outdir_polyA=$Outdir_final"polyU_reads/" #results for polyU tail longer than 8As that could potentially not map to the target sequence"
else
	Outdir_polyA=$Outdir_final"polyA_reads/" #results for polyU tail longer than 8As that could potentially not map to the target sequence"
fi

###Program and scripts:
removeduplicate=$Root"removeduplicate_Q10.py"
removeduplicate_minus=$Root"removeduplicate_Q20_minus.py"

targetresearch=$Root"target_research.py" #this script is set-up to analyze three target but can be easily modified for more or less targets
RNAdelimiter=$Root"RNA_delimiter.py"
RNAdelimiter_minus=$Root"RNA_delimiter_minus.py"
mappingloops=$Root"mapping_loops.py"
mappingloops_minus=$Root"mapping_loops_minus.py"
trimming=$Root"trimming.py"

extension_analysis_A=$Root"extension_analysis_blockA.py"
extension_analysis_B1=$Root"extension_analysis_blockB1.py"
extension_analysis_B2=$Root"extension_analysis_blockB2.py"

extension_analysis_B1_minus=$Root"extension_analysis_blockB1_minus.py"
extension_analysis_B2_minus=$Root"extension_analysis_blockB2_minus.py"

##########################################
# Create Output directories
##########################################
if [ ! -d $Outdir ]
then
	mkdir $Outdir
	echo $Outdir": doesn't exist. Created now"
else
	echo $Outdir": exists"
fi

if [ ! -d $Outdir_d ]
then
	mkdir $Outdir_d
	echo $Outdir_d": doesn't exist. Created now"
else
	echo $Outdir_d": exists"
fi

if [ ! -d $Outdir_dR1 ]
then
	mkdir $Outdir_dR1
	echo $Outdir_dR1": doesn't exist. Created now"
else
	echo $Outdir_dR1": exists"
fi

if [ ! -d $Outdir_dR2 ]
then
	mkdir $Outdir_dR2
	echo $Outdir_dR2": doesn't exist. Created now"
else
	echo $Outdir_dR2": exists"
fi
if [ ! -d $Outdir_targetR1 ]
then
	mkdir $Outdir_targetR1
	echo $Outdir_targetR1": doesn't exist. Created now"
else
	echo $Outdir_targetR1": exists"
fi
if [ ! -d $Outdir_targetR2 ]
then
	mkdir $Outdir_targetR2
	echo $Outdir_targetR2": doesn't exist. Created now"
else
	echo $Outdir_targetR2": exists"
fi
if [ ! -d $Outdir_trimmedR2 ]
then
	mkdir $Outdir_trimmedR2
	echo $Outdir_trimmedR2": doesn't exist. Created now"
else
	echo $Outdir_trimmedR2": exists"
fi

if [ $Option == 'A' ] || [ $Option == 'B1' ] || [ $Option == 'B1B2' ] || [ $Option == '5RACE' ]; then
	if [ ! -d $Outdir_ext ]
	then
		mkdir $Outdir_ext
		echo $Outdir_ext": doesn't exist. Created now"
	else
		echo $Outdir_ext": exists"
	fi
	if [ ! -d $Outdir_other ]
	then
		mkdir $Outdir_other
		echo $Outdir_other": doesn't exist. Created now"
	else
		echo $Outdir_other": exists"
	fi
	if [ ! -d $Outdir_results ]
	then
		mkdir $Outdir_results
		echo $Outdir_results": doesn't exist. Created now"
	else
		echo $Outdir_results": exists"
	fi
fi


if [ ! -d $Outdir_final ]
then
	mkdir $Outdir_final
	echo $Outdir_final": doesn't exist. Created now"
else
	echo $Outdir_final": exists"
fi

if [ ! -d $Outdir_lost ]
then
	mkdir $Outdir_lost
	echo $Outdir_lost": doesn't exist. Created now"
else
	echo $Outdir_lost": exists"
fi

if [ $Option == 'B1' ] || [ $Option == 'B1B2' ] || [ $Option == 'B2' ] || [ $Option == '5RACE' ]; then
	if [ ! -d $Outdir_polyA ]
	then
		mkdir $Outdir_polyA
		echo $Outdir_polyA": doesn't exist. Created now"
	else
		echo $Outdir_polyA": exists"
	fi
fi

if [ $Option == 'B1' ] || [ $Option == 'B1B2' ] || [ $Option == '5RACE' ]; then
	if [ ! -d $Outdir_mapped ]
	then
		mkdir $Outdir_mapped
		echo $Outdir_mapped": doesn't exist. Created now"
	else
		echo $Outdir_mapped": exists"
	fi
fi

if [ $Option == 'B1' ] || [ $Option == 'B1B2' ] || [ $Option == 'A' ] || [ $Option == '5RACE' ]; then
	if [ ! -d $Outdir_mapped_PA ]
	then
		mkdir $Outdir_mapped_PA
		echo $Outdir_mapped_PA": doesn't exist. Created now"
	else
		echo $Outdir_mapped_PA": exists"
	fi
fi

if [ $Option == 'B1B2' ] || [ $Option == '5RACE' ]; then
	if [ ! -d $Outdir_combined ]
	then
		mkdir $Outdir_combined
		echo $Outdir_combined": doesn't exist. Created now"
	else
		echo $Outdir_combined": exists"
	fi
fi

####################################################################################
##Pipeline
####################################################################################
for i in ${Datadir_R1}*.fastq.gz;do
	echo $i
	prefix=$(basename $i _L001_R1_001.fastq.gz) # prefix to be used for all output names
	for j in ${Datadir_R2}*.fastq.gz;do
		prefix_R2=$(basename $j _L001_R2_001.fastq.gz) #to be adjusted according to the name of your fastq file
		if [ $prefix_R2 == $prefix ]; then #analyse R1 file and its corresponding read 2 file
			echo "*******************RACEseq_Virusv3.1*******************"
			echo $(basename $i)
			echo $(basename $j)

			if [ $Option == '5RACE' ]; then
				#step 1: "removeduplicate.py" deduplication of sequences with identical nucleotides
				#in read 1 (insert) and the 1st to 8th cycle in read 2
				#(randomized bases in 3’ adapter)
				Outname_dR1=${Outdir_dR1}${prefix}.fastq.gz
				Outname_dR2=${Outdir_dR2}${prefix}.fastq.gz
				echo "_______Step1_deduplication_______"
				time $removeduplicate_minus $i $j $Outname_dR1 $Outname_dR2
			
			else
				#step 1: "removeduplicate.py" deduplication of sequences with identical nucleotides
				#in read 1 (insert) and the 1st to 15th cycle in read 2
				#(randomized bases in 3’ adapter)
				Outname_dR1=${Outdir_dR1}${prefix}.fastq.gz
				Outname_dR2=${Outdir_dR2}${prefix}.fastq.gz
				echo "_______Step1_deduplication_______"
				time $removeduplicate $i $j $Outname_dR1 $Outname_dR2
			fi

			#step 2: "target_research.py" extraction of read1 that match to target
			Infilename=$Outname_dR1
			Outname=${Outdir_targetR1}${prefix}.fastq.gz
			echo "_______Step2_target research_______"
			time $targetresearch $Infilename $Outname $Nb_mmR1 $list_target
			
			if [ $Option == '5RACE' ]; then
				#step 3: "RNA_delimiter.py" extracting read 2 with delimiter
				#and removing delimiter + randomized sequence
				Infilename_R1=$Outname
				Infilename_R2=$Outname_dR2
				Outname=${Outdir_targetR2}${prefix}.fastq.gz
				echo "_______Step3_delimiter research_______"
				time $RNAdelimiter_minus $Infilename_R1 $Infilename_R2 $Outname
			
			else
				#step 3: "RNA_delimiter.py" extracting read 2 with delimiter
				#and removing delimiter + randomized sequence
				Infilename_R1=$Outname
				Infilename_R2=$Outname_dR2
				Outname=${Outdir_targetR2}${prefix}.fastq.gz
				echo "_______Step3_delimiter research_______"
				time $RNAdelimiter $Infilename_R1 $Infilename_R2 $Outname
			fi
			
			#step 4: "trimming.py" removing 5' PCR primer sequence in read 2  
			Infilename=$Outname
			Outname=${Outdir_trimmedR2}/${prefix}.fastq.gz
			echo "_______Step4_trimming of the PCR primer sequence_______"
			time $trimming $Infilename $Outname 
			
			echo "_______Step5_Read2_analysis_______"
			echo "Selected option: ${Option}"
			
			if [ $Option == 'A' ]; then
				#BlockA.1: "mapping_loops.py" mapping of reads 2 to the RNA reference sequence
				Infilename=$Outname
				Outname=${Outdir_results}/${prefix}.txt
				Outname_other=${Outdir_other}/${prefix}.fastq.gz
				echo "_______BlockA.1_mapping of reads 2 to the reference sequence_______"
				time $mappingloops $Infilename $Outname $Outname_other $Reference_seq
			
				#BlockA.2: "extension_analysing.py" filtering 3' tails : 3’ modifications longer than 6, 10 and 15 nt 
				#were considered only if they contained at least one stretch of AAA or TTT, 
				#two stretches of AAA or TTT, or three stretches of AAA or TTT, respectively
				Infilename=$Outdir_results/${prefix}.txt
				Outname_mapped_PA=${Outdir_mapped_PA}/${prefix}.txt
				Outname_lost=${Outdir_lost}/${prefix}.txt
				echo "_______BlockA.2_filtering 3' tails and reverse complementing sequences_______"
				time $extension_analysis_A $Infilename $Outname_mapped_PA $Outname_lost 
			fi
			
			if [ $Option == 'B1' ]; then
				#BlockB1.1: "mapping_loops.py" mapping of reads 2 to the viral RNA reference sequence
				Infilename=$Outname
				Outname=${Outdir_results}${prefix}.txt
				Outname_other=${Outdir_other}${prefix}.fastq
				echo "_______BlockB1.1_mapping of reads 2 to the target reference sequences_______"
				time $mappingloops $Infilename $Outname $Outname_other $Reference_seq
			
				#BlockB1.2: "extension_analysing_loop.py" Analysis of 3' tail found after read2 mapping.
				#Research of potential other modification (U or C, G) on these tails
				#Reverse complementing read 2
				Infilename=$Outdir_results${prefix}.txt
				Infilename2=$Outname_other
				Outname_mapped=${Outdir_mapped}${prefix}.txt
				Outname_mapped_PA=${Outdir_mapped_PA}${prefix}.txt
				echo "___BlockB1.2_analysing 3' tails found after read 2 mapping___"
				echo "Number of mm accepted in the polyU tail ${Nb_mm}"
				time $extension_analysis_B1 $Infilename $Outname_mapped $Outname_mapped_PA $Infilename2 $Nb_mm
			fi
			if [ $Option == 'B2' ]; then
				#BlockB2: "extension_analysing_other.py" Research of potential long A tail (> 30nt)
				#in the read 2 that didn't  map target sequence.
				#Research of potential other modification (U or C, G) on these tails and reverse complementing read 2
				Infilename=$Outname
				Outname_polyU=${Outdir_polyA}${prefix}.txt
				Outname_lost=${Outdir_lost}${prefix}.txt
				echo "___BlockB2_looking for long polyU tail >= 8 in read 2___"
				echo "Number of mm accepted in the polyU tail ${Nb_mm}"
				time $extension_analysis_B2 $Infilename $Outname_polyU $Outname_lost $Nb_mm
			fi
			
			if [ $Option == 'B1B2' ]; then
				#BlockB1.1: "mapping_loops.py" mapping of reads 2 to the viral RNA reference sequence
				Infilename=$Outname
				Outname=${Outdir_results}${prefix}.txt
				Outname_other=${Outdir_other}${prefix}.fastq.gz
				echo "_______BlockB1.1_mapping of reads 2 to the target reference sequences_______"
				time $mappingloops $Infilename $Outname $Outname_other $Reference_seq
			
				#BlockB1.2: "extension_analysing_loop.py" Analysis of 3' tail found after read2 mapping.
				#Research of potential other modification (U or C, G) on these tails
				#Reverse complementing read 2
				Infilename=$Outdir_results${prefix}.txt
				Infilename2=$Outname_other
				Outname_mapped=${Outdir_mapped}${prefix}.txt
				Outname_mapped_PA=${Outdir_mapped_PA}${prefix}.txt
				echo "___BlockB1.2_analysing 3' tails found after read 2 mapping___"
				echo "Number of mm accepted in the polyU tail ${Nb_mm}"
				time $extension_analysis_B1 $Infilename $Outname_mapped $Outname_mapped_PA $Infilename2 $Nb_mm
				
				#BlockB2: "extension_analysing_other.py" Research of potential long A tail (> 30nt)
				#in the read 2 that didn't  map target sequence.
				#Research of potential other modification (U or C, G) on these tails and reverse complementing read 2
				Infilename=$Outname_other
				Outname_polyU=${Outdir_polyA}${prefix}.txt
				Outname_lost=${Outdir_lost}${prefix}.txt
				echo "___BlockB2_looking for long polyU tail >= 8 for read 2 that didn't map target sequence: ___"
				echo "Number of mm accepted in the polyU tail ${Nb_mm}"
				time $extension_analysis_B2 $Infilename $Outname_polyU $Outname_lost $Nb_mm
				
				#BlockB1B2: Concatener output files for step6 and step7 (take line header of file step 5)
				Outname_combined=${Outdir_combined}${prefix}.txt
				echo "___BlockB1B2_Combining results from blockB1 and blockB2___"
				time head -1 $Outname_mapped > $Outname_combined; tail -n +2 -q $Outname_mapped >> $Outname_combined; tail -n +2 -q $Outname_polyU >> $Outname_combined
			fi
			
			if [ $Option == '5RACE' ]; then
				#BlockB1.1: "mapping_loops.py" mapping of reads 2 to the viral RNA reference sequence
				Infilename=$Outname
				Outname=${Outdir_results}${prefix}.txt
				Outname_other=${Outdir_other}${prefix}.fastq.gz
				echo "_______BlockB1.1 for 5RACE_mapping of reads 2 to the target reference sequences_______"
				time $mappingloops_minus $Infilename $Outname $Outname_other $Reference_seq
			
				#BlockB1.2: "extension_analysing_loop.py" Analysis of 3' tail found after read2 mapping.
				#Research of potential other modification (U or C, G) on these tails
				#Reverse complementing read 2
				Infilename=$Outdir_results${prefix}.txt
				Infilename2=$Outname_other
				Outname_mapped=${Outdir_mapped}${prefix}.txt
				Outname_mapped_PA=${Outdir_mapped_PA}${prefix}.txt
				echo "___BlockB1.2 for 5RACE_analysing 3' tails found after read 2 mapping___"
				echo "Number of mm accepted in the polyU tail ${Nb_mm}"
				time $extension_analysis_B1_minus $Infilename $Outname_mapped $Outname_mapped_PA $Infilename2 $Nb_mm
				
				#BlockB2: "extension_analysing_other.py" Research of potential long A tail (> 30nt)
				#in the read 2 that didn't  map target sequence.
				#Research of potential other modification (U or C, G) on these tails and reverse complementing read 2
				Infilename=$Outname_other
				Outname_polyU=${Outdir_polyA}${prefix}.txt
				Outname_lost=${Outdir_lost}${prefix}.txt
				echo "___BlockB2 for 5RACE_looking for long polyU tail >= 8 for read 2 that didn't map target sequence: ___"
				echo "Number of mm accepted in the polyU tail ${Nb_mm}"
				time $extension_analysis_B2_minus $Infilename $Outname_polyU $Outname_lost $Nb_mm
				
				#BlockB1B2: Concatener output files for step6 and step7 (take line header of file step 5)
				Outname_combined=${Outdir_combined}${prefix}.txt
				echo "___BlockB1B2 for 5RACE_Combining results from blockB1 and blockB2___"
				time head -1 $Outname_mapped > $Outname_combined; tail -n +2 -q $Outname_mapped >> $Outname_combined; tail -n +2 -q $Outname_polyU >> $Outname_combined
			fi
		fi
	done
done
