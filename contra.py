#!/usr/bin/env python

# ----------------------------------------------------------------------#
# Copyright (c) 2011, Richard Lupat & Jason Li.
#
# > Source License <
# This file is part of CONTRA.
#
#    CONTRA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CONTRA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CONTRA.  If not, see <http://www.gnu.org/licenses/>.
#
# 
#-----------------------------------------------------------------------#
# Last Updated : 23 July 2012 16:43PM


import os
from optparse import OptionParser
import sys
import subprocess
import shlex
from multiprocessing import Process, Manager

from scripts.assign_bin_number_v2 import *
from scripts.average_count import *
from scripts.cn_apply_threshold import *
from scripts.convert_gene_coordinate import *
from scripts.convert_targeted_regions import *
from scripts.split_chromosome import *
from scripts.vcf_out import *
from scripts.get_chr_length import *
from scripts.count_libsize import *
from scripts.target_breakdown import *

#VERSION
VERSION="2.0.8"

#Absolute Path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

class Params:
	"""
	Class for top-level system parameters
	"""
		
	def __init__(self):
		# command-line option definition
		self.parser = OptionParser()
		self.parser.add_option("-t", "--target", 
			help="Target region definition file [REQUIRED] [BED Format]",
			action="store", type="string", dest="target")
		self.parser.add_option("-s", "--test", 
			help="Alignment file for the test sample [REQUIRED] [BAM/SAM]",
			action="store", type="string", dest="test")	
		self.parser.add_option("-c", "--control", 
			help="Alignment file for the control sample [REQUIRED] [BAM/SAM]",
			action="store", type="string", dest="control") 
		self.parser.add_option("-f", "--fasta", 
			help="Reference Genome [NOT REQUIRED since v2.0.8][FASTA]",
			action="store", type="string", dest="fasta")
		self.parser.add_option("-o", "--outFolder", 
			help="the output folder path name to store the output of analysis [REQUIRED]",
			action="store", type="string", dest="outFolder") 
		self.parser.add_option("--numBin", 
			help="Numbers of bins to group regions. User can specify multiple experiments with different number of bins (comma separated) [20]",
			action="store", type="string", dest="numBin", default="20") 
		self.parser.add_option("--minReadDepth", 
			help="The threshold for minimum read depth for each bases [10]",
			action="store", type="string", dest="minReadDepth", default=10) 
		self.parser.add_option("--minNBases",
			help="The threshold for minimum number of bases for each target regions [10]",
			action="store", type="string", dest="minNBases", default= 10) 
		self.parser.add_option("--sam", 
			help="If the specified, test and control sample are in SAM [False]",
			action="store_true", dest="sam", default="False") 
		self.parser.add_option("--bed",
			help="if specified, control will be in BED format [False]",
			action="store_true", dest = "bedInput", default="False")
		self.parser.add_option("--pval", 
			help="The p-value threshold for filtering [0.05]. Applies to Adjusted P-Value.",
			action="store", type="string", dest="pval", default=0.05) 
		self.parser.add_option("--sampleName", 
			help ="The name to be appended to the front of default output name ['']",
			action="store", type="string", dest="sampleName", default='')
		self.parser.add_option("--nomultimapped", 
			help="The option to remove multi-mapped reads [False]",
			action="store_true", dest="nomultimapped",default="False")	
		self.parser.add_option("-p", "--plot", 
			help="Plots log-ratio distribution for each bin [False]", 
			action="store_true", dest="plot", default="False")
		self.parser.add_option("--minExon",
			help="Minimum number of Exons in one bin (if less than this, bin that contains small number of exons"
				+"will be moved to the adjacent bins) [2000] ",
			action="store", type="string", dest="minExon", default="2000")
		self.parser.add_option("--minControlRdForCall",
			help="Minimum control readdepth for call [5]",
			action="store", type="string", dest="minControl", default="5")

		self.parser.add_option("--minTestRdForCall",
			help="Minimum test readdepth for call [0]",
			action="store", type="string", dest="minTest", default="0")

		self.parser.add_option("--minAvgForCall",
			help="Minimum average coverage for call [20]",
			action="store", type="string", dest="minAvg", default="20")

		self.parser.add_option("--maxRegionSize",
			help="Maximum Region Size in target region (for breaking large region into smaller region. By default, maxRegionSize 0 means no breakdown) [0]",
			action="store", type="string", dest="maxRegionSize", default="0") 

		self.parser.add_option("--targetRegionSize",
			help="Target Region Size for breakdown (if maxRegionSize is non zero) [200]",
			action="store", type="string", dest="targetRegionSize", default="200")

		self.parser.add_option("-l", "--largeDeletion",
                        help="if specified, CONTRA will run large deletion analysis (CBS). User must have DNAcopy R-library installed to run the analysis. [False]",
                        action="store_true", dest = "large", default="False")

		self.parser.add_option("--smallSegment",
			help="CBS segment size for calling large variations [1]",
			action="store", type="string", dest="smallSegment", default="1")

		self.parser.add_option("--largeSegment",
			help="CBS segment size for calling large variations [25]",
			action="store", type="string", dest="largeSegment", default="25")

		self.parser.add_option("--lrCallStart",
			help="Log ratios start range that will be used to call CNV [-0.3]",
			action="store", type="string", dest="lrs", default="-0.3")

		self.parser.add_option("--lrCallEnd",
			help="Log ratios end range that will be used to call CNV [0.3]",
			action="store", type="string", dest="lre", default="0.3")

		self.parser.add_option("--passSize", 
			help="Size of exons that passed the p-value threshold compare to the original exon size [0.35]",
			action="store", type="string", dest="passSize", default="0.35")

                ###
		self.parser.add_option("--removeDups",
			help="if specified, will remove PCR duplicates [False]",
			action="store_true", dest = "removeDups", default="False")	
		self.parser.add_option("--version",
			help="Returns version",
			action="store_true", dest = "version")

		# required parameters list
		self.ERRORLIST = []

		# change system parameters based on any command line
		(options, args) = self.parser.parse_args()
		if options.target:
			self.TARGET = options.target
		else:
			#self.parser.print_help()
			#self.parser.error("--target not supplied")
			self.ERRORLIST.append("target")

		if options.test:
			self.TEST = options.test
		else:
			#self.parser.error("--test not supplied")
			self.ERRORLIST.append("test")

		if options.control:
			self.CONTROL = options.control
		else:
			#self.parser.error("--control not supplied")
			self.ERRORLIST.append("control")

#		if options.fasta:
#			self.FASTA = options.fasta
#		else:
#			#self.parser.error("--fasta not supplied")
#			self.ERRORLIST.append("fasta")

		if options.outFolder:
			self.OUTFOLDER = options.outFolder
		else:
			#self.parser.error("--outFolder not supplied")
			self.ERRORLIST.append("outfolder")

		if len(self.ERRORLIST) != 0:
			self.parser.print_help()
			self.parser.error("Missing required parameters")	

		if options.numBin:
			binsNumber = options.numBin.split(",")
			try:
				self.NUMBIN = [int(j) for j in binsNumber]
			except:
				self.NUMBIN = [20]
		if options.minReadDepth:
			self.MINREADDEPTH = int(options.minReadDepth)
		if options.minNBases:
			self.MINNBASES = int(options.minNBases)
		if options.sam:
			self.SAM = str(options.sam)
		if options.pval:
			self.PVAL = options.pval
		if options.sampleName:
			self.SAMPLENAME = options.sampleName
		else:
			self.SAMPLENAME = 'No-SampleName'
		if options.nomultimapped:
			self.NOMULTIMAPPED = str(options.nomultimapped)
		if options.plot:
			self.PLOT = str(options.plot)
		if options.bedInput:
			self.BEDINPUT = options.bedInput
		if options.minExon:
			self.MINEXON	= int(options.minExon)
		if options.minControl:
			self.MINCONTROL = options.minControl
		if options.minTest:
			self.MINTEST = options.minTest
		if options.minAvg:
			self.MINAVG  = options.minAvg
		if options.maxRegionSize:
			self.MAXREGIONSIZE = int(options.maxRegionSize)
		if options.targetRegionSize:
			self.TARGETREGIONSIZE = int(options.targetRegionSize)
		if options.large:
			self.LARGE	= str(options.large)
		if options.smallSegment:
			self.SMALLSEGMENT = options.smallSegment
		if options.largeSegment:
			self.LARGESEGMENT = options.largeSegment
		if options.lre:
			self.LRE	= options.lre
		if options.lrs:
			self.LRS	= options.lrs
		if options.passSize:
			self.PASSSIZE	= options.passSize
			
                ### either "False" or True atn
		if options.removeDups:
			self.REMOVEDUPS = str(options.removeDups)
			
	def repeat(self):
		# params test
		print "target		:", self.TARGET
		print "test		:", self.TEST
		print "control		:", self.CONTROL
#		print "fasta		:", self.FASTA
		print "outfolder	:", self.OUTFOLDER
		print "numBin		:", self.NUMBIN
		print "minreaddepth	:", self.MINREADDEPTH
		print "minNBases	:", self.MINNBASES
		print "sam		:", self.SAM
		print "pval		:", self.PVAL
		print "sampleName	:", self.SAMPLENAME
		print "nomultimapped	:", self.NOMULTIMAPPED
		print "plot		:", self.PLOT
		print "bedInput		:", self.BEDINPUT
		print "minExon		:", self.MINEXON
		print "largeDeletion	:", self.LARGE
		print "removeDups	:", self.REMOVEDUPS
def checkOutputFolder(outF):
	print "Creating Output Folder :",
 
	if outF[len(outF)-1] == "/":
		outF = outF[:len(outF)-1]

	try:
		os.mkdir(outF)
	except:
		print "output folder already exists: " + outF
#		print "cannot create folder '%s'" %outF
#		print "if folder already exist, please specify other folder"
#		sys.exit(1)
	
	try:
		os.mkdir(outF+"/table")
		os.mkdir(outF+"/plot")
		os.mkdir(outF+"/buf")
		os.mkdir(outF+"/buf/ctrData/")
		os.mkdir(outF+"/buf/testData/")
	except:
		print "[ERROR: CANNOT CREATE SUBFOLDERS]"
		sys.exit(1)

	print " Done."

	return outF


#BEDINPUT
def countTotalReads3(params, folder):
	tempFileName    = folder + "/temp.txt"
        tempReadFile    = open(tempFileName, "w")
	libsize		= get_libsize(params.BEDINPUT)
	tempReadFile.write(libsize)
        #tempReadFile.write(params.CONTROLREADCOUNT)
        tempReadFile.close()


def countTotalReads(params, folder):
	if 'testData' in folder:
		inF = params.TEST
	else:
		inF = params.CONTROL

	# Get Total ReadCount
	getreadcount = os.system("samtools view %s | wc -l > %s/temp.txt" %(inF,folder))

def samToBam(samfile, bamfile):
	args = shlex.split("samtools view -bS %s -o %s" %(samfile, bamfile))
	samtobam = subprocess.call(args)

	return bamfile

def removeMultiMapped(inF, newBAM):
        # Get New BAM Files with mapping quality > 0
        args = shlex.split("samtools view -bq 1 %s -o %s" %(inF, newBAM))
	removeMM = subprocess.call(args)
        print "Multi mapped reads removed. "
        
def removeDups(inF, newBAM):
        # Remove
        args = shlex.split("samtools view -b -F 0x400 %s -o %s" %(inF, newBAM))
	removeDupsCall = subprocess.call(args)
        print "Removed PCR duplicates. "
        
#BEDINPUT
def convertBamSimple(params, folder, targetList, genomeFile):
	if 'testData' in folder:
                inF = params.TEST
                print "Converting TEST Sample... "
        else:
                inF = params.CONTROL
                print "Converting CONTROL Sample... "
	
	#Copy file to working folder
	os.system("cp %s %s" %(inF, folder+"sample.BEDGRAPH"))

	# Split Bedgraph by its chromosomes
        splitByChromosome(folder)

        # Slice the coverage files to only cover the targeted regions
        print "Getting targeted regions DOC..."
        convertGeneCoordinate(targetList, folder)
        
	# LIBSIZE
	libsize = str(get_libsize(folder+"geneRefCoverage2.txt"))
	tempLibSize = open(folder + "/temp.txt", "w")
	tempLibSize.write(libsize)
	tempLibSize.close()

	print "Targeted regions pre-processing: Done"


def convertBam(params, folder, targetList, genomeFile):
	if 'testData' in folder:
		inF = params.TEST
		print "Converting TEST Sample... "
	else:
		inF = params.CONTROL
		print "Converting CONTROL Sample... "
	
	# Convert BAM Files to BEDGRAPH
	bedgraph = folder + "sample.BEDGRAPH"
	args = shlex.split("genomeCoverageBed -ibam %s -bga -g %s" %(inF, genomeFile))
	print "DEBUG 123 " + " ".join(args)
	#output = subprocess.Popen(args, stdout = subprocess.PIPE).communicate()[0]
	iOutFile = open(bedgraph, "w")
	#iOutFile.write(output)
	output	= subprocess.Popen(args, stdout = iOutFile).wait()
	iOutFile.close()

	# Split Bedgraph by its chromosomes
	splitByChromosome(folder)

	# Slice the coverage files to only cover the targeted regions
	print "Getting targeted regions DOC..."
	convertGeneCoordinate(targetList, folder)
	
	# LIBSIZE
        libsize = str(get_libsize(folder+"geneRefCoverage2.txt"))
	tempLibSize = open(folder + "temp.txt", "w")
	tempLibSize.write(libsize)
	tempLibSize.close()

	print "Targeted regions pre-processing: Done"		

def analysisPerBin(params, num_bin, outFolder, targetList):
	import shutil
	
	bufLoc = outFolder + "/buf"
	# Assign bin number to the median and average file
	numBin = assignBin(num_bin, bufLoc+"/average.txt", bufLoc+"/bin", targetList, params.MINEXON)
	
	#copy bin_boundary to plot folder
	#outBounFile = os.path.join(outFolder,  "plot", "bin_boundary"+str(num_bin))
	#curBounFile = os.path.join(bufLoc, "bin" + str(num_bin) + ".boundaries.txt")
	#shutil.copy(curBounFile, outBounFile)


	print "Significance Test ...  "
	rScriptName = os.path.join(scriptPath, "scripts", "cn_analysis.v3.R")
	args = shlex.split("Rscript %s %s %s %s %s %s %s %s %s %s %s" 
		%(rScriptName, num_bin, params.MINREADDEPTH, params.MINNBASES, outFolder, params.SAMPLENAME,params.PLOT, numBin, params.MINCONTROL, params.MINTEST, params.MINAVG))
	rscr = subprocess.call(args)


	print "Generating Output Files ... "
	# Analysis of CNV
	tNameList = os.listdir(outFolder+"/table/")
	if num_bin > 1:
		tNameId = str(num_bin) + "bins"
	else:
		tNameId = str(num_bin) + "bin"
        for tName in tNameList:
		if tNameId in tName:
			break

        if "CNATable" in tName:
		tName = tName[:len(tName)-4]
		tableName = outFolder + "/table/" + tName
		bufTable  = bufLoc + "/" + tName
		applyThreshold(tableName, bufTable, params.PVAL, 100000) #params.MAXGAP = 100000

		# Large Region CBS
		if (params.LARGE != "False"):
			print "DEBUG 266a"
			rScriptName2 = os.path.join(scriptPath, "scripts", "large_region_cbs.R")
			args = shlex.split("Rscript %s %s %s %s %s %s %s %s %s"
			%(rScriptName2, tableName+".txt", params.SMALLSEGMENT, params.LARGESEGMENT, params.PVAL, params.PASSSIZE, params.LRS, params.LRE, bufLoc))
			rscr2 = subprocess.call(args)
			print str(args)
		else:
			print "DEBUG 266b"

		# Generate the DNA sequence (for VCF file)
		print  "Skipping VCF generation.. use tabular file instead."
#		bedFile  = bufTable + ".BED"
#		bedFasta = bufTable + ".fastaOut.txt"
#		fastaFile = params.FASTA
#		args = shlex.split("fastaFromBed -fi %s -bed %s -fo %s -name"
#				%(fastaFile, bedFile, bedFasta))
#		print (" ".join(args))
#		fastaBED = subprocess.call(args)
#
#		# Write VCF
#		print  "Creating VCF file ... "
#		vcfFile = tableName + ".vcf"
#		vcf_out(bedFasta, vcfFile)
#
#		print "%s created. " %(vcfFile)

        else:
		print "Table not found"

def removeTempFolder(tempFolderPath):
	import shutil
	
	shutil.rmtree(tempFolderPath)

	print "Temp Folder Removed"


def main():
	if len(sys.argv) == 2 and sys.argv[1] == "--version":
		print VERSION
		sys.exit(1)
	
	# option handling
	params = Params()
	params.repeat()

	# output folder handling
	outFolder = checkOutputFolder(params.OUTFOLDER)
	bufLoc = outFolder + "/buf"

	# convert target file
	sorted_target = os.path.join(bufLoc, "target.BED")
	os.system("sort -k1,1 -k2n %s > %s" %(params.TARGET, sorted_target))	

	# target breakdown
	if params.MAXREGIONSIZE > 0:
		new_target = os.path.join(bufLoc, "target_breakdown.BED")
		target_breakdown(sorted_target, params.MAXREGIONSIZE, params.TARGETREGIONSIZE, new_target)
		sorted_target = new_target

	targetList = convertTarget(sorted_target)

	# convert sam to bam if -sam specified
	if (params.SAM == "True"):
		print "Pre-processing SAM files"

		test_bam = bufLoc + "/test.BAM"
		ctr_bam  = bufLoc + "/control.BAM"

		samTest = Process(target= samToBam, args=(params.TEST, test_bam))
		if params.BEDINPUT == "False":
			samCtr = Process(target= samToBam, args=(params.CONTROL, ctr_bam))

		samTest.start()
		if params.BEDINPUT == "False":
			samCtr.start()

		samTest.join()
		if params.BEDINPUT == "False":
			samCtr.join()

		params.TEST = test_bam
		if params.BEDINPUT == "False":
			params.CONTROL = ctr_bam

	# remove multi mapped reads if --nomultimapped is specified
	if (params.NOMULTIMAPPED == "True"):
		print "Removing multi-mapped reads"

		test_bam = bufLoc + "/test_reliable.BAM"
                ctr_bam  = bufLoc + "/control_reliable.BAM"

                bamTest = Process(target= removeMultiMapped, args=(params.TEST, test_bam))
		if params.BEDINPUT == "False":
                	bamCtr = Process(target= removeMultiMapped, args=(params.CONTROL, ctr_bam))

                bamTest.start()
		if params.BEDINPUT == "False":
	                bamCtr.start()

                bamTest.join()
		if params.BEDINPUT == "False":
	                bamCtr.join()

                params.TEST = test_bam
		if params.BEDINPUT == "False":
	                params.CONTROL = ctr_bam

        ###
	# Remove PCR duplicates if --removeDups specified
	if (params.REMOVEDUPS == "True"):
                print "Removing reads marked as duplicates (PCR)"
                
		test_bam = bufLoc + "/test_removedups.BAM"
                ctr_bam  = bufLoc + "/control_removedups.BAM"
                
                bamTest = Process(target = removeDups, args=(params.TEST, test_bam))
		if params.BEDINPUT == "False":
                	bamCtr = Process(target = removeDups, args=(params.CONTROL, ctr_bam))

                bamTest.start()
		if params.BEDINPUT == "False":
	                bamCtr.start()

                bamTest.join()
		if params.BEDINPUT == "False":
	                bamCtr.join()

                params.TEST = test_bam
		if params.BEDINPUT == "False":
	                params.CONTROL = ctr_bam
                
	
	# Get Chromosome Length
	genomeFile = bufLoc + '/sample.Genome'
	get_genome(params.TEST, genomeFile)

	# spawn bam converting scripts
	pTest = Process(target= convertBam, 
			args=(params, bufLoc+'/testData/', targetList, genomeFile))

	#BEDINPUT
	if params.BEDINPUT == "False":

		cTest = Process(target= convertBam, 
			args=(params, bufLoc+'/ctrData/' , targetList, genomeFile))
	else:
		cTest = Process(target= convertBamSimple,
			args=(params, bufLoc+'/ctrData/', targetList, genomeFile))
	# start the processes
	pTest.start()
	cTest.start()

	# wait for all the processes to finish before continuing
	pTest.join()
	cTest.join()

	# Get the read depth count from temporary folder
	for folder in [bufLoc+'/testData/', bufLoc+'/ctrData/']:	
		if 'testData' in folder:
			t1 = int(file.readlines(open(folder+"temp.txt"))[0].strip("\n"))
		else:
			n1 = int(file.readlines(open(folder+"temp.txt"))[0].strip("\n"))
	print "Test file read depth 	= ", t1
	print "Control file read depth 	= ", n1
	print "Pre-processing Completed. "
	
	# Get the Average of the Log Ratio	
	print "Getting the Log Ratio ... "
	testListName = bufLoc + '/testData/geneRefCoverage.txt'
	controlListName = bufLoc + '/ctrData/geneRefCoverage.txt'
	avOut = bufLoc + "/average.txt"
	averageCount(testListName, controlListName, avOut, t1, n1, params.MINREADDEPTH, params.MINNBASES)	

	# Analysis. [Bin, significance test, large deletion, vcf output] 	
	print "Binning ... "
	binProc = []
	for numBin in params.NUMBIN:	
		binProc.append(Process(target= analysisPerBin, args=(params,numBin,outFolder,targetList)))

	for proc in binProc:
		proc.start()

	for proc in binProc:
		proc.join()
		
	# Removed Temp Folder 
	removeTempFolder(bufLoc)
	
if __name__ == "__main__":
	main()
	print "Done... "
