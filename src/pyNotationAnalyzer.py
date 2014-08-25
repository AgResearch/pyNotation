#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
 pyNotation is a Functional Analysis software.
    Copyright (C) 2014 Clement DELESTRE (cclementddel@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import argparse
import sys
from Bio import SeqIO,pairwise2
from Bio import SeqFeature
import re
from Bio.SubsMat import MatrixInfo as matlist
import datetime
from ExternalToolsManager import ExternalToolsManager
from Visualizer import VisuManager

#########################
### Global variables ###
#########################
date = datetime.datetime.now().strftime('%Y-%m-%d')
author="Clement DELESTRE"
contact="cclementddel@gmail.com"
appliName="Py-notation Analyzer"
year=2014
version=1.0
#To count rRNA operons
operonThreshold=3000
# Alignments Parameters
matrix = matlist.blosum80 
openGap=-5
extendGap=-0.5
# threshold = 70%*ref length
thresholdMax=70 

#########################
######### Class #########
#########################

#File containing similar products
class OtherFile (object):
	def __init__(self,name):
		self.name=name
		self.otherProductsFound={}
	def addSequence(self,seqName):
		self.otherProductsFound[seqName]=[]
	def addProduct(self,prod,seqName,similarity=0):
		s=SimilarProducts(prod,similarity)
		try :
			self.otherProductsFound[seqName].append(s) 
		except KeyError:
			self.otherProductsFound[seqName]=[]
			self.otherProductsFound[seqName].append(s)
	def getNumberOfProducts(self):
		toReturn=0
		for keys in self.otherProductsFound:
			toReturn+=len(self.otherProductsFound[keys])
		return toReturn
	def __str__(self) :
		toReturn=""
		for key in self.otherProductsFound:
			for prod in self.otherProductsFound[key]:
				toReturn+="***LOCUS :\n"+str(key)+"\n***PRODUCT :\n"+str(prod)
		return toReturn

#Product similar to the target		
class SimilarProducts(object):
	def __init__(self,product,similarity=0):
		self.product=product # Suppose to be SeqFeature.SeqFeature
		self.similarity=similarity #need to be an integer
	def setSimilarity(similarity):
		self.similarity=similarity
	def __str__(self) :
		toReturn=""
		toReturn+=str(self.product)+"\n* Similarity : "+str(self.similarity)+" *\n"
		return toReturn

		
#Product we are looking for (the target)
class Product (object): 
	etm=ExternalToolsManager(None)
	def __init__(self, product, locus,firstFile):
		self.product=product # Suppose to be SeqFeature.SeqFeature
		self.locus=locus # Suppose to be SeqFeature
		self.otherFiles=[] # Suppose to be an OtherFile array
		self.sourceFile=firstFile # Suppose to be string
		self.name=self.product.qualifiers['product'][0].replace(" ","_") # Suppose to be string
	def __str__(self) :
		toReturn="** ORIGINAL LOCUS AND PRODUCT :\n"
		toReturn+=str(self.locus)+"\n"+str(self.product)+"\n"
		for file in self.otherFiles :
			nb = file.getNumberOfProducts()
			if nb==0:
				toReturn+="*** No similar products found in file "+str(file.name)+"\n"
			else :
				toReturn+="*** "+str(nb)+" similar product(s) found in file : "+str(file.name)+"\n"+str(file)
		return toReturn
		
	def getNumberOfSimilarProducts(self):
		toReturn=0
		for file in self.otherFiles :
			nb = file.getNumberOfProducts()
			toReturn+=nb
		return toReturn
	def getSummary(self):
		toReturn=""
		found=False
		for file in self.otherFiles :
			if file.getNumberOfProducts() != 0 :
				found=True
				for key in file.otherProductsFound :
					for prod in file.otherProductsFound[key]:
						toReturn+="\n"+self.product.qualifiers['product'][0].replace(",","_")+","+self.locus.name.replace(",","_")+","+prod.product.qualifiers['product'][0].replace(",","_")+","+key.name.replace(",","_")+","+str(prod.similarity)+" %,"+file.name
						
		if not found :
			toReturn+="\n"+self.product.qualifiers['product'][0].replace(",","_")+","+self.locus.name.replace(",","_")+",No similar products found"
		return toReturn
		
	def addFile(self,file):
		self.otherFiles.append(OtherFile(file))
		
	def lookingForFile(self,file):
		for other in self.otherFiles :
			if other.name==file:
				return other
		print"*** WARNING : FILE "+	file.name+" not found. None is return ***"	
		return None
	
	def addSequence(self,seqName,fileName):
		temp=self.lookingForFile(fileName)
		if temp :
			temp.addSequence(seqName)
		else :
			print "File "+fileName+" not found\nSystem will exit."
			sys.exit(1) 
			
	def addProduct(self,prod,seqName,fileName):
		temp=self.lookingForFile(fileName)
		if temp :
			self.etm.computeNeedle(prod,self.product)
			temp.addProduct(prod,seqName,self.etm.needleSimilarityPercent)
		else :
			print "File "+fileName+" not found\nSystem will exit."
			sys.exit(1) 
			
	def getSize(self):
		return self.product.location.end-self.product.location.start



#########################
###### Functions ########
#########################
def createOutputFile(suffix,ext): 
	outputFile = open(suffix.replace(" ","_")+ext.replace(" ","_"),'w')
	return outputFile

def countType(gb_record,type):
	i=0
	for feat in gb_record.features:
		if feat.type == type : 
			i=i+1
	return i

def preComputeFile(file) :
	record=[]
	for gb_record in SeqIO.parse(open(file,"r"), "genbank") :
		record.append(gb_record)	
	return record

def countCDSQualifiers(list):
	result={'db_xref':0,'EC_number':0,'translation':0,'locus_tag':0,'inference':0,'note':0,'protein_id':0,'Nb of CDS':0}
	for seq in list:
		for feat in seq.features :
			if feat.type=="CDS":
				result['Nb of CDS']+=1
				for key in result :
					try :
						feat.qualifiers[key]
						result[key]=result[key]+1
					except KeyError:
						pass
	return result

def researchProduct(record,prod,type="CDS"):
	hashMap={}
	pattern=re.compile(prod,re.IGNORECASE)
	for gb_record in record :
		CDSprod=[]
		for feat in gb_record.features:
			if feat.type!=type:
				continue
			match = re.search(pattern, feat.qualifiers['product'][0])	
			if match:
				CDSprod.append(feat)
		if len(CDSprod)!=0 :
			hashMap[gb_record]=CDSprod
	return hashMap

def CDSfound(dict) :
	for key in dict:
		print "#- In the following locus :"
		print key
		print "\n--*--\n"
		for feat in dict[key] :
			print feat
			
def nbOfProduct(dict):
	total=0
	for key in dict:
			total+=len(dict[key])
	return total
	
	
def computeAlignment(seqOne,seqTwo,feature,k=3.5): 
	threshold=len(seqOne)*k
	aln=pairwise2.align.globalds(seqOne, seqTwo, matrix,openGap, extendGap)
	if aln[0][2] < threshold :
		return False
	return True

def compareSeq(prod,feature):
	try :
		seqOne= prod.qualifiers['translation'][0]
		seqTwo= feature.qualifiers['translation'][0]
		return computeAlignment(seqOne,seqTwo,feature)
	except KeyError as ke :
		return False	

def computeLength(sizeRef,featuresTwo):
	sizeTwo=featuresTwo.location.end-featuresTwo.location.start
	if sizeRef > sizeTwo :
		diffSize=(float(sizeTwo)/sizeRef*100)
	else :
		diffSize=(float(sizeRef)/sizeTwo*100)
	return diffSize
		
def compareLength(feature,sizeRef,threshold): 
	d = computeLength(sizeRef,feature)
	if d < threshold: 
		return False
	return True
	
def compareECnumber(prod,feature):
	result=[]
	try : # Some feature could don't have EC number
		#Some feature can have severals EC number
		for ec in feature.qualifiers['EC_number']:
			for ecprod in prod.qualifiers['EC_number']:
				if ec == ecprod :
					result.append(feature)
	except KeyError:
		result.append(feature)
	return result
def countGenes(recordOne):
		total=0
		for seq in recordOne:
			total=total+countType(seq,"CDS")
		return total
def countCodingRegion(recordOne):
	total=0
	for seq in recordOne:
		start=end=0
		for feat in seq.features:
			if feat.type != "CDS" :
				continue
			start=feat.location.start
			if end >= start :
				start=end
			else :
				start=feat.location.start
			end=feat.location.end
			diff=end-start
			total+=diff
	return total
	
def countDNA(recordOne):
	total=0
	for seq in recordOne:
		total+=len(seq)
	return total
	
def countGC(recordOne):
	total=0
	for seq in recordOne:
		nbC=seq.seq.count('C')
		nbG=seq.seq.count('G')
		total=total+nbC+nbG
	return total
		
	
def countOneProduct(recordOne,prod):
	result = researchProduct(recordOne,prod,type="CDS")
	total=0
	for key in result:
		total+=len(result[key])
	result={}
	return total
	
def countHypProt(recordOne):
	total=0
	pattern=re.compile("hypothetical protein",re.IGNORECASE)
	for gb_record in recordOne :
		for feat in gb_record.features:
			if feat.type!="CDS":
				continue
			match = re.search(pattern, feat.qualifiers['product'][0])
			if match:
				total+=1
				
	return total
	
def countRrna(recordOne):
	nbRrna=0
	featTemp=None
	seqTemp=None
	for seq in recordOne :
		for feat in seq.features:
			if feat.type == "rRNA" :
				if featTemp:
					d=feat.location.start-featTemp.location.end
				if ((not featTemp) or ((not seqTemp) or seq!=seqTemp) or (feat.location.start-featTemp.location.end > operonThreshold )):
					nbRrna+=1
				seqTemp=seq
				featTemp=feat
	return nbRrna
	
def countNbSeq(recordOne):
	total=0
	for seq in recordOne:	
		total+=1
	return total
def analyseSeq (recordOne):
	nbNT=countDNA(recordOne)
	nbGC=countGC(recordOne)
	nbCodingNt=countCodingRegion(recordOne)
	result={"Genome size (bp)":nbNT,"DNA G+C content (bp)":nbGC,"DNA coding region (bp)":nbCodingNt}
	return result
	
def countRNAgenes(recordOne):
	nbRna=0
	featTemp=None
	seqTemp=None
	for seq in recordOne :
		for feat in seq.features:
			if feat.type == "rRNA" or feat.type == "tRNA":
				nbRna+=1
	return nbRna
	
def analyseGenes (recordOne,file):
	proteinCodingGene=countGenes(recordOne)
	nbRNA=countRNAgenes(recordOne)
	nbrRNA=countRrna(recordOne)
	nbHyp=countHypProt(recordOne)
	etm=ExternalToolsManager(file)
	etm.computeExternalToolsAnalysis()
	totalGene=nbRNA+proteinCodingGene
	functionGene=totalGene-nbHyp
	result={"Total gene":totalGene,"RNA gene":nbRNA,"rRNA operons":nbrRNA,"Genes with function prediction":functionGene,"Protein-coding genes":proteinCodingGene,"Genes with signal peptides":etm.nbOfSignal,"Genes assigned Pfam domains":etm.nbOfPfamDomains}
	return result
	
	
def getResults(dict,percentReference,sep=",",string="\""):
	toReturn=string+"Attribute"+string+sep+string+"Value"+string+sep+string+"% of Total"+string+"\n"
	for key in dict:
		toReturn+=string+key+string+sep+str(dict[key])+sep+str(float(dict[key])/dict[percentReference]*100)+"\n"
	return toReturn
	
#########################
### Parsing Arguments ###
#########################
parser = argparse.ArgumentParser(description="A program that analyze GBB files.",epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('file', metavar='File', nargs='+',help='file(s) you want analyzing (only in GenBank format)')
parser.add_argument("-s",metavar="product",help="looking for a product in several files (product must be in the first file, case insensitive)")

parser.add_argument("-a", "--analyze",help="analyze one or several file(s)",action="store_true")

parser.add_argument("-o",metavar="prefix output name",help="print only matching CDS in a file. Output file format is : prefix_product for -s option and prefix_analysis_file for -a option.")
parser.add_argument("-q", "--quiet",help="don't display meta data on the stdout.",action="store_true")

parser.add_argument("-x","--ext",metavar="outputfile name",help="extract all protein sequences in one fasta file.")

parser.add_argument("--pdf",metavar="prefix pdf file name",help="create a pdf containing the genetic neighborhood. Output file format is prefix_product.pdf.")
parser.add_argument("--interval",metavar="int",help="the interval of the genetic neighborhood (used with --pdf option) default : 5 .",default=5,type=int)
parser.add_argument("--summary",metavar="prefix output file name",help="print a summary of research (-s option) in a CSV file. Output file format is prefix_product.csv.")



args=parser.parse_args()
recordOne=[]
#########################
####### Compute #######
#########################
if not args.quiet :
	print "## "+appliName+" "+date+" ##"	
#Compute COG
if args.ext :
	for file in args.file :
		etm=ExternalToolsManager(file)
		etm.extAllProtSeqInOneFile(args.ext)
	if not args.quiet :
		print "## Extraction done at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
		
#Analyzing
if args.analyze :
	#COMPUTE ANALYSIS
	if not args.quiet :
		print "## Analysis launch at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
	for file in args.file :
		if not args.quiet :
			print "# Analyzing file "+file+" #"
		del recordOne[:]
		recordOne=preComputeFile(file)
		#Analyse GBK tags
		nbTag=countCDSQualifiers(recordOne)
		#Analyze Seq
		hashSeq=analyseSeq(recordOne)
		#Analyze genes
		hashGenes=analyseGenes(recordOne,file)
	#PRINT RESULTS
		if args.o :
			DNAseqFile=createOutputFile(args.o,"_DNAseq.csv")
			DNAseqFile.write(getResults(hashSeq,"Genome size (bp)"))
			DNAseqFile.close()
			if not args.quiet :
				print "## Results saved in file "+DNAseqFile.name+" at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
			CDStagsFile=createOutputFile(args.o,"_GBKtags.csv")
			CDStagsFile.write(getResults(nbTag,"Nb of CDS"))
			CDStagsFile.close()
			if not args.quiet :
				print "## Results saved in file "+CDStagsFile.name+" at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
			genesFiles=createOutputFile(args.o,"_genesInformations.csv")
			genesFiles.write(getResults(hashGenes,"Total gene"))
			genesFiles.close()
			if not args.quiet :
				print "## Results saved in file "+genesFiles.name+" at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
			
		else :
			print "# DNA sequence informations found in file : "+file+" #"
			print getResults(hashSeq,"Genome size (bp)",sep="\t",string="")
			print "# Genetic informations found in file : "+file+" #"
			print getResults(hashGenes,"Total gene",sep="\t",string="")
			print "# GenBanK Tag informations found in file : "+file+" #"
			print getResults(nbTag,"Nb of CDS",sep="\t",string="")
		
	if not args.quiet :
		print "## Analysis finish at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
#Searching
if args.s:
	allProds=[]
	#Looking for errors
	if len(args.file)==1:
		print "\n*** ERROR : -s option should use with several file. If you want to search this product in this first file please write it twice. ***\n"
		sys.exit(1)
	del recordOne[:]
	recordOne=preComputeFile(args.file[0])
	tableOfProducts=researchProduct(recordOne,args.s)
	if len(tableOfProducts)==0:
		print "\n*** ERROR : product "+args.s+" not found in file "+args.file[0]+" ***\n"
		sys.exit(2)
	#Compute searching
	if not args.quiet :
		print "## Searching launch at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
		print "# "+str(nbOfProduct(tableOfProducts))+" product(s) found in the first file : "+args.file[0]+" #"
		CDSfound(tableOfProducts)
	firstFile=args.file.pop(0)
	for key in tableOfProducts : #Original locus
		for prod in tableOfProducts[key] : #original product(s)	
			currentProd=Product(prod,key,firstFile)
			prodSize=currentProd.getSize()
			allProds.append(currentProd)
			if not args.quiet :
					print "#- Computing the following product : "
					print prod
			for otherFile in args.file :
				recordTwo=[]
				recordTwo=preComputeFile(otherFile)
				currentProd.addFile(otherFile)
				if not args.quiet :
					print "#- Computing the following file : "
					print otherFile
				for records in recordTwo : # Sequences
					currentProd.addSequence(records,otherFile)
					similarProducts=[]
					for otherProduct in records.features :
						if otherProduct.type=="CDS":
							temp=compareECnumber(prod,otherProduct)
							if len(temp)!=0:
								similarProducts.extend(temp)
					# http://stackoverflow.com/questions/1207406/remove-items-from-a-list-while-iterating-in-python/1207500	
					similarProducts[:] = [product for product in similarProducts if compareLength(product,prodSize,thresholdMax)]
					similarProducts[:] = [product for product in similarProducts if compareSeq(prod,product)]
					if len(similarProducts)!=0 :
						for productFound in similarProducts :
							currentProd.addProduct(productFound,records,otherFile)
							
	#Finally, print results	
	if args.summary :
		summaryFile=createOutputFile(args.summary,".csv")
		summaryFile.write("Original Product,Ref Locus,Similar Product,Query Locus,% similarity,File name")
	if args.o:
		outputFile=createOutputFile(args.o,"")
	for currentProd in allProds :
		if not args.o:
			print("**** Results from input file : "+firstFile+"\n")
			print(currentProd)
		else :
			outputFile.write("**** Results from input file : "+firstFile+"\n\n")
			outputFile.write(str(currentProd))
		if args.summary :
			summaryFile.write(currentProd.getSummary())
	#Save outputfile
	if args.o:
		outputFile.close()
		if not args.quiet :
			print "## Results saved in file "+outputFile.name+" at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
	#Save Summary File		
	if args.summary :			
		summaryFile.close()
		if not args.quiet :
			print "## Summary saved in file "+summaryFile.name+" at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
	if not args.quiet :
		print "## Search finish at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"
		

#Draw the PDF			
if args.pdf :
	allNames=[]
	vm = VisuManager(args.interval,args.pdf)
	for currentProd in allProds :
		if currentProd.getNumberOfSimilarProducts()==0:
			if not args.quiet :
				print "## No product found for "+currentProd.name+" so no PDF will be draw ##"
		else :
			vm.outputfileName=args.pdf+"_"+currentProd.name+".pdf"
			i=2
			while vm.outputfileName in allNames :
				vm.outputfileName=args.pdf+"_"+currentProd.name+"_"+str(i)+".pdf"
				i+=1
			allNames.append(vm.outputfileName)
			vm.addRefCDS(currentProd.locus,currentProd.product)
			for file in currentProd.otherFiles :
				for key in file.otherProductsFound :
					for otherProd in file.otherProductsFound[key] :
						b=vm.addQueryCDS(key,otherProd.product)
						#print b #debug : b should be True if everything OK
			vm.draw()
			vm.clean()
				
			if not args.quiet :
				print "## PDF saved in file "+vm.outputfileName+" at "+datetime.datetime.now().strftime('%H:%M:%S')+" ##"