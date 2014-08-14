#!/usr/bin/python
'''
 pyNotation is a Functional Analysis software.
    Copyright (C) 2014 Clément DELESTRE (cclementddel@gmail.com)

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
import re
import os 
from Bio import SeqIO,pairwise2
from Bio import Application
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Emboss.Applications import NeedleCommandline
#########################
######### Class #########
#########################

#TODO : 
# Check exit statut


class ExternalToolsManager (object):

	def __init__(self,GBKfile,remove=True,signalP=True,Hmmer=True,cog=True):
		self.GBKfile=GBKfile
		self.remove=remove
		self.signalP=signalP
		self.Hmmer=Hmmer
		self.cog=cog
		
		#Fasta File of temporary seq extracted
		self.fastaFile=""
		#Do we remove the fasta file after used ?
		self.remove=True
		#TemporaryFile
		self.tempFile="tmp"
		
		#Needle parameters
		self.needleIdentity=0 
		self.needleIdentityPercent=0
		self.needleSimilarity=0
		self.needleSimilarityPercent=0
		self.needleOutputFile="needle.txt"
		
		#HMMER parameters
		self.HMMscanbinCommand="hmmscan"
		self.HMMdb="Pfam-A.hmm"
		self.nbOfPfamDomains=0
		self.hmmScanOutput="hmmscan.output"
		self.indHmmerOutputFile=3
		
		
		#Signal P parameters 
		# bin command
		self.SignalPbinCommand="signalp"
		# Indice of information to parse in output files
		self.isAsignalPeptidInd=9
		#Number of signal peptides found
		self.nbOfSignal=0
		self.organismType="gram-"
		
	def setNewGBKFile(self,GBKfile) :
		self.GBKfile=GBKfile
	def flip(self):
		self.nbOfSignal=self.nbOfPfamDomains=0
		
	#Extract all sequence from a Locus	
	def extAllProtFromOneLocus(self,gb_record):
		fasTab=[]
		#Let s load proteins FASTA sequences
		for feat in gb_record.features: 
			if feat.type=="CDS":
				record=SeqRecord(Seq(feat.qualifiers['translation'][0],IUPAC.protein),id=gb_record.name.replace("|","_")+"_"+str(feat.location.start)+"_"+str(feat.location.end),description="")
				fasTab.append(record)
		#Let s write them
		if len(fasTab)!=0:
			self.fastaFile=gb_record.name.replace("|","_")+"_prot.fasta"
			file=open(self.fastaFile,"w")
			for seq in fasTab :
				file.write(seq.format("fasta"))
			file.close()
		else :
			self.fastaFile=None
			
	def extAllProtSeqInOneFile(self,outputName):
		tempname=""
		for gb_record in SeqIO.parse(open(self.GBKfile,"r"), "genbank") :
			self.extAllProtFromOneLocus(gb_record)
			if (self.fastaFile) and tempname != self.fastaFile :
				tempname=self.fastaFile
				os.popen("cat "+self.fastaFile+" >> "+outputName)
				os.remove(self.fastaFile)
	
	#Extract one sequence from a CDS
	def extProtSeqFromCDS(self,CDS,suffixe="_protSeq.fa"):
		nameWithoutSpace=CDS.qualifiers['product'][0].replace(" ","_")
		fileName=nameWithoutSpace.replace(";","_")
		fileName=fileName.replace("(","_")
		fileName=fileName.replace(")","_")
		fileName=fileName.replace("/","_")+suffixe
		record=SeqRecord(Seq(CDS.qualifiers['translation'][0],IUPAC.protein),id=nameWithoutSpace.replace(";","_"),description="")
		file=open(fileName,"w")
		file.write(record.format("fasta"))
		file.close()
		return fileName

		
			
	#Main method		
	def computeExternalToolsAnalysis(self):
		for gb_record in SeqIO.parse(open(self.GBKfile,"r"), "genbank") :
			self.extAllProtFromOneLocus(gb_record)
			if self.fastaFile :
				if self.signalP :
					check=self.computeSignap()
					if not check:
						print "Error with SignalP, -maybe is not installed on your computer- computing continue to the next step- "
						self.signalP=False
				if self.Hmmer:
					check=self.computeHmmer()
					if not check:
						print "Error with HMMER, -maybe is not installed on your computer- computing continue to the next step- "
						self.Hmmer=False
				if self.remove :
					check=os.remove(self.fastaFile)				
		
	#SignalP methods	
	def computeSignap(self):
		file=os.popen(self.SignalPbinCommand+" -t "+self.organismType+" "+self.fastaFile+" > "+self.tempFile)
		status=file.close()
		if not status :
			self.parseSignalPResults()	
		return False
	def parseSignalPResults(self):
		file=open(self.tempFile,'r')
		contents=file.read()
		lines=re.split('\n',contents)
		header=re.compile("^#(.*)")
		for i,line in enumerate(lines) :
			if i != len(lines)-1 :
				match = re.search(header,line)
				if not match :
					fields=re.split('\s+',line)
					if fields[self.isAsignalPeptidInd] =="Y":
						self.nbOfSignal+=1
		file.close()
	
		
	#Hmm Methods
	def computeHmmer(self):
		file=os.popen(self.HMMscanbinCommand+" --domtblout="+self.tempFile+" -o "+self.hmmScanOutput+" --cut_tc "+self.HMMdb+" "+self.fastaFile)
		status=file.close()
		if not status :
			self.parseHMMERResults()
			os.remove(self.hmmScanOutput)
			return True
		else :
			return False
	def parseHMMERResults(self):
		file=open(self.tempFile,'r')
		contents=file.read()
		lines=re.split('\n',contents)
		header=re.compile("^#(.*)")
		tmpID=""
		for i,line in enumerate(lines) :
			if i != len(lines)-1 :
				match = re.search(header,line)
				if not match : #skip the headers
					fields=re.split('\s+',line)
					queryID=fields[self.indHmmerOutputFile]
					if queryID != tmpID:
						self.nbOfPfamDomains+=1
						tmpID=queryID		
		file.close()
		
		
	#Needle methods
	def computeNeedle(self,CDSone,CDStwo,clean=True):
		fileOne=self.extProtSeqFromCDS(CDSone,suffixe="_protOne.fasta")
		fileTwo=self.extProtSeqFromCDS(CDStwo,suffixe="_protTwo.fasta")
		#Execute Needle
		try:
			needle_cline = NeedleCommandline(asequence=fileOne, bsequence=fileTwo,gapopen=10, gapextend=0.5, outfile=self.needleOutputFile)
			needle_cline()
			#Parse Results
			needle = open(self.needleOutputFile)
			for line in needle:
				if line.find("Identity") > 0:
					self.needleIdentity = int(line[line.find(":") + 1:line.find("/")].strip())
					self.needleIdentityPercent = float(line[line.find("(") + 1:line.find("%")])
				if line.find("Similarity") > 0:
					self.needleSimilarity = int(line[line.find(":") + 1:line.find("/")].strip())
					self.needleSimilarityPercent = float(line[line.find("(") + 1:line.find("%")])
			needle.close()
			# https://www.biostars.org/p/74060/ thanks
			if clean : #Delete files
				os.remove(fileOne)
				os.remove(fileTwo)
				os.remove(self.needleOutputFile)
		except Application.ApplicationError :
			print "Warning : Error with Needle, the similarity will be 0."
