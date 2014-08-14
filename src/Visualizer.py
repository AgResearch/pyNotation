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
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.colors import red,black,orange,green,blue,brown,purple,pink,white,grey
from Bio import SeqIO
from ExternalToolsManager import ExternalToolsManager
import sys

minSimilarityScore=15
redMinValue=0.90
greenMinValue=0.90
blueMinValue=0.90
redMaxValue=0.98
greenMaxValue=0.98
blueMaxValue=0.98

#########################
######### Class #########
#########################
class Similarity(object):
	def __init__(self,cdsAlpha,cdsBeta,similarity,betaTrack):
		self.cdsAlpha=cdsAlpha
		self.cdsBeta=cdsBeta
		if similarity < minSimilarityScore :
			similarity=minSimilarityScore+0.05
		self.similarity=similarity
		self.betaTrack=betaTrack
	def __str__(self):
		return str(self.cdsAlpha)+"\n-\n"+str(self.cdsBeta)+"-> "+str(self.similarity)
		

class Track(object): 
	backgroundColor=colors.Color(0.96,0.96,0.96)
	maxSize=0
	minSize=sys.maxint
	maxLabelWord=3
	angle=[0,180]
	strand=[1,-1]
	
	def __init__(self,trackName):
		self.trackName=trackName
		self.gdTrack=GenomeDiagram.Track(greytrack=True,name=self.trackName,greytrack_labels=1,greytrack_fontsize=33)
		self.gdFeature=GenomeDiagram.FeatureSet()
		self.size=0
		self.nbFeats=0
	def __str__(self):
		return self.name+" "+str(len(self.gdTrack.get_sets()))+" sets."
	
class RefTrack(Track):	
	allColors=[red,black,orange,green,blue,brown,purple,pink]
	def __init__(self,trackName):
		Track.__init__(self,trackName[0:10])
		self.c=0
		self.colorHash={}
		
	def addFeat(self,feat):
		tempCol=self.allColors[self.c%len(self.allColors)]
		self.c+=1
		self.colorHash[feat]=tempCol
		labelTab=[]
		try :
			labelTab=feat.qualifiers['product'][0].split(" ")
		except KeyError:
			print " NO PRODUCT FOUND FOR "
			print feat
			labelTab[0]="No name"
		labelName=""
		if len(labelTab)<= Track.maxLabelWord :
			for word in labelTab :
				labelName+=word+" "
		else :
			labelName=labelTab[0]+" "+labelTab[1]+" "+labelTab[2]+" "
				
		labelName=labelName[0:len(labelName)-1]+" \n "+str(feat.location.start)+" - "+str(feat.location.end) #skip the final space and add location
		if self.nbFeats==0:
			self.gdFeature.add_feature(feat,color=Track.backgroundColor,sigil="ARROW",name=labelName,label_position="start",label_angle=Track.angle[self.nbFeats%2],label=True,strand=Track.strand[self.nbFeats%2])
		else :
			self.gdFeature.add_feature(feat,color=Track.backgroundColor,sigil="ARROW",name=labelName,label_position="middle",label_angle=Track.angle[self.nbFeats%2],label=True,strand=Track.strand[self.nbFeats%2])
		
		if feat.location.end > Track.maxSize :
			Track.maxSize=feat.location.end
		if feat.location.start < Track.minSize :
			Track.minSize=feat.location.start 
		self.nbFeats+=1
		if feat.strand==1:
			self.gdFeature.add_feature(feat,color=tempCol,sigil="ARROW",label_position="middle",label_angle=0,label=False)
		else :
			self.gdFeature.add_feature(feat,color=tempCol,sigil="ARROW",label_position="middle",label_angle=180,label=False)
		self.gdTrack.add_set(self.gdFeature)
			
	
class QueryTrack(Track):
	Level=1
	def __init__(self,name):
		Track.__init__(self,name)
		self.similarities=[]
		self.LevelToDraw=QueryTrack.Level
		QueryTrack.Level+=1
		self.diff=0

	def flip(self):
		QueryTrack.Level=1
		self.diff=0
		Track.minSize=sys.maxint
		Track.maxSize=0
		
	def scale(self,firstFeat):
		self.diff=firstFeat.location.start-Track.minSize
		
	def addFeat(self,feat,colorMax,percent):
		labelTab=[]
		try :
			labelTab=feat.qualifiers['product'][0].split(" ")
		except KeyError:
			print " NO PRODUCT FOUND FOR "
			print feat
			labelTab[0]="No name"
		labelName=""
		if len(labelTab)<= Track.maxLabelWord :
			for word in labelTab :
				labelName+=word+" "
		else :
			labelName=labelTab[0]+" "+labelTab[1]+" "+labelTab[2]+" "
				
		labelName=labelName[0:len(labelName)-1]+" \n "+str(feat.location.start)+" - "+str(feat.location.end) #skip the final space and add location
		
		#change location
		newStart=feat.location.start-self.diff
		newEnd=feat.location.end-self.diff
		if newEnd > Track.maxSize :
			Track.maxSize=newEnd
		newLocation = FeatureLocation(newStart,newEnd,feat.strand)
		feat.location=newLocation
		if self.nbFeats==0:
			self.gdFeature.add_feature(feat,color=Track.backgroundColor,sigil="ARROW",name=labelName,label_position="start",label_angle=Track.angle[self.nbFeats%2],label=True,strand=Track.strand[self.nbFeats%2])
		else :
			self.gdFeature.add_feature(feat,color=Track.backgroundColor,sigil="ARROW",name=labelName,label_position="middle",label_angle=Track.angle[self.nbFeats%2],label=True,strand=Track.strand[self.nbFeats%2])
		if feat.strand==1:			
			self.gdFeature.add_feature(feat,border=colorMax,color=colors.linearlyInterpolatedColor(white,colorMax,minSimilarityScore,100,percent),sigil="ARROW",name=feat.qualifiers['product'][0][0:11].replace(" ","_")+" \n "+str(feat.location.start)+" - "+str(feat.location.end),label_position="middle",label_angle=0,label=False)
		else :
			self.gdFeature.add_feature(feat,border=colorMax,color=colors.linearlyInterpolatedColor(white,colorMax,minSimilarityScore,100,percent),sigil="ARROW",label_position="middle",label_angle=180,label=False)
		self.nbFeats+=1
		self.gdTrack.add_set(self.gdFeature)
	def addSimilarity(self,similarity):
		self.similarities.append(similarity)
	def to_string(self):
		return self.gdTrack.to_string()

class TrackManager(object):
	def __init__(self):
		self.tracks=[]
		self.gdd = GenomeDiagram.Diagram('Diagram') #Is this name useful ?
	def setRef(self,trackName): 
		self.refTrack=RefTrack(trackName)
	def addTrack(self,name): 
		track=QueryTrack(name)
		self.tracks.append(track)
		return track
	def addRefFeature(self,feat): #feat supposed to be a Feature object from BioPython
		self.refTrack.addFeat(feat)
	def addQueryFeature(self,queryTrack,similarity):
		for track in self.tracks:
			if track==queryTrack:
				track.addFeat(similarity.cdsBeta,self.refTrack.colorHash[similarity.cdsAlpha],similarity.similarity)
				track.addSimilarity(similarity) 
				return True
		print "The track "+queryTrack.to_string()+" was not found.\nSystem will exit" 
		sys.exit(1) 
	def draw(self,outputname):
		#self.refTrack could be None
		if not self.refTrack :
			print "Ref track must be define before drawing.\nSystem will exit" 
			sys.exit(1)
		self.gdd.add_track(self.refTrack.gdTrack,len(self.tracks)+1) 
		for track in self.tracks :
			self.gdd.add_track(track.gdTrack,track.LevelToDraw)
		self.gdd.draw(start=Track.minSize-100,end=Track.maxSize+100,pagesize="A3",format="linear",fragments=1) #options to add
		self.gdd.write(outputname,"PDF")

	def scale(self,t,feat):
		for track in self.tracks:
			if track==t:
				t.scale(feat)
				
	def removeTracks(self):
		self.gdd.del_track(len(self.tracks)+1) #remove the ref Track
		for i,track in enumerate(self.tracks) :
			track.flip()
			self.gdd.del_track(track.LevelToDraw) #remove the query
		t=self.gdd.get_levels()
		self.tracks=[]
		
		
class VisuManager(object):
	def __init__(self,interval,outputfileName):
		self.outputfileName=outputfileName
		if interval %2 ==0: #interval must be odd
			self.interval=interval+1
		else :
			self.interval=interval
		self.refCDSarray=[]
		self.tm=TrackManager()
		self.etm=ExternalToolsManager(None)
	
	def clean(self):
		self.tm.removeTracks()

	def preComputeRecord(self,record) :
		temp=[]
		for feat in record.features :
			if feat.type=="CDS":
				temp.append(feat)
		return temp

	def initializeRefTrack(self,name):
		self.tm.setRef(name)
		for feat in self.refCDSarray:
			self.tm.addRefFeature(feat)
		
	def draw(self):
		self.tm.draw(self.outputfileName)
			
	def addRefCDS(self,locus,cds): #return False if locus or cds are not found in file
	#Locus should be Bio.SeqRecord.SeqRecord
	#CDS should be Bio.SeqFeature.SeqFeature
		self.refCDS=cds
		feaTab=self.preComputeRecord(locus)
		if ( len(feaTab)< self.interval ):
			for feat in feaTab:
				self.refCDSarray.append(feat)
			self.initializeRefTrack(locus.name) 
			return True
		for i,feature in enumerate(feaTab):
			j=0
			del self.refCDSarray[:]
			# Stocking the CDS neighborhood
			while j<self.interval :
				self.refCDSarray.append(locus.features[i+j])
				j+=1
			#Checking if the CDS is on the table refCDSarray	
			if self.refCDSarray[self.interval/2]==cds:
				self.initializeRefTrack(locus.name)
				return True
			elif i==0 : #The CDS could be in [0:self.interval/2] i.e. at the beginning of the locus
				for k in range(0,self.interval/2+1):
					if self.refCDSarray[k]==cds:
						self.initializeRefTrack(locus.name) 
						return True
			elif i + self.interval == len(feaTab): #The CDS could be in [self.interval/2:len(self.refCDSarray)] i.e. at the end of the locus
				if cds in self.refCDSarray :
					self.initializeRefTrack(locus.name) 
					return True
		return False
		
	def addTrack(self,name,tab):	
		t=self.tm.addTrack(name)
		for i,otherTrack in enumerate(tab):
			max=0
			for reFeat in self.refCDSarray:
				self.etm.computeNeedle(reFeat,otherTrack)
				if ( self.etm.needleSimilarityPercent  > max ):
					max = self.etm.needleSimilarityPercent
					refTemp=reFeat
				if (max==100.0): #optimization
					break 
			if i==0:
				self.tm.scale(t,otherTrack)
			s=Similarity(refTemp,otherTrack,max,t)
			self.tm.addQueryFeature(t,s)
		
	def addQueryCDS(self,locus,cds):
		feaTab=self.preComputeRecord(locus)
		neighborhood=[]
		if ( feaTab< self.interval ):
			for feat in feaTab:
				neighborhood.append(feat)
			self.addTrack(locus.name,neighborhood) 
			return True
		for i,feature in enumerate(feaTab):
			j=0
			del neighborhood[:]
			# Stocking the CDS neighborhood
			while j<self.interval :
				neighborhood.append(feaTab[i+j])
				j+=1
			if neighborhood[self.interval/2]==cds :
				self.addTrack(locus.name,neighborhood)
				return True
			elif i==0 : #The CDS could be in [0:self.interval/2] i.e. at the beginning of the locus
				for k in range(0,self.interval/2+1):
					if neighborhood[k]==cds:
						self.addTrack(locus.name,neighborhood) 
						return True
			elif i + self.interval==len(feaTab): #The CDS could be in [self.interval/2:len(neighborhood)] i.e. at the end of the locus
				if cds in neighborhood :
					self.addTrack(locus.name,neighborhood) 
					return True
				else :
					print "CDS"
					print str(cds)
					print " not found in"
					print str(locus)
					return False
			if len(feaTab)==self.interval:
				self.addTrack(locus.name,neighborhood) 
				return True
		return False
