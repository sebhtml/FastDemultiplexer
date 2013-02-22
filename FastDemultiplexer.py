#!/usr/bin/python
# encoding: UTF-8
# author: Sébastien Boisvert
# this is GPL code
'''
	FastDemultiplexer: a better demultiplexer for Illumina HiSeq
sequencers
	Copyright (C) 2011, 2012 Sébastien Boisvert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#import profile
import gzip
import zlib
import sys
import os
import os.path

if len(sys.argv)!=4:
	print("usage")
	print("FastDemultilexer SampleSheet.csv Project_XYZ/Sample_lane1 Demultiplexed")
	sys.exit(0)

class GzFileReader:
	def __init__(self,file):
		self.m_file=gzip.open(file)

	def readline(self):
		return self.m_file.readline()
	def close(self):
		self.close()

class Entry:
	def __init__(self,project,sample,index1,index2):
		self.m_project=project
		self.m_sample=sample
		self.m_index1=index1
		self.m_index2=index2

	def getProject(self):
		return self.m_project
	def getSample(self):
		return self.m_sample
	def getIndex1(self):
		return self.m_index1
	def getIndex2(self):
		return self.m_index2

class SampleSheet:
	def __init__(self,sampleSheet,lane):
		# C0947ACXX,4,CQDM1-1,No,TAAGGCGA-TAGATCGC,P2J0-1,N,PE_indexing,LR,CQDM
		projectColumn=9
		sampleColumn=2
		indexColumn=4
		laneColumn=1

		self.m_entries=[]

		for line in open(sampleSheet):
			tokens=line.split(",")
			if len(tokens)<4:
				continue

			theLane=tokens[laneColumn]

			if lane!=theLane:
				continue

			project=tokens[projectColumn].strip()
			sample=tokens[sampleColumn]
			index=tokens[indexColumn]

			tokens2=index.split("-")
			index1=tokens2[0]
			index2=""
			if len(tokens2)==2:
				index2=tokens2[1]

			entry=Entry(project,sample,index1,index2)
			
			self.m_entries.append(entry)

			self.m_maxMismatches1=len(index1)/2
			self.m_index1Length=len(index1)
			self.m_maxMismatches2=len(index2)/2
			self.m_index2Length=len(index2)


		if len(self.m_entries)==0:
			print("Error: the SampleSheet does not contain entries for the lane provided.")
			
		self.makeIndex()

	def getErrorList(self,base):
		list=[]

		i=0
		changes=['A','T','C','G','N']

		theLength=len(base)

		while i<theLength:
			actual=base[i]
			for j in changes:
				if j==actual:
					continue
				before=base[0:i]
				after=base[(i+1):(theLength)]
				newSequence=before+j+after
				list.append(newSequence)

			i+=1

		return list

	def makeIndex(self):
		self.m_index={}
		for entry in self.m_entries:
			key=entry.getIndex1()+entry.getIndex2()
			value=[entry.getProject(),entry.getSample()]
			self.m_index[key]=value

			# generate things with 1 error in index1
			index1ErrorList=self.getErrorList(entry.getIndex1())
	
			for i in index1ErrorList:
				key=i+entry.getIndex2()
				self.m_index[key]=value

			# generate things with 1 error in index2
			index2ErrorList=self.getErrorList(entry.getIndex2())

			for i in index2ErrorList:
				key=entry.getIndex1()+i
				self.m_index[key]=value

			# generate things with 1 error in index1 and 1 error in index2
			for i in index1ErrorList:
				for j in index2ErrorList:
					key=i+j
					self.m_index[key]=value

		print("IndexSize= "+str(len(self.m_index)))
		print("Index1Length= "+str(self.m_index1Length))
		print("AllowedMismatchesInIndex1= "+str(self.m_maxMismatches1))
		print("Index2Length= "+str(self.m_index2Length))
		print("AllowedMismatchesInIndex2= "+str(self.m_maxMismatches2))

	def compare(self,sequence1,sequence2):
		score=0
		i=0
		len1=len(sequence1)
		len2=len(sequence2)

		while i<len1 and i<len2:
			if sequence1[i]!=sequence2[i]:
				score+=1
			i+=1
		
		return score

	def classify(self,index1,index2,lane):
		key=index1+index2

		if key in self.m_index:
			return self.m_index[key]

		for entry in self.m_entries:
			score1=self.compare(entry.getIndex1(),index1)
			score2=self.compare(entry.getIndex2(),index2)


			if score1<=self.m_maxMismatches1 and score2<=self.m_maxMismatches2:
				return [entry.getProject(),entry.getSample()]

		return ["Undetermined_indices","Sample_lane"+lane]

class Sequence:
	def __init__(self,line1,line2,line3,line4):
		self.m_line1=line1
		self.m_line2=line2
		self.m_line3=line3
		self.m_line4=line4

	def getLine1(self):
		return self.m_line1
	def getLine2(self):
		return self.m_line2
	def getLine3(self):
		return self.m_line3
	def getLine4(self):
		return self.m_line4

class FileReader:
	def __init__(self,filePath):
		print("Opening "+filePath)
		if filePath.find(".gz")>=0:
			self.m_file=GzFileReader(filePath)
		else:
			self.m_file=open(filePath)
		
		self.m_buffer=self.m_file.readline().strip()

	def hasNext(self):
		return len(self.m_buffer)>0

	def getNext(self):
		sequence=Sequence(self.m_buffer,self.m_file.readline().strip(),self.m_file.readline().strip(),self.m_file.readline().strip())
		
		self.m_buffer=self.m_file.readline().strip()

		return sequence

class InputDirectory:
	def __init__(self,laneDirectory):
		self.m_directory=laneDirectory

		self.m_r1Files=[]
		self.m_r2Files=[]
		self.m_r3Files=[]
		self.m_r4Files=[]

		for i in os.listdir(self.m_directory):
			if i.find("_R1_")>=0:
				self.m_r1Files.append(i)
				self.m_r2Files.append(i.replace("_R1_","_R2_"))
				self.m_r3Files.append(i.replace("_R1_","_R3_"))
				self.m_r4Files.append(i.replace("_R1_","_R4_"))

		self.m_current=0

		self.setReadersToCurrent()

		while not self.m_reader1.hasNext() and self.m_current<len(self.m_r1Files):
			self.m_current+=1
			self.setReadersToCurrent()

	def setReadersToCurrent(self):
		if not self.m_current<len(self.m_r1Files):
			return

		self.m_reader1=FileReader(self.m_directory+"/"+self.m_r1Files[self.m_current])
		self.m_reader2=FileReader(self.m_directory+"/"+self.m_r2Files[self.m_current])
		self.m_reader3=FileReader(self.m_directory+"/"+self.m_r3Files[self.m_current])
		self.m_reader4=FileReader(self.m_directory+"/"+self.m_r4Files[self.m_current])
			

	def hasNext(self):
		if self.m_current>=len(self.m_r1Files):
			return False

		while not self.m_reader1.hasNext() and self.m_current<len(self.m_r1Files):
			print("opening next files")
			self.m_current+=1
			self.setReadersToCurrent()

		return self.m_reader1.hasNext()

	def getNext(self):
		return [self.m_reader1.getNext(),self.m_reader2.getNext(),self.m_reader3.getNext(),self.m_reader4.getNext()]

class FileWriter:
	def __init__(self,name):
		if name.find(".gz")>=0:
			self.m_file=gzip.open(name,"w")
		else:
			self.m_file=open(name,"w")
	def write(self,data):
		self.m_file.write(data)
	def close(self):
		self.m_file.close()

class OutputDirectory:
	def __init__(self,outputDirectory,maxInFile):
		self.m_directory=outputDirectory
		self.m_max=maxInFile

		self.makeDirectory(self.m_directory)
		self.m_files1={}
		self.m_files2={}

		self.m_counts={}
		self.m_currentNumbers={}

	def makeDirectory(self,name):
		if not os.path.exists(name):
			os.mkdir(name)

	def closeFiles(self):
		for i in self.m_files1.items():
			i[1].close()

		for i in self.m_files2.items():
			i[1].close()

	def write(self,project,sample,lane,sequenceTuple):

		key=project+sample+lane

		if (key not in self.m_files1) or self.m_counts[key]==self.m_max:

			#close the old files
			if key in self.m_files1 and self.m_counts[key]==self.m_max:
				self.m_files1[key].close()
				self.m_files2[key].close()
				self.m_currentNumbers[key]+=1
			else:
				self.m_currentNumbers[key]=1

			self.m_counts[key]=0

			projectDir=project
			sampleDir=sample

			if project!="Undetermined_indices":
				projectDir="Project_"+project
				sampleDir="Sample_"+sample

			self.makeDirectory(self.m_directory+"/"+projectDir)
			self.makeDirectory(self.m_directory+"/"+projectDir+"/"+sampleDir)

			file1=self.m_directory+"/"+projectDir+"/"+sampleDir+"/"+sample+"_Lane"+lane+"_R1_"+str(self.m_currentNumbers[key])+".fastq"
			file2=self.m_directory+"/"+projectDir+"/"+sampleDir+"/"+sample+"_Lane"+lane+"_R2_"+str(self.m_currentNumbers[key])+".fastq"

			self.m_files1[key]=FileWriter(file1)
			self.m_files2[key]=FileWriter(file2)

		f1=self.m_files1[key]
		f1.write(sequenceTuple[0].getLine1()+"\n"+sequenceTuple[0].getLine2()+"\n"+sequenceTuple[0].getLine3()+"\n"+sequenceTuple[0].getLine4()+"\n")

		f2=self.m_files2[key]
		f2.write(sequenceTuple[1].getLine1()+"\n"+sequenceTuple[1].getLine2()+"\n"+sequenceTuple[1].getLine3()+"\n"+sequenceTuple[1].getLine4()+"\n")

		self.m_counts[key]+=1



class Demultiplexer:
	def __init__(self,sampleSheet,inputDirectoryPath,outputDirectoryPath,lane):
		sheet=SampleSheet(sampleSheet,lane)
		inputDirectory=InputDirectory(inputDirectoryPath)

		maxInFile=4000000

		outputDirectory=OutputDirectory(outputDirectoryPath,maxInFile)

		self.m_processed=0

		self.m_stats={}

		while inputDirectory.hasNext():
			sequenceTuple=inputDirectory.getNext()
			index1=sequenceTuple[1].getLine2()
			index2=sequenceTuple[2].getLine2()

			[project,sample]=sheet.classify(index1,index2,lane)

			outputDirectory.write(project,sample,lane,[sequenceTuple[0],sequenceTuple[3]])
		
			self.m_processed+=1

			projectDir=project
			sampleDir=sample

			if project!="Undetermined_indices":
				projectDir="Project_"+project
				sampleDir="Sample_"+sample

			if projectDir not in self.m_stats:
				self.m_stats[projectDir]={}

			if sampleDir not in self.m_stats[projectDir]:
				self.m_stats[projectDir][sampleDir]=0

			self.m_stats[projectDir][sampleDir]+=1

			if self.m_processed%10000==0:
				self.printStatus()

		outputDirectory.closeFiles()

		self.printStatus()

	def printStatus(self):
		print("[Status]")

		for i in self.m_stats.items():
			for j in i[1].items():
				percent=100.0*j[1]/self.m_processed
				print(i[0]+"	"+j[0]+"	"+str(j[1])+"	"+str(percent)+"%")

		print("All	All	"+str(self.m_processed)+"	100.00%")
		sys.stdout.flush()

def main():
	sheet=sys.argv[1]
	inputDir=sys.argv[2]
	outputDir=sys.argv[3]

	lane=inputDir.split("lane")[1].strip().replace("/","")

	demultiplexer=Demultiplexer(sheet,inputDir,outputDir,lane)

if __name__=="__main__":
#	doProfiling=False

#	if doProfiling:
#		profile.run('main()')
#	else:

	main()


