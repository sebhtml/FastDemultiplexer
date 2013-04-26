#!/usr/bin/python
# encoding: UTF-8
# author: Sébastien Boisvert
# this is GPL code
'''
	FastDemultiplexer: a better demultiplexer for Illumina HiSeq
sequencers
	Copyright (C) 2011, 2012, 2013 Sébastien Boisvert

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

import profile
import gzip
import zlib
import sys
import os
import os.path

class GzFileReader:
	def __init__(self,file):
		self.m_file=gzip.open(file)

	def readline(self):
		return self.m_file.readline()
	def close(self):
		self.m_file.close()

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
		self.m_debug = False
		self.m_conservativeMode = False
		#self.m_debug = True
		self.m_error = False
		# C0947ACXX,4,CQDM1-1,No,TAAGGCGA-TAGATCGC,P2J0-1,N,PE_indexing,LR,CQDM
		projectColumn=9
		sampleColumn=2
		indexColumn=4
		laneColumn=1

		self.m_entries=[]

		allocatedKeys = {}

		for line in open(sampleSheet):
			if line[0] == '#':
				continue

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

			lengthOfIndex1 = len(index1)
			lengthOfIndex2 = len(index2)

			if lengthOfIndex1 == 0 and lengthOfIndex2 == 0:
				print("Warning: " + sample + " has no index in sheet")
				continue

			key = index1 + index2

			if key in allocatedKeys:
				print("Warning: " + sample + " uses a key already in use by another sample")
				# don't continue on this
	
			allocatedKeys[key]=1
			
			self.m_entries.append(entry)

			self.m_index1Length = lengthOfIndex1
			self.m_index2Length = lengthOfIndex2

		if len(self.m_entries)==0:
			print("Error: the SampleSheet does not contain entries for the lane provided.")
			print("Lane is: " + lane)
			self.m_error = True
			return

		self.makeIndex()

	def hasError(self):
		return self.m_error

# returns a list of sequence based on template with mismatches when compared to origin
# this code actually returns duplicates too.
# we don't really care about that.
	def getErrorList(self, origin, template, mismatches, list):

		if mismatches == 0:
			list.append(template)
			return

		position = 0
		changes=['A','T','C','G','N']

		theLength=len(template)

		while position < theLength:
			actual=origin[position]
			for j in changes:

				# the change restore the original thing...
				if j==actual:
					continue

				before=template[0:position]
				after=template[(position+1):(theLength)]
				newSequence=before+j+after

				self.getErrorList(origin, newSequence, mismatches - 1, list)

			position += 1

	def makeIndex(self):

		print("[makeIndex] Index1Length= "+str(self.m_index1Length))
		print("[makeIndex] Index2Length= "+str(self.m_index2Length))

		self.m_index={}
		for entry in self.m_entries:
			rule = [0, 1, 2, 3]

			for i in rule:
				for j in rule:
					if i > 1 and j > 1:
						continue

					self.addEntriesInIndex(self.m_index, entry.getIndex1(), entry.getIndex2(), i, j, entry)

			# 2, 2 is too prohibitive
			#addEntriesInIndex(self.m_index, entry.getIndex1(), entry.getIndex2(), 2, 2, entry)


		print("[makeIndex] IndexSize= "+str(len(self.m_index)))
		if self.m_debug:
			for i in self.m_index.items():
				print("[makeIndex] " + i[0] + " ---> " + str(i[1]))

		return True

	def addEntriesInIndex(self, index, sequence1, sequence2, mismatches1, mismatches2, target):
		objects1 = []
		self.getErrorList(sequence1, sequence1, mismatches1, objects1)
		objects2 = []
		self.getErrorList(sequence2, sequence2, mismatches2, objects2)

		count = 0
		for i in objects1:
			for j in objects2:
				key = i + j
				index[key] = target
				count += 1

		print("[addEntriesInIndex] added " + str(count) + " entries with mismatch configuration " + str(mismatches1)+ "," + str(mismatches2))

	def getMismatches(self,sequence1,sequence2):

		score=0
		i=0
		len1=len(sequence1)
		len2=len(sequence2)

		while i<len1 and i<len2:
			if sequence1[i]!=sequence2[i]:
				score+=1
			i+=1

		if self.m_debug:
			print("[getMismatches] " + sequence1 + " " + sequence2 + " " + str(score))

		return score

	def classify(self,index1,index2,lane):

		result = self.classifyWithTheIndex(index1, index2, lane)

		if result != None:
			return result

		return self.classifyWithBruteForce(index1, index2, lane)

	def classifyWithTheIndex(self, index1, index2, lane):
		if self.m_debug:
			print("[classify] index1: " + index1 + " index2: " + index2)

		key=index1+index2

		if self.m_index2Length == 0:
			key = index1

		# some kits will include an extra A between the tag and the sequence
		if len(index1) == self.m_index1Length + 1:
			key = key[0:self.m_index1Length]

		# use the hash table to classify it in a
		# fast way
		if key in self.m_index:
			if self.m_debug:
				print("[classify] using fast-path with index")
			return self.m_index[key]

		return None

	def classifyWithBruteForce(self, index1, index2, lane):

		best1 = 9999
		best2 = 9999
		bestEntry = None

		for entry in self.m_entries:
			score1=self.getMismatches(entry.getIndex1(),index1)
			score2=self.getMismatches(entry.getIndex2(),index2)

			if self.m_index2Length == 0:
				score2 = 0

			# both most be at least as good
			# and at least one of them must be better
			if ( ( score1 <= best1 and score2 <= best2 ) and 
				( score1 < best1 or score2 < best2 )):

				best1 = score1
				best2 = score2
				bestEntry = entry
				
				if self.m_debug:
					print("[classify] new best is " + entry.getSample())

			# at least two entries have the same number of
			# mismatches
			elif score1 == best1 and score2 == best2:
				bestEntry = None

		if best1 >= self.m_index1Length/2 and self.m_index1Length != 0 and self.m_conservativeMode:
			bestEntry = None

		if best2 >= self.m_index2Length/2 and self.m_index2Length != 0 and self.m_conservativeMode:
			bestEntry = None

		if self.m_debug:
			print("best1: " + str(best1) + " best2: " + str(best2) + " with " + entry.getSample())

		return bestEntry

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
		if filePath.find(".gz")>=0:
			self.m_file=GzFileReader(filePath)
		else:
			self.m_file=open(filePath)

		self.m_buffer=self.m_file.readline().strip()

	def hasNext(self):
		result = len(self.m_buffer)>0
		if not result:
			self.m_file.close()

		return result

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

		self.m_hasError = False
		self.setReadersToCurrent()

		if self.m_hasError:
			print("Error: No sequence file found !")
			return

		while not self.m_reader1.hasNext() and self.m_current<len(self.m_r1Files):
			self.m_current+=1
			self.setReadersToCurrent()

	def hasError(self):
		return self.m_hasError

	def setReadersToCurrent(self):
		if not self.m_current<len(self.m_r1Files):
			self.m_hasError = True
			return

		self.m_reader1=FileReader(self.m_directory+"/"+self.m_r1Files[self.m_current])
		self.m_reader2=FileReader(self.m_directory+"/"+self.m_r2Files[self.m_current])
		self.m_reader3=FileReader(self.m_directory+"/"+self.m_r3Files[self.m_current])
		self.m_reader4=FileReader(self.m_directory+"/"+self.m_r4Files[self.m_current])


	def hasNext(self):
		if self.m_current>=len(self.m_r1Files):
			return False

		while not self.m_reader1.hasNext() and self.m_current<len(self.m_r1Files):
			self.m_current+=1
			self.setReadersToCurrent()

		return self.m_reader1.hasNext()

	def getNext(self):
		return [self.m_reader1.getNext(),self.m_reader2.getNext(),self.m_reader3.getNext(),self.m_reader4.getNext()]

class FileWriter:
	def __init__(self,name):
		if name.find(".fastq.gz")>=0:
			self.m_file=gzip.open(name,"w")
		else:
			self.m_file=open(name,"w")
	def write(self,data):
		self.m_file.write(data)
	def close(self):
		self.m_file.close()

class OutputDirectory:
	def __init__(self,outputDirectory):
		self.m_debug = False
		self.m_directory=outputDirectory
		self.m_maximumNumberOfSequencesPerFile = 4000000

		self.m_maximumNumberOfStagedObjects = 50000

		self.makeDirectory(self.m_directory)
		self.m_files1={}
		self.m_files2={}

		self.m_stagingArea1 = {}
		self.m_stagingArea2 = {}

		self.m_counts={}
		self.m_currentNumbers={}

	def makeDirectory(self,name):
		if not os.path.exists(name):
			os.mkdir(name)

	def closeFiles(self):
		for i in self.m_files1.items():
			key = i[0]
			self.flushWriteOperationsForKey(key, True)
			self.m_files1[key].close()
			self.m_files2[key].close()


	def write(self,project,sample,lane,sequenceTuple):

		key=project+sample+lane

		if (key not in self.m_files1) or self.m_counts[key]==self.m_maximumNumberOfSequencesPerFile:

			#close the old files
			if key in self.m_files1 and self.m_counts[key]==self.m_maximumNumberOfSequencesPerFile:
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

			compressFiles = True

			if compressFiles:
				file1 += ".gz"
				file2 += ".gz"

			self.m_files1[key]=FileWriter(file1)

			if self.m_debug:
				print("[write] opening " + key + " --> " + file1)
			self.m_files2[key]=FileWriter(file2)

			self.m_stagingArea1[key] = []
			self.m_stagingArea2[key] = []

		self.m_stagingArea1[key].append(sequenceTuple[0])
		self.m_stagingArea2[key].append(sequenceTuple[1])

		self.flushWriteOperationsForKey(key, False)

	def flushWriteOperationsForKey(self, key, forceOperation):

		entryIterator = 0

		f1=self.m_files1[key]
		f2=self.m_files2[key]

		stagedEntries = len(self.m_stagingArea1[key])

		proceed = False

		if stagedEntries  >= self.m_maximumNumberOfStagedObjects or forceOperation:
			proceed = True

		if stagedEntries == 0:
			proceed = False

		if not proceed:
			if self.m_debug and False:
				print("[flushWriteOperationsForKey] key: " + key + " stagedEntries: " + str(stagedEntries))
			return False

		buffer1 = ""
		buffer2 = ""

		if self.m_debug:
			print("[flushWriteOperationsForKey] flushing " + str(stagedEntries) + " for " + key)

		while entryIterator < stagedEntries:
			entry1 = self.m_stagingArea1[key][entryIterator]
			entry2 = self.m_stagingArea2[key][entryIterator]
			line1 = entry1.getLine1()+"\n"+entry1.getLine2()+"\n"+entry1.getLine3()+"\n"+entry1.getLine4()+"\n"
			buffer1 += line1 
			line2 = entry2.getLine1()+"\n"+entry2.getLine2()+"\n"+entry2.getLine3()+"\n"+entry2.getLine4()+"\n"
			buffer2 += line2

			self.m_counts[key]+=1
			entryIterator += 1

		f1.write(buffer1)
		f2.write(buffer2)

		self.m_stagingArea1[key] = []
		self.m_stagingArea2[key] = []

		return True

class Demultiplexer:
	def __init__(self,sampleSheet,inputDirectoryPath,outputDirectoryPath,lane):
		sheet=SampleSheet(sampleSheet,lane)
		self.m_debug = False
		#self.m_debug = True

		if sheet.hasError():
			return

		inputDirectory=InputDirectory(inputDirectoryPath)

		if inputDirectory.hasError():
			return

		outputDirectory=OutputDirectory(outputDirectoryPath)

		self.m_processed=0

		self.m_stats={}

		while inputDirectory.hasNext():
			sequenceTuple=inputDirectory.getNext()
			index1=sequenceTuple[1].getLine2()
			index2=sequenceTuple[2].getLine2()

			project = "Undetermined_indices"
			sample = "Sample_lane" + lane
			entry = sheet.classify(index1,index2,lane)

			if entry != None:
				project = entry.getProject()
				sample = entry.getSample()

			if self.m_debug:
				print("index1= " + index1 + " index2= " + index2 + " lane= " + lane + " result: " + sample)

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

		print("Project	Sample	Count	Percentage")

		for i in self.m_stats.items():
			for j in i[1].items():
				percent=100.0*j[1]/self.m_processed
				print(i[0]+"	"+j[0]+"	"+str(j[1])+"	"+str(percent)+"%")

		print("*	*	"+str(self.m_processed)+"	100.00%")
		sys.stdout.flush()

def main():
	if len(sys.argv)!=5:
		print("usage")
		print("FastDemultilexer SampleSheet.csv lane Project_XYZ/Sample_lane1 Demultiplexed")
		sys.exit(0)

	arguments = sys.argv
	parameterIndex = 0
	parameterIndex += 1

	sheet = arguments[parameterIndex]
	parameterIndex += 1

	lane = arguments[parameterIndex]
	parameterIndex += 1

	inputDir = arguments[parameterIndex]
	parameterIndex += 1

	outputDir = arguments[parameterIndex]
	parameterIndex += 1

	demultiplexer=Demultiplexer(sheet,inputDir,outputDir,lane)

if __name__=="__main__":
	doProfiling = False

	if doProfiling:
		profile.run('main()')
	else:
		main()


