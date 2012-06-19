#!/usr/bin/env python

import sys

file=sys.argv[1]

objects=[]

for line in open(file):
	tokens=line.split(",")

	index=tokens[5-1].strip()

	if index=="Index":
		continue

	objects.append(index)

i=0

while i<len(objects):

	j=0

	while j<len(objects):

		if i<j:
			matches1=0
			k=0
			width=8

			while k<width:
				if objects[i][k:k+1]==objects[j][k:k+1]:
					if objects[i][k:k+1]!='-':
						matches1+=1
				k+=1

			k=8
			matches2=0
			while k<2*width:
				if objects[i][k:k+1]==objects[j][k:k+1]:
					if objects[i][k:k+1]!='-':
						matches2+=1
				k+=1

			print objects[i]+"	"+objects[j]+"	Matches1="+str(matches1)+"	Matches2="+str(matches2),
			if matches1>=width/2 and matches2>=width/2:
				print "Warning: collisions",

			print ""
		j+=1
	i+=1

