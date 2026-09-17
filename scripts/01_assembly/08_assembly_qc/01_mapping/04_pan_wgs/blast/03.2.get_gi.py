# -*- coding:utf-8 -*-
import re
from collections import defaultdict

xmlfile=open("MG20_subsample.blast.xml","r")
outfile=open("MG20_subsample_tiqu_gi.txt","w")

dict1=defaultdict(list)
for lines in xmlfile:
	line=lines.strip()
	read_id = re.match('<Iteration_query-def>.*</Iteration_query-def>',line)
	Hit_id = re.match('<Hit_id>.*</Hit_id>',line)
	Hit_def = re.match('<Hit_def>.*</Hit_def>',line)
	
	if read_id !=None:
		read_id=read_id.group()
		read_id = read_id.split("<")[1].split(">")[1]
		key=read_id

	elif Hit_id !=None:
		Hit_id = Hit_id.group()
		Hit_id = Hit_id.split("<")[1].split(">")[1]
		dict1[key].append(Hit_id)

	elif Hit_def !=None:
		Hit_def = Hit_def.group()
		Hit_def = Hit_def.split("<")[1].split(">")[1]
		dict1[key].append(Hit_def)

for key in dict1:
	outfile.write(key + "\t" + "\t".join(dict1[key])+"\n")


