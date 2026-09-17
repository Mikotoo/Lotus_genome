# -*- coding:utf-8 -*-
from collections import Counter

scientific_name=open("scientific_name.txt","r")
final_result =open("final_result.txt","w")

name_list_all=[]
for lines in scientific_name:
	line = lines.strip().split("\t")
	name = line[-1]
	name_list_all.append(name)

count_result = Counter(name_list_all)
count_list = list(count_result.items())
count_list.sort(key=lambda x:x[1],reverse=True)

final_result.write("Name\tHit_reads\tpercentage1\tpercentage2\n")
for i in count_list:
	name = i[0]
	number = i[1]
	reads_num = 100000
	percentage1 = "%.2f%%"%(100*float(number)/float(reads_num))
	percentage2 ="%.2f%%"%(100*float(number)/float(len(name_list_all)))
	final_result.write(name+"\t"+str(number)+"\t"+str(percentage1)+"\t"+str(percentage2)+"\n")


