# -*- coding:utf-8 -*-

tiqu_gi = open("tiqu_gi.txt","r")
gi2taxid = open("cutted_nucl_gb.accession2taxid","r")
taxid2name = open("names.dmp","r")
get_name = open("scientific_name.txt","w")
from collections import defaultdict

taxid_name_dict={}
for lines in taxid2name:
	if "scientific name" in lines:
		line = lines.strip().split("|")
		taxid = line[0].strip()
		name = line[1].strip()
		taxid_name_dict[taxid]=name
		
tiqu_dict=defaultdict(list)
for lines in tiqu_gi:
	line = lines.strip().split("\t")
	gi = line[1].split("|")[1]
	tiqu_dict[gi].append("\t".join(line))


gi_taxid_dict={}
for lines in gi2taxid:
	line = lines.strip().split("\t")
	GI = line[1]
	taxid = line[0]
	gi_taxid_dict[GI]=taxid

jiaoji=set(tiqu_dict.keys())&set(gi_taxid_dict.keys())

tax_list=taxid_name_dict.keys()

tiqu_gi = open("tiqu_gi.txt","r")
for lines in tiqu_gi:
	line = lines.strip().split("\t")
	gi = line[1].split("|")[1]
	if gi in jiaoji:
		taxid=gi_taxid_dict[gi]
		if taxid in tax_list:
			get_name.write("\t".join(line)+"\t"+taxid_name_dict[taxid]+"\n")




