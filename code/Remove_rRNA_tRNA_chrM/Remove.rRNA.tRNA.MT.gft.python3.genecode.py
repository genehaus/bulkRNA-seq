## Written by H.J.Kim
## Date: May 31 2021 
##
## To Remove rRNA, tRNA, MT RNA on gencode.v38.annotation.gtf 
## [Step1. to download the gtf file
## curl -O http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
## gzip -d *.gz
## [Step2. to remove rRNA,tRNA and MT RNA]
## [Run & Python 3 & screen]$python Remove.rRNA.tRNA.MT.gft.python3.genecode.py gencode.v38.annotation.gtf
## [Output & current dir] ./gencode.v38.annotation.gtf.no.rRNA.tRNA.gtf
## --------------------------------- ## 

import os,sys

def read(a):
	col1=list(map(lambda x:x.split('\n')[0], open(a, 'r')))
	col2=list(map(lambda x:x.split('\t'), col1))
	return (col2)

input_data=read(sys.argv[1])
input_data_noMT=[ sub_i for sub_i in input_data if not sub_i[0].strip()=='chrM' ]

filter_out_biotype=['Mt_tRNA','rRNA','Mt_rRNA', 'rRNA_pseudogene']

genename=[ '\t'.join(i)
		for i in input_data_noMT if len(i) > 1
		if 'gene_type' in ''.join(i) 
		if not ''.join(i).split('gene_type ')[1].split('"')[1].strip() in filter_out_biotype
		]

fw=open(sys.argv[1]+'.no.rRNA.tRNA.gtf', 'w')
for sub in genename:
	print (sub)
	fw.write(sub+'\n')
fw.close()

