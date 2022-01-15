# Written by H.J.Kim
# bulk RNA seq data analysis 

<br>
<br>

1. Install nf-core/rnaseq<br>

to map reads to reference genomes then make count matrix file

https://nf-co.re/rnaseq/3.5


```
curl -s https://get.nextflow.io | bash
conda create --name nf-core python=3.7 nf-core nextflow
conda activate nf-core
##(nf_core)$ nf-core download
```



2. Build genome index<br> 
http://www.regulatory-genomics.org/hint/introduction/

```
cd /xxx/references/genecode_db_mice
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/GRCm39.primary_assembly.genome.fa.gz -O GRCm39.primary_assembly.genome.fa.gz 
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz -O gencode.vM28.annotation.gtf.gz
mkdir genecode_db_index_mice
```

```
/xxx/STAR --runThreadN 20 \
		--runMode genomeGenerate \
		--genomeDir /xxx/references/genecode_db_index_mice/ \
		--genomeFastaFiles /xxx/references/genecode_db_mice/GRCm39.primary_assembly.genome.fa \
		--sjdbGTFfile /xxx/references/genecode_db_mice/gencode.vM28.annotation.gtf \
		--sjdbOverhang 100 > genecode_db_index_mice.log.txt
```

output in genecode_db_index_mice.log.txt

```
        STAR version: 2.7.9a   compiled: 2021-05-04T09:43:56-0400 vega:/home/dobin/data/STAR/STARcode/STAR.master/source
Jun 07 10:54:11 ..... started STAR run
Jun 07 10:54:11 ... starting to generate Genome files
Jun 07 10:54:56 ..... processing annotations GTF
Jun 07 10:55:19 ... starting to sort Suffix Array. This may take a long time...
Jun 07 10:55:34 ... sorting Suffix Array chunks and saving them to disk...
Jun 07 11:03:37 ... loading chunks from disk, packing SA...
Jun 07 11:05:09 ... finished generating suffix array
Jun 07 11:05:09 ... generating Suffix Array index
Jun 07 11:08:37 ... completed Suffix Array index
Jun 07 11:08:38 ..... inserting junctions into the genome indices
Jun 07 11:11:06 ... writing Genome to disk ...
Jun 07 11:11:09 ... writing Suffix Array to disk ...
Jun 07 11:11:30 ... writing SAindex to disk
Jun 07 11:11:32 ..... finished successfully
``` 





3. Run nf-core/atacseq (Python 3) 


	3-1. to make output dir

	
		mkdir /xxx/7_nextflow_out
	
	
	3-2. to prepare the shell script to run nf-core/atacseq


		#!/usr/XXX/bin/zsh
		#SBATCH -J run_nf-core_atac
		#SBATCH -t 100:00:00
		#SBATCH --ntasks-per-node=20
		#NXF_OPTS='-Xms1g -Xmx4g'
		#SBATCH --output=output.%J.txt
	
		source /xxx/anaconda3/bin/activate nf-core
	
		dir_fasta=/xxx/references/genecode_db/GRCh38.primary_assembly.genome.fa
		dir_star=/xxx/references/genecode_db_index
		dir_fq=/xxx/bulk_RNA/input
	
		cd $dir_fq
		/xxx/nextflow run nf-core/rnaseq --input $dir_fq/mydata.csv \
						 --gencode \
						 --star_index $dir_star \
						 --fasta $dir_fasta \
						 --gtf $dir_gtf \
						 --aligner star_salmon \
						--outdir '/xxx/7_nextflow_out'


	(1) Here, 'xxx' should be your full directory 
	
	(2) how to make "mydata.csv" 

	https://nf-co.re/rnaseq/usage#multiple-runs-of-the-same-sample





4. Run HINT


	4-1. to download the genome data 
	
	```
	# The directory "rgtdata" would be downloaded in your HOME when you install "RGT" (No.2)
	cd ~/rgtdata
	python setupGenomicData.py --mm10
	
	# you can check your genome name just by typing "python setupGenomicData.py"
	```





	4-2. To change chromosome name (optional)
	
	https://www.biostars.org/p/13462/
		

	```
	#!/usr/XXX/bin/zsh
	#SBATCH -J samtools
	#SBATCH -t 100:00:00
	#SBATCH --output=output.%J.txt
	
	source /xxx/anaconda3/bin/activate nf-core
	
	bam_dir="/xxx/7_nextflow_out/bwa/mergedLibrary"
	
	for file in $bam_dir/*sorted.bam
	do
	  filename=`echo $file | cut -d "." -f 1`
	  samtools view -H $file | \
	      sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | \
	      sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | \
	      sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | \
	      sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | \
	      sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \
	      sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | \
	      sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | \
	      sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | \
	      sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | \
	      sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | \
	      sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | \
	      sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
	      sed -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > ${filename}_chr.bam
	done
	
	
	``` 




	4-3. To make index file for the ${filename}_chr.bam (optional, followed by step 4-2)
	
	in Python3
	
	```
	import os,sys
	import glob
	import pysam
		
	bam_dir="/xxx/7_nextflow_out/bwa/mergedLibrary"
	files = glob.glob(bam_dir+"/*_chr.bam")
	for a_file in files :
	        pysam.index(a_file)
	
	```






	4-4. To run HINT
	
	https://www.regulatory-genomics.org/hint/tutorial/
	
	```
	#!/usr/local_rwth/bin/zsh
	#SBATCH -J HINT
	#SBATCH -t 100:00:00
	#SBATCH --output=output.%J.txt
	
	bam_dir="/xxx/7_nextflow_out/bwa/mergedLibrary"
	cd ${bam_dir}
	for file in *chr.bam
	do
        	filename=`echo $file | cut -d "_" -f 1`
        	/xxx/.local/bin/rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-location=$dir1/Footprints --output-prefix=${filename} ${filename}_R1.bam $peak_dir/${filename}_R1.mLb.clN_peaks.narrowPeak
	done
	```	
	
	
	
	


5. acknowledgement

5-1.team

http://www.costalab.org<br>
https://www.medizin.rwth-aachen.de/cms/Medizin/Die-Fakultaet/Einrichtungen/~dgun/IZKF-Aachen/<br>

5-2. website

https://nf-co.re/atacseq/1.2.1<br>
http://www.regulatory-genomics.org/hint/introduction/<br>
https://www.biostars.org/p/13462/<br>
https://www.regulatory-genomics.org/hint/tutorial/<br>
	


	
