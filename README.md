# Written by H.J.Kim
# bulk RNA seq data analysis 

<br>
<br>

1. Install nf-core/rnaseq<br>
https://nf-co.re/rnaseq/3.5

to map reads to reference genomes then make count matrix file

brief explanation <br> 
(1) input : fastq files  <br> 
(2) output : count matrix file <br>


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





3. Run nf-core/rnaseq (Python 3) 


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
		dir_star=/xxx/references/genecode_db_index_mice
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





4. Process


	strategy : counts (Normalize) -> cpm (Normalize) -> TMM (inter-sample normalize) -> EdgeR (DEG)

	4-1. call or install libraries

	```
	library(ggplot2)
	library(reshape)
	library(pheatmap)
	library(gridExtra)
	library(grid)
	library(dplyr)
	library(cowplot)
	library(ggrepel)
	library(hexbin)
	library(textshape)
	library(HTSFilter)
	library(tidyverse)
	library(tibble)
	library(edgeR)
	library(glmpca)
	library(readr)
	library(vsn)
	library(data.table)
	library(goseq)
	library(GO.db)
	library(rlist)
	library(erer)
	library(progeny)
	library(dorothea)
	library(pheatmap)
	library(readr)
	library(ggrepel)
	library(matrixStats)
	library(EnhancedVolcano)
	```













5. acknowledgement

5-1.team

https://www.medizin.rwth-aachen.de/cms/Medizin/Die-Fakultaet/Einrichtungen/~dgun/IZKF-Aachen/<br>

5-2. website

https://davetang.org/muse/2011/01/24/normalisation-methods-for-dge-data/

https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

	


	
