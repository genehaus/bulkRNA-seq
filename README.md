# Written by Hyojin Kim
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
	'%ni%' = Negate('%in%')
	```


	4-2. Read file 

	
	``` 
	indir <- "/xxx/7_nextflow_output/star_salmon/"
	outdir <- "/xxx/out/"
	if(!dir.exists(outdir)) dir.create(outdir);
	
	D <- "salmon.merged.gene_counts.tsv"
	counts_df <- as.data.frame(read_delim(paste0(indir, D),
                                "\t",
                                escape_double = FALSE,
                                trim_ws = TRUE))
	```
	
	
	4-3. Filter out rRNA, mtRNA, non-expresed genes across samples 

	
	```
	clear_G_db <- as.data.frame(read_delim("/xxx/references/genecode_db_mice/gencode.vM28.annotation.gtf.no.rRNA.tRNA.gtf",
                                'gene_name "',
                                col_names = FALSE,
                                escape_double = FALSE,
                                trim_ws = TRUE,
                                ))
	
	print ("it'll take some time")
	clear_G <- strsplit(clear_G_db$X2, '["\"]' ) %>% lapply(., function(x) x[1]) %>% do.call("rbind", .) %>% unique()
	rm(clear_G_db)
	counts_df <- counts_df %>%
                dplyr::select(-gene_id) %>%
                filter(gene_name %in% clear_G)


	counts_df_orig <- counts_df
	counts_df_colnames <- colnames(counts_df[2:length(colnames(counts_df))])
	
	counts_df <- counts_df %>%
        	        rowwise() %>%
        	        mutate(sumVar = sum(eval(parse(text=counts_df_colnames)))) %>%
        	        filter(sumVar > 0) %>%
        	        dplyr::select(-sumVar) %>%
        	        as.data.frame()
	```
	
	
	4-4. Get maximum gene exp for the identical gene id 

	
	```
	order_G <- counts_df$gene_name %>% unique()
	max_col <- length(colnames(counts_df))
	columns <- colnames(counts_df)[2:max_col]
	counts_df_max <- c()
	for ( column in columns ) {
        	counts_df_max[[column]] <- counts_df %>%
        	                        dplyr::select(gene_name, column) %>%
        	                        group_by(gene_name) %>%
        	                        mutate(max = max(eval(parse(text=column)))) %>%
        	                        dplyr::select(-column) %>%
        	                        unique() %>%
        	                        as.data.frame() %>%
        	                        mutate(gene_name = fct_reorder(gene_name, order_G))

        	rownames(counts_df_max[[column]]) <- counts_df_max[[column]]$gene_name
        	counts_df_max[[column]] <- counts_df_max[[column]] %>% dplyr::select(-gene_name)
        	colnames(counts_df_max[[column]]) <- column
	
	        }
	T <- do.call("cbind", counts_df_max)

	```


	4-5. Make density plot before normalization 
	

	```
	pdf(paste0(outdir, "density.pre.pdf"),width = 5, height = 5)
	d <- density(log2(as.matrix(T)+1)) # returns the density data 
	print(plot(d)) # plots the result
	dev.off()
	```

	4-6. Prepare meta of samples

		
	```
	condition <- gsub( paste0("_","[0-9]"),"", colnames(T))
	targets <- data.frame(sample=colnames(T), condition=condition)
	design <- model.matrix(~0+condition, data=targets$samples)
	rownames(design) <- targets$sample
	```
	
	
	4-7. run HTSFilter in order to remove lowly expressed genes 

	https://github.com/andreamrau/HTSFilter
	
	
	```
	T <- HTSFilter(as.matrix(T), condition, s.min = 0.1, normalization ="TMM", s.len=25, plot=TRUE)
	T <- T$filteredData	
	```
	
	
	4-8. Make a density plot after gene filtering 

	to check if the data distribution is more closely to normal distribution
	
	
	```
	pdf(paste0(outdir, "density.post.pdf"),width = 5, height = 5)
	d <- density(log2(as.matrix(T)+1)) # returns the density data
	print(plot(d)) # plots the result
	dev.off()
	```


	4-9. Add data to EdgeR object

	```
	T <- DGEList(T, group=condition)
	T <- calcNormFactors(T, method="TMM")
	
	plotMD(T, column=1)
	abline(h=0, col="red", lty=2, lwd=2)
	```


	4-10. Make a box plot

	```
	pdf(paste0(outdir, "density.boxplot.without.left.with.right.TMM.pdf"),width = 5, height = 5)

		pre_TMM <- cpm(T, normalized.lib.sizes=FALSE)
		p1 <- log2(pre_TMM+1) %>%
        
		melt %>%
        		ggplot(aes(x=Var2, y=value)) +
        		geom_boxplot()

		post_TMM <- cpm(T, normalized.lib.sizes=TRUE)
		p2 <- log2(post_TMM+1) %>%
			        melt %>%
			        ggplot(aes(x=Var2, y=value)) +
	        		geom_boxplot()
		print(plot_grid(p1, p2))
	dev.off()
	
	write.table(post_TMM, file=paste0(outdir, "post_TMM.txt"),row.names = TRUE, quote = FALSE, sep="\t")	
	```
	
	
	4-11. Make PCA plot with post TMM-treated data
	
	
	```
	gpca <- glmpca(post_TMM, L=2)
	gpca.dat <- gpca$factors
	gpca.dat$dex <- targets$condition
	gpca.dat$cell <- targets$sample
	
	pdf(paste0(outdir, "pca.post.TMM.pdf"),width = 5, height = 5)
		ggplot(gpca.dat, aes(x = dim1, y = dim2, color = factor(dex))) +
			  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
	dev.off()
	```
	
	
	4-12. Get DE gene

	
	```
	plotMDS(T) # MDA plot like PCA 
	T <- estimateDisp(T, design, robust=TRUE)
	plotBCV(T)
	
	fit <- glmFit(T, design)
	Qfit <-  glmQLFit(T, design)
	my.contrasts <- makeContrasts(KO.vs.sham = conditionXKO-conditionXsham,
                                	levels=design)
	contrasts <- c("KO.vs.sham")


	## likelihood ratio test. -> glmFIT
	## quasi-likelihood (QL) F-test  -> glmGLFIT
	
	
	for (contrast in contrasts) {

        	lrt <- glmLRT(fit, contrast=my.contrasts[,contrast])
        	Qlrt <- glmQLFTest(Qfit, contrast=my.contrasts[,contrast])

	
        	lrt_fdr <- topTags(lrt, n=dim(post_TMM)[1])
        	Qlrt_fdr <- topTags(Qlrt, n=dim(post_TMM)[1])
        	

		## I chose to use "quasi-likelihood (QL) F-test" to utilize F value 
		## F value can be used to measure Z score with Log2FC later 
	

        	Qlrt_fdr$table %>% write.table(file=paste0(outdir, "top_KOvsWT.FDR.edgeR.QLFtest.", contrast, ".txt"),
        	                        row.names = TRUE,
        	                        quote = FALSE,
        	                        sep="\t")

	
	        Qlrt$df.total %>% write.table(file=paste0(outdir, "top_KOvsWT.FDR.edgeR.df.total.QLFtest.", contrast, ".txt"),
	                                row.names = TRUE,
	                                quote = FALSE,
	                                sep="\t")

	
	        pdf(paste0(outdir, "pca.post.counts.", contrast, ".pdf"),width = 5, height = 5)
	                summary(decideTests(lrt, p.value=0.05))
	                print(plotMD(lrt))
	                abline(h=c(-1, 1), col="blue")
	        dev.off()
	
	        pdf(paste0(outdir, "pca.post.counts.QLFtest.", contrast, ".pdf"),width = 5, height = 5)
	                summary(decideTests(Qlrt, p.value=0.05))
	                print(plotMD(Qlrt))
	                abline(h=c(-1, 1), col="blue")
	        dev.off()

        	}
	```
	
	
	4-13. Make a volcano plot

	https://www.bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
	
	```
	filenames <- Sys.glob( paste0(outdir, "top_KOvsWT.FDR.edgeR.QLFtest*txt"))
	
	for ( F in filenames) {

        	F <- strsplit(F, outdir) %>% unlist() %>% .[2]
        	DEG_o <- read.csv(paste0(outdir, F), sep="\t")
	
	        pdf(paste0(outdir, F, ".vocano.Log2FC.0.5.pdf"),width = 15, height = 10)
	        p <- EnhancedVolcano(DEG_o,
	                lab = rownames(DEG_o),
	                x = 'logFC',
	                y = 'FDR',
	                ylab = bquote('-' ~Log[10]~ 'FDR'),
	                pCutoff = 0.05,
	                FCcutoff = 0.5, ## 2.0
	                pointSize = 4.0,
	                labSize = 6.0,
	                labCol = 'black',
	                labFace = 'bold',
	                boxedLabels = TRUE,
	                colAlpha = 4/5,
	                    legendLabels=c('Not sig.',
	                                'Log (base 2) FC',
	                                'FDR',
	                                'FDR & Log (base 2) FC'),
	                legendPosition = 'right',
	                legendLabSize = 14,
	                legendIconSize = 4.0,
	                drawConnectors = TRUE,
	                widthConnectors = 1.0,
	                colConnectors = 'black')
	
	        print(p)
	        dev.off()
	        }
	```

	


	
	
	
	



5. acknowledgement

	5-1.team

	https://www.medizin.rwth-aachen.de/cms/Medizin/Die-Fakultaet/Einrichtungen/~dgun/IZKF-Aachen/<br>
	
	5-2. website
	
	https://davetang.org/muse/2011/01/24/normalisation-methods-for-dge-data/
	
	https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

	https://www.bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

	https://github.com/saezlab/transcriptutorial/blob/master/scripts/04_TranscriptionFactor_activity_with_Dorothea.md








	
