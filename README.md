## Bulk RNA seq data analysis <br> Written by Hyojin Kim <br> June 2022~June 2023

<br>
<br>

For human data analysis, you can check my script through this link. <br> 
https://github.com/genehaus/bulkRNA-seq_2021/blob/main/code/run_process.new.R <br>
This pipeline used for this project. <br>

	```
	“Platelet-instructed SPP1(+) macrophages drive myofibroblast activation in fibrosis in a CXCL4-dependent manner”. Cell reports (2022) 
	“Dissecting CD8+ T cell pathology of severe SARS-CoV-2 infection by single-cell immunoprofiling.” Frontiers in immunology (2022)  
	```

<br>

##

This analysis is based on the data from mice. <br>

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
	curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/GRCm39.primary_assembly.genome.fa.gz 
	curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz
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
	
	
	3-2. to prepare the shell script to run nf-core/rnaseq in SLURM


		#!/usr/XXX/bin/zsh
		#SBATCH -J run_nf-core_rnaseq
		#SBATCH -t 100:00:00
		#SBATCH --ntasks-per-node=20
		#NXF_OPTS='-Xms1g -Xmx4g'
		#SBATCH --output=output.%J.txt
	
		source /xxx/anaconda3/bin/activate nf-core
	
		dir_gtf=/xxx/references/genecode_db_mice/gencode.vM27.annotation.gtf
		dir_fasta=/xxx/references/genecode_db_mice/GRCm39.primary_assembly.genome.fa
		dir_star=/xxx/references/genecode_db_index_mice
		dir_fq=/xxx/bulk_RNA/input
	
		cd $dir_fq
		cd $dir_fq
		/xxx/nextflow run nf-core/rnaseq --input $dir_fq/mydata..csv --gencode --star_index $dir_star --fasta $dir_fasta --gtf $dir_gtf --aligner star_salmon --outdir '/xxx/7_nextflow_out'


	(1) Here, 'xxx' should be your full directory 
	
	(2) how to make "mydata.csv" 

	go to : https://nf-co.re/rnaseq/usage#multiple-runs-of-the-same-sample





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
	cd /code/Remove_rRNA_tRNA_chrM/ # in this git directory 
	python Remove.rRNA.tRNA.MT.gft.python3.genecode.py gencode.vM28.annotation.gtf
	```

	this example is about Mouse data 	

	
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
	colnames(T) <- gsub("__", "_", colnames(T))
	condition <- gsub( paste0("_REP","[0-9]"),"", colnames(T))
	targets <- data.frame(sample=colnames(T), condition=condition)
	```
	
	
	4-7. run HTSFilter in order to remove lowly expressed genes 
	
	ref: https://github.com/andreamrau/HTSFilter
	
	if you have replicates for all conditions	
	
	```
	T <- HTSFilter(as.matrix(T), condition, s.min = 0.1, normalization ="TMM", s.len=25, plot=TRUE)
	T <- T$filteredData
	```
	
	If you don't have replicates, <br>
	try to remove lowely expressed genes <br>
	based on the 35 qualite of all read count. <br> 
	You can choose the quantile-based cutoff like 10%, 20%, etc. 
	 
	``` 
	cutoff <- quantile(as.matrix(T), probs = c(0.35))[[1]]
	print ("------------------------------")
	print ("print 35% quantile for cutoff")
	print (cutoff)
	print ("------------------------------")
	keep <- rowSums(as.matrix(T)) > cutoff
	T <- T[keep,]
	```
	
	
	4-8. Make a density plot after gene filtering <br>

	to check if the data distribution is more closely to normal distribution
	
	```
	pdf(paste0(outdir, "density.post.pdf"),width = 5, height = 5)
	d <- density(log2(as.matrix(T)+1)) # returns the density data
	print(plot(d)) # plots the result
	dev.off()
	```


	4-9. Add data to EdgeR object

	EdgeR : https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf <br>
	TMM : https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25 <br> 


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
	
	
	4-11. Make PCA plot with filtered count matrix data
	
	
	```
	gpca <- glmpca(T$counts, L=2)
	gpca.dat <- gpca$factors	
	gpca.dat$dex <- targets$condition
	gpca.dat$cell <- targets$condition
	
	pdf(paste0(outdir, "pca.post.counts.pdf"),width = 5, height = 5)
	ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell )) +
		geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
	dev.off()
	```


	
	OF-NOTE. If you want to jump to DESeq2, <br>
	please stop here and save outputs. <br>
	
	Go to : http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
	
	```	
	T$counts %>% write.table(file=paste0(outdir, "filtered.count.txt"),
					row.names = TRUE,
					col.names = TRUE,
					quote = FALSE,
					sep="\t")
	```
	
	
	
	4-12. Get DE gene by using EdgeR <br>

	what is glmGLFTest -> ref: https://rdrr.io/bioc/edgeR/man/glmQLFTest.html
	
	how to make a design matrix -> ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/

	Below is just an exmaple. <br>
	
	```
	design <- model.matrix(~0+condition, data=targets$samples)
	rownames(design) <- targets$sample
	```
	
	
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
		## https://rdrr.io/bioc/edgeR/man/glmQLFTest.html	

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




<br>
<br>




6. SessionInfo



```
> sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] stats4    parallel  grid      stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] EnhancedVolcano_1.12.0 matrixStats_0.61.0     dorothea_1.5.0        
 [4] progeny_1.14.0         erer_3.0               lmtest_0.9-39         
 [7] zoo_1.8-9              rlist_0.4.6.2          GO.db_3.13.0          
[10] AnnotationDbi_1.54.1   IRanges_2.26.0         S4Vectors_0.31.5      
[13] goseq_1.44.0           geneLenDataBase_1.28.0 BiasedUrn_1.07        
[16] data.table_1.14.2      vsn_3.60.0             Biobase_2.52.0        
[19] BiocGenerics_0.38.0    glmpca_0.2.0.9000      edgeR_3.34.1          
[22] limma_3.48.3           forcats_0.5.1          stringr_1.4.0         
[25] purrr_0.3.4            readr_2.1.1            tidyr_1.1.4           
[28] tibble_3.1.6           tidyverse_1.3.1        HTSFilter_1.32.0      
[31] textshape_1.7.3        hexbin_1.28.2          ggrepel_0.9.1         
[34] cowplot_1.1.1          dplyr_1.0.7            gridExtra_2.3         
[37] pheatmap_1.0.12        reshape_0.8.8          ggplot2_3.3.5
```




