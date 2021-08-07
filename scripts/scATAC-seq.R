setwd("~/Desktop/-/School/scottlab")
install.packages('Seurat')
devtools::install_github("satijalab/seurat-data")
install.packages('SeuratData')
install.packages('BiocManager')
BiocManager::install('Rsamtools')
BiocManager::install('biovizBase')
BiocManager::install('multtest')
BiocManager::install('limma')
BiocManager::install("EnsDb.Hsapiens.v75")
install.packages('metap')
devtools::install_github("hhoeflin/hdf5r")
BiocManager::install("rhdf5")
install.packages("vroom")


library(BiocManager)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(readr)
set.seed(1234)
library(rhdf5)
library(vroom)
library(Matrix)
library(biovizBase)

### update: 
# finished pipeline, but can't get it to work with our data 
# I was able to extract counts matrix
# no fragments file or metadata file (which is optional)


# file <- read.delim("GSM3034633_PreFrontalCortex_62216.5kbwindowmatrix.txt", sep = "\t", header = TRUE)


h5file <- h5ls("atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")

# read count matrix from 10x CellRanger hdf5 file 
# GSM3034633_PreFrontalCortex_62216.peakmatrix.txt.gz
counts
counts_file <- readline(prompt="What is the counts file? ")
counts_file <- "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5" 
#counts_file2 <- "GSM3034633_PreFrontalCortex_62216.peakmatrix.txt"
#data = read.delim(counts_file2)
#data
metadata_file <- readline(prompt="What is the metadata file? ")
metadata_file <- "atac_v1_pbmc_10k_singlecell.csv"
# https://timoast.github.io/sinto/basic_usage.html -- how to create fragments file 
fragment_file <- readline(prompt="What is the fragment file? ")
fragment_file <- "atac_v1_pbmc_10k_fragments.tsv.gz" 
library(utils)
fread(fragment_file)
read.csv(fragment_file, sep="\t")

project_name = "trial"
scATAC-seq <- function(project_name, counts_file, metadata_file, fragment_file){
  if(grepl("\\.h5", counts_file)){ # if it's an .h5
    counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5") 
    
    # counts in peak matrix 
  }
  else{
    counts2 <- vroom("GSM3034633_PreFrontalCortex_62216.peakmatrix.txt", delim = "\t", col_names=TRUE) # A tibble: 436,206 x 6,572
    # potentially read_table  
    #counts2
    #counts2$annot2 = str_extract(counts2$annot, "(chr[0-21]|chr[X]|chr[Y])")
    # calls <- calls %>% separate(hap_cn, c("Maj", "Min"), sep="([a|b])" )
    # counts2$annot2 <- counts2 %>% separate(annot, c(""))
    # head(counts2 %>% select(chr, start, end, annot, annot2))
    
    # change format to chr:start-end
    counts2 <- counts2 %>% mutate(
      annot2 = paste(chr, start, sep=":")
    )
    counts2 <- counts2 %>% mutate(
      annot2 = paste(annot2, end, sep="-")
    )
    
    library(tidyverse)
    # put annot2 first 
    counts2 <- counts2 %>% dplyr::select(annot2, everything())
    
    counts2
    
    counts2 <- counts2 %>% dplyr::select(-c(chr, start, end, annot))
    
    counts2
    
    ### try to make it a sparse matrix
    counts3 <- Reduce(cbind, lapply(counts2[,-1], Matrix, sparse = TRUE))
    # or 
    counts4 <- as(counts2, "sparseMatrix")
    # or 
    counts3 <- Matrix(counts2, sparse = TRUE)
    rownames(counts2)
    counts4 <- counts2[-1]
    Matrix(as.matrix(counts2), sparse = TRUE)
    # or 
    counts10 <- map(counts2, Matrix, sparse = TRUE) %>% reduce(cbind2)
    
    
  }
  
  dir.name = paste(project_name, "_scATAC-seq", sep = "")
  dir.create(dir.name)
  metadata <- read.csv(
    file = metadata_file, 
    header = TRUE, 
    row.names = 1
  )
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg19',
    fragments = fragment_file,
    min.cells = 10,
    min.features = 200
  )
  #   some values in the "start" column cannot be turned into numeric values
  
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  ) # only need counts
  
  # pbmc 
  # - atac-seq data stored using custom assay - Chromatin Assay 
  # pbmc[['peaks']] contains additional information 
  
  granges(pbmc) # see genomic ranges associted with each feature 
  
  ## add gene annotations to pbmc object for the human genome
  
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75) # HERE
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg19"
  
  # add the gene information to the object
  Annotation(pbmc) <- annotations
  
  # compute nucleosome signal score per cell 
  pbmc <- NucleosomeSignal(object = pbmc)
  # compute TSS enrichment score per cell
  pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
  # add blacklist ratio and fraction of reads in peaks
  pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
  pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
  
  # inspect TSS enrichment scores by grouping cells based on score and plotting accessibility over TSS sites
  name <- paste(dir.name, "/TSSenrichment.pdf", sep = "")
  pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
  pdf(name)
  print(TSSPlot(pbmc, group.by = 'high.tss') + NoLegend())
  dev.off()
  
  # look at fragment length periodicity for all cells, group by cells with high or low nucleosomal signal strength
  name <- paste(dir.name, "/counts_per_fragment_length.pdf", sep = "")
  pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  pdf(name)
  print(FragmentHistogram(object = pbmc, group.by = 'nucleosome_group'))
  dev.off()
  
  name <- paste(dir.name, "/VLNplots.pdf", sep = "")
  pdf(name)
  print(VlnPlot(
    object = pbmc,
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  ))
  dev.off()
  
  ## remove cells that are outliers for the QC metrics 
  pbmc <- subset(
    x = pbmc,
    subset = peak_region_fragments > 3000 &
      peak_region_fragments < 20000 &
      pct_reads_in_peaks > 15 &
      blacklist_ratio < 0.05 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )
  
  ## normalization and linear dimensional reduction 
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  pbmc <- RunSVD(pbmc)
  
  name <- paste(dir.name, "/Depth_Correlation.pdf", sep = "")
  print(DepthCor(pbmc))
  dev.off()
  
  ## non-linear dimension reduction and clustering 
  pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
  name <- paste(dir.name, "/DimPlot.pdf", sep = "")
  print(DimPlot(object = pbmc, label = TRUE) + NoLegend())
  dev.off()
  
  ## create gene acitvity matrix 
  
  
  
}
