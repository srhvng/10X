# load libraries
suppressWarnings(suppressMessages(library("dplyr", quietly = T)))
library(Seurat)
library(patchwork)
library(getopt)
library(purrr)

#### functions that will be used
merge.all <- function(x, y){
  merge(x,y, merge.data=T)
}

#### set up the arguments to get from docker call in Seurat_dockerwraper.R ####
spec <- matrix(c(
  'TenXdir'              , 'i', 1, "character", "10X directory that contains the Sample directories that will be overlayed in the tSNE plot. Required",
  'SampleIDs'            , 's', 1, "character", "Sample IDs separated by a space, example: 27301994A1 27301994A2. Required",
  'doublet'              , 'd', 1, "character", "Directory that contains each samples' directory of Doublet/Barcode output table from Scrublet. Required",
  'outdir'               , 'o', 2, "character", "Alternate output directory, default is current working directory",
  'help'                 , 'h', 0, "logical"  , "Help request"
),ncol=5,byrow=T)

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
}, simplify=FALSE)
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)

params <- options.args

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(params$help)){
  cat(getopt(spec, usage =T))
  q(status=1)
}
if(is.null(params$TenXdir)){
  cat('\n--TenXdir and --SampleIDs are both required, use the --help flag for more information.\n')
  q(status=1)
}
if(is.null(params$SampleIDs)){
  cat('\n--TenXdir and --SampleIDs are both required, use the --help flag for more information.\n')
  q(status=1)
}
if(is.null(params$doublet)){
  cat('\n--doublet is required, use the --help flag for more information.\n')
  q(status=1)
}
if(is.null(params$outdir)){
  params$outdir = getwd()
} else {
  if(substr(x = params$outdir,start = nchar(params$outdir),stop = nchar(params$outdir))=="/"){
    params$outdir = substr(x = params$outdir,start = 1,stop = nchar(params$outdir)-1)
  }
}

# Read in 10x matrix data using lapply to get all the samples in the TenXdir directory
tenxdata <- lapply(params$SampleIDs, function(x) Read10X(data.dir=paste0(params$TenXdir,"/",x,"/filtered_feature_bc_matrix"))) 
names(tenxdata) <- params$SampleIDs
# Read using the h5 file
#tenxdata <- lapply(params$SampleIDs, function(x) Read10X_h5(filename = paste0(params$TenXdir,"/",x,"/outs/filtered_feature_bc_matrix.h5"),
#                                                            use.names = TRUE, unique.features = TRUE)) 

# Create Seurat objects
tenxdata_seurat <- lapply(seq_along(tenxdata), function(x) CreateSeuratObject(counts=tenxdata[[x]], project = names(tenxdata)[[x]]))

# Normalize and find variable features
tenxdata_seurat <- lapply(tenxdata_seurat, function(x){
  NormalizeData(x, verbose = F)
  FindVariableFeatures(x, selection.method = "vst", nfeatures=2000, verbose=F)
})
names(tenxdata_seurat) <- params$SampleIDs

# Create reference and anchors
reference.list <- tenxdata_seurat[params$SampleIDs]
my.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

# Make integrated expression matrix
mydata.integrated <- IntegrateData(anchorset = my.anchors, dims = 1:30)
# Scale data and run PCA
mydata.integrated <- ScaleData(mydata.integrated, verbose = FALSE)
mydata.integrated <- RunPCA(mydata.integrated, npcs = 30, verbose = FALSE)
# Create UMAP plots
mydata.integrated <- RunUMAP(mydata.integrated, reduction = "pca", dims = 1:30)

paste0(params$SampleIDs, collapse = "_")

# Save plot
png(file.path(params$outdir, paste0(paste0(params$SampleIDs, collapse = "_"),"_overlay_integrated_UMAP.png")), width=1000, height=800)
DimPlot(mydata.integrated, reduction = "umap", group.by = "orig.ident")
dev.off()

###################################################################
# Metrics for QC table
###################################################################
# Add/extract qc metrics per UMI for each cell to metadata
new_dat<-lapply(tenxdata_seurat, function(x) {
  x[["log10GenesPerUMI"]] <- (log10(x[["nFeature_RNA"]]) / log10(x[["nCount_RNA"]]))
  x[["nGene"]] <- x[["nFeature_RNA"]]
  x[["nUMI"]] <- x[["nCount_RNA"]]
  x[["seq_folder"]] <- x[["orig.ident"]]
  mitoRatio<- PercentageFeatureSet(object = x, pattern = "^MT-")
  x[["mitoRatio"]] <- mitoRatio / 100
  
  data.frame(x[["seq_folder"]], x[["log10GenesPerUMI"]], x[["nGene"]], x[["nUMI"]], x[["mitoRatio"]])
})

# Create metadata dataframe
metadata <- lapply(new_dat, function(x) {
  Barcode <- rownames(x)
  predicted_dead <- x %>% mutate(predicted_dead=ifelse(mitoRatio>0.1,"True","False"))
  meta.data<- cbind(Barcode, predicted_dead) %>% select(seq_folder, Barcode, log10GenesPerUMI, nGene, nUMI, mitoRatio, predicted_dead)
  
  data.frame(meta.data)
})

# Merge barcode from scurblet to table
doublet<- lapply(params$SampleIDs, function(x) {
  read.csv(paste0(params$doublet,"/",x,"/","scrublet_output_barcode_table.csv"))
})
# assign sample name to doublet tables for consistency
names(doublet)<- params$SampleIDs

qc.tab<-map2(metadata, doublet, left_join)

# save per barcode qc table
for(i in names(qc.tab)){
  write.csv(qc.tab[[i]], file.path(params$outdir,paste0(paste0(i, collapse = "_"),"_barcode_qc.csv")), row.names=FALSE)
}

new_dat<-lapply(qc.tab, function(x) {
  x[["sample"]]<- x[["seq_folder"]]
  
  deadT <- x %>% filter(predicted_dead=="True")
  x[["ratio.predicted_dead"]]<- round((nrow(deadT)/nrow(x)), digits=4) 
  doubletT <- x %>% filter(predicted_doublet=="True")
  x[["ratio.predicted_doublet"]]<- round((nrow(doubletT)/nrow(x)), digits = 4) 
  
  dead0.2<- x %>% filter(mitoRatio<0.2)
  x[["nMitoRatio<0.2"]]<- nrow(dead0.2)
  doublet0.25<- x %>% filter(doublet_score<0.25)
  x[["ndoublet<0.25"]]<- nrow(doublet0.25)
  
  data.frame(x[["sample"]],x[["ratio.predicted_dead"]], x[["ratio.predicted_doublet"]], x[["nMitoRatio<0.2"]], x[["ndoublet<0.25"]])[1,]
})

# save per sample qc table
qc<-do.call("rbind",new_dat)
colnames(qc)<-c("Sample","Ratio.predicted_dead","Ratio.predicted_doublet","nMitoRatio<0.2","ndoublet<0.25")
rownames(qc)<-NULL

write.csv(qc, file.path(params$outdir,"samples_qc.csv"), row.names=FALSE)

print("Successfully completed running Seurat_QC.R\n")
system('grep MemTotal /proc/meminfo')
