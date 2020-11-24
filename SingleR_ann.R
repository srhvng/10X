# Load libraries
library(Seurat)
library(SingleR)
library(Biobase)
library(org.Hs.eg.db)
library(dplyr)
library(data.table)
library(getopt)

#### functions that will be used
merge.all <- function(x, y){
  merge(x,y, merge.data=T)
}

#### set up the arguments to get from docker call in Seurat_dockerwraper.R ####
spec <- matrix(c(
  'TenXdir'              , 'i', 1, "character", "10X directory that contains the Sample directories that will be overlayed in the tSNE/UMAP plot. Required",
  'SampleIDs'            , 's', 1, "character", "Sample IDs separated by a space, example: 27301994A1 27301994A2. Required",
  'outdir'               , 'o', 2, "character", "alternate output directory, default is current working directory",
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
if(is.null(params$outdir)){
  params$outdir = getwd()
} else {
  if(substr(x = params$outdir,start = nchar(params$outdir),stop = nchar(params$outdir))=="/"){
    params$outdir = substr(x = params$outdir,start = 1,stop = nchar(params$outdir)-1)
  }
}

# Read in 10X matrix data
tenxdata <- lapply(params$SampleIDs, function(x) Read10X(data.dir=paste0(params$TenXdir,"/",x,"/filtered_feature_bc_matrix"))) 
#Read H5 file
#tenxdata <- lapply(params$SampleIDs, function(x) Read10X_h5(filename = paste0(params$TenXdir,"/",x,"/filtered_feature_bc_matrix.h5"),
#                                                            use.names = TRUE, unique.features = TRUE)) 

names(tenxdata) <- params$SampleIDs

# Create single Seurat object, normalize data and run PCA
mydata.combined <- CreateSeuratObject(counts = do.call(cbind,tenxdata), project = "PBMC")
mydata.combined <- NormalizeData(mydata.combined, verbose = FALSE)
mydata.combined <- FindVariableFeatures(mydata.combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
mydata.combined <- ScaleData(mydata.combined, verbose = FALSE)
mydata.combined <- RunPCA(mydata.combined, npcs = 20, verbose = FALSE)

# Define datasets and plot components
sample <- NULL
for(i in 1:length(tenxdata)){
  samp<-c(rep(paste(params$SampleIDs[i]),ncol(tenxdata[[i]])))
  samp<-data.frame(samp)
  sample<-rbind(sample,samp)
}
mydata.combined@meta.data$sample <- sample$samp

# Find neighbors, clusters and run UMAP
mydata.combined <- RunUMAP(mydata.combined, reduction = "pca", dims = 1:20)
mydata.combined <- FindNeighbors(mydata.combined, dims = 1:10)
mydata.combined <- FindClusters(mydata.combined, resolution = 0.5)

# Add cell type annotations
bp.se <- BlueprintEncodeData()
myresults <- SingleR(method = "single", test = mydata.combined@assays$RNA@counts, ref = list(BP=bp.se), labels = list(bp.se$label.fine), clusters= mydata.combined@active.ident)
mydata.combined <- FindClusters(mydata.combined, resolution = 0.8)

# Expand assignments and make cluster UMAP
cluster.results <- SingleR(method = "cluster", test = mydata.combined@assays$RNA@counts, ref = list(BP=bp.se), labels = list(bp.se$label.fine), clusters= mydata.combined@active.ident)
expand.assignments <- c()
for (i in mydata.combined@meta.data$seurat_clusters)
{
  ctypes <- cluster.results$labels[as.integer(i)+1]
  expand.assignments <- c(expand.assignments, ctypes)
}
mydata.combined@meta.data$singler.predictions.cluster <- expand.assignments

#extract cluster table
metadata<-data.frame(mydata.combined@meta.data)
clusters<-setDT(metadata)[, .(Freq = .N), by = .(sample, singler.predictions.cluster)]
clusters<-as.data.frame.matrix(xtabs(Freq ~ sample + singler.predictions.cluster, clusters))
clusters[,"Samples"]<-rownames(clusters)
clusters<- clusters %>% select(Samples,everything())

#Loop for generating percentage across genes
sample<-NULL
for (i in 2:length(clusters)){
  ratio<- round(clusters[,i]/rowSums(clusters[,2:length(clusters)]),3)
  sample<-data.frame(cbind(sample,ratio))
}
#rename percentage table
colnames(sample)<-paste0("ratio.",colnames(clusters[,2:length(clusters)]))

#Combine percent table to count table
neworder <- order(c(2*(seq_along(clusters) - 1) + 1,
                    2*seq_along(sample)))
clusters<-cbind(clusters,sample)[,neworder]

#save table
write.csv(clusters, (file.path(params$outdir, paste0("clusters.csv"))), row.names=FALSE)

#save plot
png(file.path(params$outdir, paste0("SingleR_UMAP.png")), width=1000, height=800)
DimPlot(mydata.combined, reduction = "umap", group.by = 'singler.predictions.cluster', split.by = "sample", pt.size=3, label = TRUE, repel = TRUE)
dev.off()

png(file.path(params$outdir, paste0("SingleR_VlnPC.png")), width=1000, height=800)
VlnPlot(object = mydata.combined, features = "PC_1", pt.size = .1, group.by = "sample")
dev.off()

print("Successfully completed running SingleR_ann.R\n")
system('grep MemTotal /proc/meminfo')
