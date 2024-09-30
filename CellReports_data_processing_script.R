library(Seurat) 
library(Matrix)
library(dplyr)
library(metap)
library(stringr)
library(ggplot2)
library(openxlsx)
library(reshape2)
library(ggrepel)
library(cowplot)

# Read in cell ranger outputs
#directories = list.dirs(path = "./EWAT", full.names = TRUE, recursive = TRUE)
directories = list.dirs(path = "./IWAT", full.names = TRUE, recursive = TRUE)
directories =directories[-1]

#iterate through the samples
for(sample in directories){
  #sample=directories[2]
  print(sample)
  splitSample = unlist(strsplit(sample,"/"))
  sampleName = splitSample[3]
  mid=unlist(strsplit(sampleName,"-"))
  #tissue="eWAT"
  tissue="iWAT"
  geneticbackground=mid[1]
  treatment=mid[2]
  group=str_sub(sampleName,1,-2)
  dropEST.data <- Read10X(data.dir = sample) 
  dropEST.data@Dimnames[[2]] = paste0(sampleName,"_",dropEST.data@Dimnames[[2]])
  dropEST.seurat <- CreateSeuratObject(counts = dropEST.data, project = sampleName)
  print(sampleName)
  rm(dropEST.data)
  dropEST.seurat@meta.data$data.geneticbackground = geneticbackground
  dropEST.seurat@meta.data$data.sampleName = sampleName
  dropEST.seurat@meta.data$data.tissue = tissue
  dropEST.seurat@meta.data$data.treatment = treatment
  dropEST.seurat@meta.data$data.group = group
  if(!exists("dropEST.combined")){
    dropEST.combined <- dropEST.seurat
    firstSampleName = sampleName
    firstSample = TRUE
  } else{
    if(firstSample==TRUE){
      dropEST.combined <- merge(x = dropEST.combined, y = dropEST.seurat, project = "Aging")
      firstSample = FALSE
    } else{
      dropEST.combined <- merge(x = dropEST.combined, y = dropEST.seurat, project = "Aging")
    }
  }
  rm(dropEST.seurat)
}

#saveRDS(dropEST.combined, file = "./eWATraw.rds")
saveRDS(dropEST.combined, file = "./iWATraw.rds")

#####QC
mito.features <- grep(pattern = "^mt-", x = rownames(x = dropEST.combined), value = TRUE) 
percent.mito <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[mito.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")) #mito.features/全体

# get ribosome %
ribo.features <- grep(pattern = "^Rps", x = rownames(x = dropEST.combined), value = TRUE)
ribo.features <- c(ribo.features, grep(pattern = "^Rpl", x = rownames(x = dropEST.combined), value = TRUE))
percent.ribo <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[ribo.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# get predicted genes % 
pred.features <- grep(pattern = "^Gm1", x = rownames(x = dropEST.combined), value = TRUE)
pred.features <- c(pred.features,grep(pattern = "^Gm2", x = rownames(x = dropEST.combined), value = TRUE)) #c添加内容
pred.features <- c(pred.features,grep(pattern = "^Gm3", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm4", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm5", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm6", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm7", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm8", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm9", x = rownames(x = dropEST.combined), value = TRUE))
percent.pred <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[pred.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# Add % mito, ribo and pred to the meta data
dropEST.combined[["percent.mito"]] <- percent.mito
dropEST.combined[["percent.ribo"]] <- percent.ribo
dropEST.combined[["percent.pred"]] <- percent.pred

table(dropEST.combined @meta.data$ data.sampleName)
ifelse(!dir.exists(file.path("./Plots/")), dir.create(file.path("./Plots/")), FALSE)
ifelse(!dir.exists(file.path("./Plots/","./QCPlots")), dir.create(file.path("./Plots/","./QCPlots")), FALSE)
ifelse(!dir.exists(file.path("./Plots/QCPlots","./ViolinPlots")), dir.create(file.path("./Plots/QCPlots","./ViolinPlots")), FALSE)
ifelse(!dir.exists(file.path("./Plots/QCPlots/ViolinPlots","./PreFilter")), dir.create(file.path("./Plots/QCPlots/ViolinPlots","./PreFilter")), FALSE)

source("./HelperFunctions.R")
cellQualityPlot(seuratObject=dropEST.combined,fileName="./Plots/QCPlots/ViolinPlots/PreFilter/SamplesQuality.pdf",H=9,W=40,
                featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo","percent.pred"),identPlot = "orig.ident",
                pointSize=0.5)
#Generate QCPlots
dropEST.combined = SetIdent(dropEST.combined, value = "data.sampleName")
pdf(file="./Plots/QCPlots/ViolinPlots/PreFilter/iWAT_PreFilter.pdf",height=9,width=40)
VlnPlot(object = dropEST.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.pred"),pt.size=0)
dev.off()

setThresholds(MinCells=3,MinGenes=200,MinUMIs=700,MinPercentMT=-Inf,MaxGenes=6000,
              MaxPercentMT=0.25,MaxUMIs=22000,MaxPercentRibo=1)
dropEST.combined.filtered = subset(x = dropEST.combined, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
                                   & percent.ribo < maxPercentRibo & nCount_RNA < maxUMIs)
table(dropEST.combined.filtered@meta.data$ data.sampleName)
ifelse(!dir.exists(file.path("./Plots/QCPlots/ViolinPlots","./PostFilter")), dir.create(file.path("./Plots/QCPlots/ViolinPlots","./PostFilter")), FALSE)

#Generate QCPlots
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "data.sampleName")
pdf(file="./Plots/QCPlots/ViolinPlots/PostFilter/iWAT_PostFilterchoice.pdf",height=9,width=40)
VlnPlot(object = dropEST.combined.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.pred"),pt.size = 0)
dev.off()

dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)
dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
length(dropEST.combined.filtered@assays$RNA@var.features) 
dropEST.combined.filtered <- ScaleData(dropEST.combined.filtered, verbose = FALSE)

#Perform PCA on the scaled data (uses the highly var genes)
dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)

#Plot PCA Plots
ifelse(!dir.exists(file.path("./Plots/QCPlots","./PCA")), dir.create(file.path("./Plots/QCPlots","./PCA")), FALSE)
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "condition")
pdf(file="./Plots/QCPlots/PCA/PCA_HighlyVarGenes.pdf")
DimPlot(object = dropEST.combined.filtered, dims = c(1,2), reduction = "pca")
dev.off()

#Plot Heatmaps of PCs
ifelse(!dir.exists(file.path("./Plots/QCPlots","./PCHeatmap")), dir.create(file.path("./Plots/QCPlots","./PCHeatmap75")), FALSE)
pdf(file="./Plots/QCPlots/PCHeatmap75/PCA_Heatmaps.pdf",height=20,width=20)
DimHeatmap(object = dropEST.combined.filtered, dims = 1:75, balanced = TRUE, cells = 100, reduction = "pca")
dev.off()

#Jackstraw permutation to determine the number of "Significant" PCs to use for tSNE 2D projection
dropEST.combined.filtered <- JackStraw(object = dropEST.combined.filtered, reduction = "pca", num.replicate = 50, 
                                       verbose = TRUE, dims = 50)
dropEST.combined.filtered = ScoreJackStraw(object = dropEST.combined.filtered, dims = 1:50, reduction = "pca")

#Visualize the Jackstraw permutations
ifelse(!dir.exists(file.path("./Plots/QCPlots","./JackStraw")), dir.create(file.path("./Plots/QCPlots","./JackStraw")), FALSE)

pdf(file="./Plots/QCPlots/JackStraw/ASC_JackStraw.pdf",height=45,width=10) #50 Sig PCs
JackStrawPlot(object = dropEST.combined.filtered, dims = 1:50)
dev.off()

dropEST.combined.filtered <- FindNeighbors(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50, k.param = 25)
dropEST.combined.filtered <- FindClusters(object = dropEST.combined.filtered, resolution = c(0.1,0.2,0.3,0.5,0.7,0.9,1.0,1.2,1.4),verbose = T, reduction = "pca")

dropEST.combined.filtered <- RunUMAP(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50)
dropEST.combined.filtered <- RunTSNE(object = dropEST.combined.filtered, reduction = "pca", dims = 1:50,check_duplicates = FALSE)
resolutions <- c("0.9","0.5","0.7","1.2","1.4")
resolutions <- c("0.1","0.2","0.3")
for(reso in resolutions){
  pdf(file=paste0("/Volumes/data/VlnPlot/umap、",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("RNA_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "umap"))
  dev.off()
}

for(reso in resolutions){
  pdf(file=paste0("./Plots/resoTSNE/iWATTSNEres",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("RNA_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "tsne"))
  dev.off()
}

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.sampleName")
p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.05, reduction = "umap", cols = c("red","pink","orange","yellow","cyan","lightblue"))+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("iWAT_sampleumap",".jpeg"), plot=p, device="jpeg",
       path="./Plots/", height=4,width=7, units="in", dpi=300)



dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.group")
my_levels <- c("WT-veh-iWAT", "obob-veh-iWAT","obob-rosi-iWAT")
dropEST.combined.filtered@meta.data$ data.group <- factor(Idents(dropEST.combined.filtered), levels= my_levels)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.sampleName")
my_levels <- c("WT-veh-iWAT1", "WT-veh-iWAT2","obob-veh-iWAT1","obob-veh-iWAT2","obob-rosi-iWAT1","obob-rosi-iWAT2")
dropEST.combined.filtered@meta.data$ data.sampleName <- factor(Idents(dropEST.combined.filtered), levels= my_levels)

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.group")
my_levels <- c("WT-veh-eWAT", "obob-veh-eWAT","obob-rosi-eWAT")
dropEST.combined.filtered@meta.data$ data.group <- factor(Idents(dropEST.combined.filtered), levels= my_levels)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.sampleName")
my_levels <- c("WT-veh-eWAT1", "WT-veh-eWAT2","obob-veh-eWAT1","obob-veh-eWAT2","obob-rosi-eWAT1","obob-rosi-eWAT2")
dropEST.combined.filtered@meta.data$ data.sampleName <- factor(Idents(dropEST.combined.filtered), levels= my_levels)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.group")
p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.05, reduction = "umap")+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("eWAT_groupumap",".jpeg"), plot=p, device="jpeg",
       path="/Volumes/data/eWATplots", height=4,width=7, units="in", dpi=300)

p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.1, reduction = "umap", split.by = "data.sampleName",repel = T)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("eWAT_samplesplit_umap",".jpeg"), plot=p, device="jpeg",
       path="/Volumes/data/eWATplots", height=4,width=25, units="in", dpi=500)

p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.1, reduction = "umap", split.by = "data.group",repel = T)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("eWAT_groupsplit_celltype",".jpeg"), plot=p, device="jpeg",
       path="/Volumes/data/mirian/eWATplots", height=4,width=15, units="in", dpi=500)

p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.1, reduction = "tsne", split.by = "data.sampleName",repel = T)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("eWAT_samplesplit_tsne",".jpeg"), plot=p, device="jpeg",
       path="./Plots/", height=4,width=25, units="in", dpi=500)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "cell.type.general")
p=DimPlot(dropEST.combined.filtered,label = F,pt.size = 0.1, reduction = "umap", split.by = "data.group",repel = T)+
  theme(text=element_text(size=22,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("iWAT_celltype_group_umap",".jpeg"), plot=p, device="jpeg",
       path="/Volumes/data/mirian/iWATplots", height=4,width=15, units="in", dpi=500)


###Do integration among groups within tissue and between tissues to further harmonize/refine cell labels 
### --- 1. integration by group for eWAT and iWAT, respectively --- ###
eWAT_Seuratobject <- readRDS("./raw_data/new_eWAT_metadata-twotype_fromGL.rds")
iWAT_Seuratobject <- readRDS("./raw_data/new_iWAT_metadata-twotype_fromGL.rds")
## Integrate by group for eWAT
eWAT_Seuratobject.list <- SplitObject(eWAT_Seuratobject, split.by = "data.newgroup")
eWAT_Seuratobject.list <- lapply(X = eWAT_Seuratobject.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
eWAT_features <- SelectIntegrationFeatures(object.list = eWAT_Seuratobject.list)
eWAT_anchors <- FindIntegrationAnchors(object.list = eWAT_Seuratobject.list, anchor.features = eWAT_features)
saveRDS(eWAT_anchors,"intermediate_data/eWAT_inegrateByGroup_eWAT_anchors.rds")
eWAT_combined <- IntegrateData(anchorset = eWAT_anchors)
DefaultAssay(eWAT_combined) <- "integrated"
eWAT_combined <- ScaleData(eWAT_combined, verbose = FALSE)
eWAT_combined <- RunPCA(eWAT_combined, npcs = 50, verbose = FALSE)
eWAT_combined <- RunUMAP(eWAT_combined, reduction = "pca", dims = 1:25)
eWAT_combined <- FindNeighbors(eWAT_combined, reduction = "pca", dims = 1:25)
eWAT_combined <- FindClusters(eWAT_combined, resolution = 0.5)
eWAT_colors <- c("APC" = "#56C4CA",
                 "BC" = "#F8A51D",
                 "EC" = "#F16523",
                 "LAM" = "#779FD4",
                 "Neu" = "#FFE4B6",
                 "NKT" = "#0F8D44",
                 "NPVM" = "#9083BD",
                 "PLAM" = "#F8E406",
                 "PVM" = "#99CA3F",
                 "SMC" = "#F29CC1")
unique(eWAT_combined$data.newgroup)
eWAT_combined$data.newgroup <- factor(eWAT_combined$data.newgroup, levels = c("WT-Veh","ob/ob-Veh","ob/ob-Rosi"))
unique(eWAT_combined$new.celltype)
eWAT_combined$new.celltype <- factor(eWAT_combined$new.celltype,levels = c("APC","LAM","PVM","NPVM","PLAM","BC","Neu","NKT","EC","SMC"))
pdf("output_plot/eWAT_integrateByGroup_umapByGroup.pdf", width = 12, height = 5)
DimPlot(eWAT_combined,reduction = "umap", split.by = "data.newgroup", group.by = "new.celltype", cols = eWAT_colors)
dev.off()
pdf("output_plot/eWAT_integrateByGroup_umap.pdf", width = 5, height = 5)
DimPlot(eWAT_combined,reduction = "umap",  group.by = "new.celltype", cols = eWAT_colors, label = T, repel = T)
dev.off()
saveRDS(eWAT_combined,"intermediate_data/Seuratobject_eWAT_integratedByGroup.rds")

## Integrate by group for iWAT  
iWAT_Seuratobject.list <- SplitObject(iWAT_Seuratobject, split.by = "data.newgroup")
iWAT_Seuratobject.list <- lapply(X = iWAT_Seuratobject.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
iWAT_features <- SelectIntegrationFeatures(object.list = iWAT_Seuratobject.list)
iWAT_anchors <- FindIntegrationAnchors(object.list = iWAT_Seuratobject.list, anchor.features = iWAT_features)
saveRDS(iWAT_anchors,"intermediate_data/iWAT_inegrateByGroup_iWAT_anchors.rds")
iWAT_combined <- IntegrateData(anchorset = iWAT_anchors)
DefaultAssay(iWAT_combined) <- "integrated"
iWAT_combined <- ScaleData(iWAT_combined, verbose = FALSE)
iWAT_combined <- RunPCA(iWAT_combined, npcs = 50, verbose = FALSE)
iWAT_combined <- RunUMAP(iWAT_combined, reduction = "pca", dims = 1:25)
iWAT_combined <- FindNeighbors(iWAT_combined, reduction = "pca", dims = 1:25)
iWAT_combined <- FindClusters(iWAT_combined, resolution = 0.5)
iWAT_colors <- c("APC" = "#56C4CA",
                 "BC" = "#F8A51D",
                 "EC" = "#F16523",
                 "LAM" = "#779FD4",
                 "Neu" = "#FFE4B6",
                 "NKT" = "#0F8D44",
                 "NPVM" = "#9083BD",
                 "DC" = "#9DD08B",
                 "PVM" = "#99CA3F")
unique(iWAT_combined$data.newgroup)
iWAT_combined$data.newgroup <- factor(iWAT_combined$data.newgroup, levels = c("WT-Veh","ob/ob-Veh","ob/ob-Rosi"))
unique(iWAT_combined$new.celltype)
iWAT_combined$new.celltype <- factor(iWAT_combined$new.celltype,levels = c("APC","LAM","PVM","NPVM","DC","BC","Neu","NKT","EC"))
pdf("output_plot/iWAT_integrateByGroup_umapByGroup.pdf", width = 12, height = 5)
DimPlot(iWAT_combined,reduction = "umap", split.by = "data.newgroup", group.by = "new.celltype", cols = iWAT_colors)
dev.off()
pdf("output_plot/iWAT_integrateByGroup_umap.pdf", width = 5, height = 5)
DimPlot(iWAT_combined,reduction = "umap",  group.by = "new.celltype", cols = iWAT_colors, label = T, repel = T)
dev.off()
saveRDS(iWAT_combined,"intermediate_data/Seuratobject_iWAT_integratedByGroup.rds")

### --- 2. integration eWAT and iWAT  --- ###
eWAT_Seuratobject <- readRDS("./raw_data/new_eWAT_metadata-twotype_fromGL.rds")
iWAT_Seuratobject <- readRDS("./raw_data/new_iWAT_metadata-twotype_fromGL.rds")
eWAT_Seuratobject$tissue.celltype <- paste(eWAT_Seuratobject$data.tissue, eWAT_Seuratobject$new.celltype, sep = "_")
unique(eWAT_Seuratobject$tissue.celltype)
iWAT_Seuratobject$tissue.celltype <- paste(iWAT_Seuratobject$data.tissue, iWAT_Seuratobject$new.celltype, sep = "_")
unique(iWAT_Seuratobject$tissue.celltype)
eWATiWAT_Seuratobject.list <- list(eWAT_Seuratobject,iWAT_Seuratobject)
names(eWATiWAT_Seuratobject.list) <- c("eWAT_Seurat","iWAT_Seurat")
eWATiWAT_Seuratobject.list <- lapply(X = eWATiWAT_Seuratobject.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
eWATiWAT_features <- SelectIntegrationFeatures(object.list = eWATiWAT_Seuratobject.list)
eWATiWAT_anchors <- FindIntegrationAnchors(object.list = eWATiWAT_Seuratobject.list, anchor.features = eWATiWAT_features)
saveRDS(eWATiWAT_anchors,"intermediate_data/eWATiWAT_inegrateByTissue_eWATiWAT_anchors.rds")
eWATiWAT_combined <- IntegrateData(anchorset = eWATiWAT_anchors)
DefaultAssay(eWATiWAT_combined) <- "integrated"
eWATiWAT_combined <- ScaleData(eWATiWAT_combined, verbose = FALSE)
eWATiWAT_combined <- RunPCA(eWATiWAT_combined, npcs = 50, verbose = FALSE)
ElbowPlot(eWATiWAT_combined,ndims = 50)
eWATiWAT_combined <- RunUMAP(eWATiWAT_combined, reduction = "pca", dims = 1:50)
eWATiWAT_combined <- FindNeighbors(eWATiWAT_combined, reduction = "pca", dims = 1:50)
eWATiWAT_combined <- FindClusters(eWATiWAT_combined, resolution = 0.5)
eWATiWAT_colors <- c("APC" = "#56C4CA",
                     "BC" = "#F8A51D",
                     "EC" = "#F16523",
                     "LAM" = "#779FD4",
                     "Neu" = "#FFE4B6",
                     "NKT" = "#0F8D44",
                     "NPVM" = "#9083BD",
                     "PLAM" = "#F8E406",
                     "PVM" = "#99CA3F",
                     "SMC" = "#F29CC1",
                     "DC" = "#9DD08B")
unique(eWATiWAT_combined$data.newgroup)
eWATiWAT_combined$data.newgroup <- factor(eWATiWAT_combined$data.newgroup, levels = c("WT-Veh","ob/ob-Veh","ob/ob-Rosi"))
unique(eWATiWAT_combined$new.celltype)
eWATiWAT_combined$new.celltype <- factor(eWATiWAT_combined$new.celltype,levels = c("APC","LAM","PVM","NPVM","PLAM","BC","Neu","NKT","EC","SMC","DC"))
pdf("output_plot/eWATiWAT_integrateByTissue_umapByTissue.pdf", width = 10, height = 5)
DimPlot(eWATiWAT_combined,reduction = "umap", split.by = "data.tissue", group.by = "new.celltype", cols = eWATiWAT_colors, raster = F)
dev.off()
pdf("output_plot/eWATiWAT_integrateByTissue_umap.pdf", width = 5, height = 5)
DimPlot(eWATiWAT_combined,reduction = "umap",  group.by = "new.celltype", cols = eWATiWAT_colors, label = T, repel = T, raster = F)
dev.off()
saveRDS(eWATiWAT_combined,"intermediate_data/Seuratobject_eWATiWAT_integratedByTissue.rds")
eWATiWAT_combined$tissue.newgroup <- paste(eWATiWAT_combined$data.tissue, eWATiWAT_combined$data.newgroup, sep = "_")
unique(eWATiWAT_combined$tissue.newgroup)
eWATiWAT_combined$tissue.newgroup <- factor(eWATiWAT_combined$tissue.newgroup,levels = c("eWAT_WT-Veh","eWAT_ob/ob-Veh","eWAT_ob/ob-Rosi","iWAT_WT-Veh","iWAT_ob/ob-Veh","iWAT_ob/ob-Rosi"))
pdf("output_plot/eWATiWAT_integrateByTissue_umapbyGroup.pdf", width = 22, height = 13)
DimPlot(eWATiWAT_combined,reduction = "umap",  group.by = "new.celltype", split.by = "tissue.newgroup", cols = eWATiWAT_colors, label = T, repel = T, raster = F, ncol = 3)
dev.off()
saveRDS(eWATiWAT_combined,"intermediate_data/Seuratobject_eWATiWAT_integratedByTissue.rds")
eWATiWAT_combined <- FindClusters(eWATiWAT_combined, resolution = 1.5)
saveRDS(eWATiWAT_combined,"intermediate_data/Seuratobject_eWATiWAT_integratedByTissue.res.1.5.rds")

### --- 3. Refine cell labels in iWAT based on integration results --- ### 
iWAT_Seuratobject <- readRDS("./raw_data/new_iWAT_metadata-twotype_fromGL.rds")
eWATiWAT_combined <- readRDS("intermediate_data/Seuratobject_eWATiWAT_integratedByTissue.res.1.5.rds")
Idents(eWATiWAT_combined) <- "integrated_snn_res.1.5"
iWAT_LAMs_inClust17.44 <- subset(eWATiWAT_combined, idents = c("17","44"))
iWAT_LAMs_inClust17.44 <- subset(iWAT_LAMs_inClust17.44, subset = tissue.celltype == "iWAT_LAM")
DefaultAssay(iWAT_LAMs_inClust17.44) <- "RNA"
unique(iWAT_LAMs_inClust17.44$new.celltype)
unique(iWAT_LAMs_inClust17.44$data.tissue)
unique(iWAT_LAMs_inClust17.44$twotype)
VlnPlot(iWAT_LAMs_inClust17.44,"Kif11", group.by = "data.tissue")
iWAT_LAMs_inClust17.44_cellid <- colnames(iWAT_LAMs_inClust17.44)
length(iWAT_LAMs_inClust17.44_cellid) # 663
DimPlot(eWATiWAT_combined, reduction = "umap",  group.by = "integrated_snn_res.1.5", cells.highlight = iWAT_LAMs_inClust17.44_cellid, raster = F, sizes.highlight = 0.1)
iWAT_Seuratobject_updated <- iWAT_Seuratobject
## NEW CELL TYPE will be stored in `new.celltype_YZ` column of metadata
## NOTE: FROM NOW ONLY, ONLY REFER TO `new.celltype_YZ` and `tissue.new.celltype_YZ` COLUMNS for updated cell type info in iWAT. 
iWAT_Seuratobject_updated$new.celltype_YZ <- iWAT_Seuratobject_updated$new.celltype
head(iWAT_Seuratobject_updated@meta.data)
iWAT_Seuratobject_updated@meta.data[rownames(iWAT_Seuratobject_updated@meta.data) %in% iWAT_LAMs_inClust17.44_cellid,]$new.celltype_YZ <- "PLAM"
iWAT_colors <- c("APC" = "#56C4CA",
                 "BC" = "#F8A51D",
                 "EC" = "#F16523",
                 "LAM" = "#779FD4",
                 "Neu" = "#FFE4B6",
                 "NKT" = "#0F8D44",
                 "NPVM" = "#9083BD",
                 "DC" = "#9DD08B",
                 "PVM" = "#99CA3F",
                 "PLAM" = "#F8E406")
unique(iWAT_Seuratobject_updated$data.newgroup)
iWAT_Seuratobject_updated$data.newgroup <- factor(iWAT_Seuratobject_updated$data.newgroup, levels = c("WT-Veh","ob/ob-Veh","ob/ob-Rosi"))
unique(iWAT_Seuratobject_updated$new.celltype_YZ)
iWAT_Seuratobject_updated$new.celltype_YZ <- factor(iWAT_Seuratobject_updated$new.celltype_YZ,levels = c("APC","LAM","PVM","NPVM","PLAM","DC","BC","Neu","NKT","EC"))
iWAT_Seuratobject_updated$tissue.new.celltype_YZ <- paste(iWAT_Seuratobject_updated$data.tissue, iWAT_Seuratobject_updated$new.celltype_YZ, sep = "_")
saveRDS(iWAT_Seuratobject_updated,"~/project-xyang123/Mirian_RosiObob_scRNAseq/ct.updated_scRNAseq_reanaysis/raw_data/iWAT_Seuratobject_ct.updated.rds")
pdf("~/project-xyang123/Mirian_RosiObob_scRNAseq/ct.updated_scRNAseq_reanaysis/output_plot/iWAT_ct.updated_umap.pdf", width = 5, height = 5)
DimPlot(iWAT_Seuratobject_updated,reduction = "umap",  group.by = "new.celltype_YZ", cols = iWAT_colors, label = T, repel = T)
dev.off()

### --- 4. Refine cell labels in eWAT based on integration results --- ### 
eWAT_Seuratobject <- readRDS("./raw_data/new_eWAT_metadata-twotype_fromGL.rds")
eWATiWAT_combined <- readRDS("intermediate_data/Seuratobject_eWATiWAT_integratedByTissue.res.1.5.rds")
Idents(eWATiWAT_combined) <- "integrated_snn_res.1.5"
eWAT_NPVMs_inClust25 <- subset(eWATiWAT_combined, idents = c("25"))
eWAT_NPVMs_inClust25 <- subset(eWAT_NPVMs_inClust25, subset = tissue.celltype == "eWAT_NPVM")
DefaultAssay(eWAT_NPVMs_inClust25) <- "RNA"
unique(eWAT_NPVMs_inClust25$new.celltype)
unique(eWAT_NPVMs_inClust25$data.tissue)
unique(eWAT_NPVMs_inClust25$twotype)
VlnPlot(eWAT_NPVMs_inClust25,"Flt3", group.by = "data.tissue") # confirmed expression
eWAT_NPVMs_inClust25_cellid <- colnames(eWAT_NPVMs_inClust25)
length(eWAT_NPVMs_inClust25_cellid) # 1117
DimPlot(eWATiWAT_combined, reduction = "umap",  group.by = "integrated_snn_res.1.5", cells.highlight = eWAT_NPVMs_inClust25_cellid, raster = F, sizes.highlight = 0.1)
eWAT_Seuratobject_updated <- eWAT_Seuratobject
## NEW CELL TYPE will be stored in `new.celltype_YZ` column of metadata
## NOTE: FROM NOW ONLY, ONLY REFER TO `new.celltype_YZ` and `tissue.new.celltype_YZ` COLUMNS for updated cell type info in iWAT. 
eWAT_Seuratobject_updated$new.celltype_YZ <- eWAT_Seuratobject_updated$new.celltype
head(eWAT_Seuratobject_updated@meta.data)
eWAT_Seuratobject_updated@meta.data[rownames(eWAT_Seuratobject_updated@meta.data) %in% eWAT_NPVMs_inClust25_cellid,]$new.celltype_YZ <- "DC"
eWAT_colors <- c("APC" = "#56C4CA",
                 "BC" = "#F8A51D",
                 "EC" = "#F16523",
                 "LAM" = "#779FD4",
                 "Neu" = "#FFE4B6",
                 "NKT" = "#0F8D44",
                 "DC" = "#9DD08B",
                 "NPVM" = "#9083BD",
                 "PLAM" = "#F8E406",
                 "PVM" = "#99CA3F",
                 "SMC" = "#F29CC1")
unique(eWAT_Seuratobject_updated$data.newgroup)
eWAT_Seuratobject_updated$data.newgroup <- factor(eWAT_Seuratobject_updated$data.newgroup, levels = c("WT-Veh","ob/ob-Veh","ob/ob-Rosi"))
unique(eWAT_Seuratobject_updated$new.celltype_YZ)
eWAT_Seuratobject_updated$new.celltype_YZ <- factor(eWAT_Seuratobject_updated$new.celltype_YZ,levels = c("APC","LAM","PVM","NPVM","PLAM","BC","Neu","NKT","EC","SMC","DC"))
eWAT_Seuratobject_updated$tissue.new.celltype_YZ <- paste(eWAT_Seuratobject_updated$data.tissue, eWAT_Seuratobject_updated$new.celltype_YZ, sep = "_")
unique(eWAT_Seuratobject_updated$tissue.new.celltype_YZ)
saveRDS(eWAT_Seuratobject_updated,"~/project-xyang123/Mirian_RosiObob_scRNAseq/ct.updated_scRNAseq_reanaysis/raw_data/eWAT_Seuratobject_ct.updated.rds")
pdf("~/project-xyang123/Mirian_RosiObob_scRNAseq/ct.updated_scRNAseq_reanaysis/output_plot/eWAT_ct.updated_umap.pdf", width = 5, height = 5)
DimPlot(eWAT_Seuratobject_updated,reduction = "umap",  group.by = "new.celltype_YZ", cols = eWAT_colors, label = T, repel = T)
dev.off()

###Redo UMAP for each tissue
eWAT_Seuratobject_ct.updated <- readRDS("./raw_data/eWAT_Seuratobject_ct.updated.rds")
eWAT_colors <- c("APC" = "#56C4CA",
                 "BC" = "#F8A51D",
                 "EC" = "#F16523",
                 "LAM" = "#779FD4",
                 "Neu" = "#FFE4B6",
                 "NKT" = "#0F8D44",
                 "DC" = "#9DD08B",
                 "NPVM" = "#9083BD",
                 "PLAM" = "#F8E406",
                 "PVM" = "#99CA3F",
                 "SMC" = "#F29CC1")
pdf("output_plot/eWAT_ct.updated_umap_splitbyGroup.pdf", width = 13, height = 5)
DimPlot(eWAT_Seuratobject_ct.updated,reduction = "umap",  group.by = "new.celltype_YZ", split.by = "data.newgroup",cols = eWAT_colors, label = F, repel = T)
dev.off()
pdf("./output_plot/eWAT_ct.updated_umap.pdf", width = 5, height = 5)
DimPlot(eWAT_Seuratobject_ct.updated,reduction = "umap",  group.by = "new.celltype_YZ", cols = eWAT_colors, label = T, repel = T)
dev.off()
iWAT_Seuratobject_ct.updated <- readRDS("./raw_data/iWAT_Seuratobject_ct.updated.rds")
iWAT_colors <- c("APC" = "#56C4CA",
                 "BC" = "#F8A51D",
                 "EC" = "#F16523",
                 "LAM" = "#779FD4",
                 "Neu" = "#FFE4B6",
                 "NKT" = "#0F8D44",
                 "NPVM" = "#9083BD",
                 "DC" = "#9DD08B",
                 "PVM" = "#99CA3F",
                 "PLAM" = "#F8E406")
pdf("output_plot/iWAT_ct.updated_umap_splitbyGroup.pdf", width = 13, height = 5)
DimPlot(iWAT_Seuratobject_ct.updated,reduction = "umap",  group.by = "new.celltype_YZ", split.by = "data.newgroup",cols = iWAT_colors, label = F, repel = T)
dev.off()
pdf("./output_plot/iWAT_ct.updated_umap.pdf", width = 5, height = 5)
DimPlot(iWAT_Seuratobject_ct.updated,reduction = "umap",  group.by = "new.celltype_YZ", cols = iWAT_colors, label = T, repel = T)
dev.off()
##VlnPlot for canonical markers
#e
cellTypeOrder = c("APC","BC","DC","EC","LAM","Neu","NKT","NPVM" ,"PLAM","PVM","SMC")
genes.to.use = c("Pdgfra","Cd79a","Flt3","Jam2","Trem2","S100a8","Cd3d","Ear2","Kif11","Lyve1","Msln","Adgre1")
my_colors=c("#56C4CA","#F8A51D","#9DD08B","#F16523","#779FD4","#FFE4B6","#0F8D44","#9083BD","#F8E406","#99CA3F","#F29CC1")
Idents(eWAT_Seuratobject_ct.updated) <- "new.celltype_YZ"
unique(Idents(eWAT_Seuratobject_ct.updated))
violinMatrix = NULL
for(celltype in cellTypeOrder){
  cellTypeSubset = subset(eWAT_Seuratobject_ct.updated, idents = celltype)
  for(gene in genes.to.use){
    expression = cellTypeSubset@assays$RNA@data[gene,]
    ExpAdd = cbind(rep(gene,length(expression)),expression,rep(celltype,length(expression)))
    if(is.null(violinMatrix)){
      violinMatrix = ExpAdd
    } 
    else{
      violinMatrix = rbind(violinMatrix,ExpAdd)
    }
  }
}
violinMatrix = as.data.frame(violinMatrix)
violinMatrix[,2]=as.numeric(violinMatrix[,2])  
colnames(violinMatrix) = c("Gene","Value","CellType")
rownames(violinMatrix) = seq(1,nrow(violinMatrix))
violinMatrix$Gene <- factor(violinMatrix$Gene, levels = genes.to.use)
violinMatrix$CellType <- factor(violinMatrix$CellType, levels = cellTypeOrder)
print(levels(violinMatrix$Gene))
print(levels(violinMatrix$CellType))
#lines between genes
line1 = data.frame(Gene=unique(violinMatrix[,'Gene']),Line=rep(-Inf,length(unique(violinMatrix[,'Gene']))))
line2 = data.frame(Gene=unique(violinMatrix[,'Gene']),Line=rep(Inf,length(unique(violinMatrix[,'Gene']))))
#lines around the outside
firstCellType = levels(violinMatrix$CellType)[1]
lastCellType = levels(violinMatrix$CellType)[length(levels(violinMatrix$CellType))]
line3 = data.frame(CellType=c(firstCellType,lastCellType),Line=c(Inf,-Inf))
pdf(file="./output_plot/eWAT_markers_vlnplot_GaoyanScript.pdf",height=10, width=13)
ggplot(violinMatrix, aes(x = CellType, y = Value, fill = CellType)) +
  facet_grid(CellType~Gene, switch = "y", scales = "free") +
  # theme(text=element_text(size=24,family="Arial"))+
  theme(text=element_text(size=24))+
  geom_violin(adjust=1) +
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  geom_hline(data = line1, aes(yintercept = Line)) +
  geom_hline(data = line2, aes(yintercept = Line)) +
  geom_vline(data = line3, aes(xintercept = Line)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 180, hjust = 1),
        strip.text.x = element_text(angle = 270, vjust = 0, face = "italic"), 
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        legend.position = "none"
  )
dev.off()
#i
cellTypeOrder = c("APC","BC","DC","EC","LAM","Neu","NKT","NPVM","PLAM","PVM")
genes.to.use = c("Pdgfra","Cd79a","Flt3","Jam2","Trem2","S100a8","Cd3e","Ear2","Kif11","Lyve1","Adgre1")
my_colors=c("#56C4CA","#F8A51D","#9DD08B","#F16523","#779FD4","#FFE4B6","#0F8D44","#9083BD","#F8E406","#99CA3F")
Idents(iWAT_Seuratobject_ct.updated) <- "new.celltype_YZ"
unique(Idents(iWAT_Seuratobject_ct.updated))
violinMatrix = NULL
for(celltype in cellTypeOrder){
  cellTypeSubset = subset(iWAT_Seuratobject_ct.updated, idents = celltype)
  for(gene in genes.to.use){
    expression = cellTypeSubset@assays$RNA@data[gene,]
    ExpAdd = cbind(rep(gene,length(expression)),expression,rep(celltype,length(expression)))
    if(is.null(violinMatrix)){
      violinMatrix = ExpAdd
    } 
    else{
      violinMatrix = rbind(violinMatrix,ExpAdd)
    }
  }
}
violinMatrix = as.data.frame(violinMatrix)
violinMatrix[,2]=as.numeric(violinMatrix[,2])  
#violinMatrix[,2] = unfactor(violinMatrix[,2])
colnames(violinMatrix) = c("Gene","Value","CellType")
rownames(violinMatrix) = seq(1,nrow(violinMatrix))
violinMatrix$Gene <- factor(violinMatrix$Gene, levels = genes.to.use)
violinMatrix$CellType <- factor(violinMatrix$CellType, levels = cellTypeOrder)
print(levels(violinMatrix$Gene))
print(levels(violinMatrix$CellType))
#lines between genes
line1 = data.frame(Gene=unique(violinMatrix[,'Gene']),Line=rep(-Inf,length(unique(violinMatrix[,'Gene']))))
line2 = data.frame(Gene=unique(violinMatrix[,'Gene']),Line=rep(Inf,length(unique(violinMatrix[,'Gene']))))
#lines around the outside
firstCellType = levels(violinMatrix$CellType)[1]
lastCellType = levels(violinMatrix$CellType)[length(levels(violinMatrix$CellType))]
line3 = data.frame(CellType=c(firstCellType,lastCellType),Line=c(Inf,-Inf))
pdf(file="./output_plot/iWAT_markers_vlnplot_GaoyanScript.pdf",height=10, width=13)
ggplot(violinMatrix, aes(x = CellType, y = Value, fill = CellType)) +
  facet_grid(CellType~Gene, switch = "y", scales = "free") +
  # theme(text=element_text(size=24,family="Arial"))+
  theme(text=element_text(size=24))+
  geom_violin(adjust=1) +
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  geom_hline(data = line1, aes(yintercept = Line)) +
  geom_hline(data = line2, aes(yintercept = Line)) +
  geom_vline(data = line3, aes(xintercept = Line)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 180, hjust = 1),
        strip.text.x = element_text(angle = 270, vjust = 0, face = "italic"), 
        strip.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        legend.position = "none"
  )
dev.off()

#test markers
genes=c("Lyz1","Cd14","Adgre1","Itgam","Itgax","Tnf","Il6","Mrc1","Ptgs1","Retnla","Siglech","Flt3","Cd24a","Csf2","Arg1","Ncr1","Klrb1c",
        "Cd3e","Klra1","Cd4","Cd8a","Foxp3","Cd163l1","Il18r1","Cd79a","Cd79b","Cd19","Ly6g","S100a8","S100a9","Siglecf","Il5ra")
p=DotPlot(eWAT,features = genes,cols = c("lightgrey","red"))+ RotatedAxis()+
  theme(text=element_text(size=20,family="Arial"),plot.title=element_text(size=24,face="italic"),
        axis.text.x=element_text(size=18),axis.text.y=element_text(size=18))
ggsave(filename=paste0("eWAT_Immune_dotplot",".jpeg"), plot=p, device="jpeg",
       path="./Plots/", height=7,width=15, units="in", dpi=300)

###Find DEG
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "new.celltype")
cellTypeNames = dropEST.combined.filtered@active.ident
cellTypeNames = as.vector(cellTypeNames)
names(cellTypeNames) = names(dropEST.combined.filtered@active.ident)
cellTypeNamesGeneral = cellTypeNames
cellTypeDEGsStandard = list()
for(cellType in levels(dropEST.combined.filtered@active.ident)){
  #for(cellType in celltypes){
  print(cellType)
  cellTypeSubset = subset(dropEST.combined.filtered,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = "data.newgroup")
  #StandardDEGs = FindMarkers(cellTypeSubset, ident.1 = "ob/ob-Veh", ident.2 = "WT-Veh")
  StandardDEGs = FindMarkers(cellTypeSubset, ident.1 = "ob/ob-Rosi", ident.2 = "ob/ob-Veh")
  cellTypeDEGsStandard[[cellType]] = StandardDEGs
}  
#save DEGs
saveRDS(cellTypeDEGsStandard, file = "/Volumes/data/iWATrosiDEGs.rds")
DEGs=cellTypeDEGsStandard
celltypes=names(DEGs)
all_DEGs = data.frame(stringsAsFactors = FALSE)
for(celltype in celltypes){
  print(celltype)
  temp = DEGs[[celltype]]
  temp2 <- temp
  GENE <- rownames(temp2)
  rownames(temp2) <- NULL
  temp3 <- cbind(GENE,temp2)
  temp3$`CellType` = celltype
  temp3$`Compare` = "Rosi_effect"
  temp3$`Tissue` = "iWAT"
  temp3$`Tissue_Compare` = "iWAT_Rosi_effect"
  all_DEGs = rbind(all_DEGs, temp3)
}
all_DEGs=all_DEGs[all_DEGs$p_val_adj < 0.05,]
up=all_DEGs[all_DEGs$avg_log2FC > log2(1.1),]
str(up)
down=all_DEGs[all_DEGs$avg_log2FC < -log2(1.1),]
str(down)
iROSIup<-up 
iROSIdown<-down
all_DEGs = rbind(iOBup,iOBdown,iROSIup,iROSIdown,eOBup,eOBdown,eROSIup,eROSIdown)
write.table(all_DEGs,file="/Volumes/data/eandiDEGall.txt",sep='\t',row.names=TRUE) 

####Make DEG split figure
removeQuotes <- function(x) gsub("\"", "", x)
library(ComplexHeatmap)
library(UpSetR)
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  dat<- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  dat[is.na(dat)] <- as.integer(0)
  dat[dat != 0] <- as.integer(1)
  dat <- data.frame(matrix(dat, ncol = length(input), byrow = F))
  dat <- dat[which(rowSums(dat) != 0), ]
  names(dat) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(dat) <- elements
  return(dat)
}

iOBup <- iOBup %>%mutate_if(is.character, removeQuotes)
iOBdown <- iOBdown %>%mutate_if(is.character, removeQuotes)
iROSIup <- iROSIup %>%mutate_if(is.character, removeQuotes)
iROSIdown <- iROSIdown %>%mutate_if(is.character, removeQuotes)
eOBup <- eOBup %>%mutate_if(is.character, removeQuotes)
eOBdown <- eOBdown %>%mutate_if(is.character, removeQuotes)
eROSIup <- eROSIup %>%mutate_if(is.character, removeQuotes)
eROSIdown <- eROSIdown %>%mutate_if(is.character, removeQuotes)

#for(celltype in celltypes){
setdiff(unique(eOBup$CellType),unique(iOBup$CellType))
setdiff(unique(iOBup$CellType),unique(eOBup$CellType))

celltypes=union(unique(eOBup$CellType),unique(iOBup$CellType))
celltype=celltypes[11]
print(celltype)
tempeOBup=eOBup[eOBup$CellType == celltype,]$GENE
tempeOBdown=eOBdown[eOBdown$CellType == celltype,]$GENE
tempeROSIup=eROSIup[eROSIup$CellType == celltype,]$GENE
tempeROSIdown=eROSIdown[eROSIdown$CellType == celltype,]$GENE

tempiOBup=iOBup[iOBup$CellType == celltype,]$GENE
tempiOBdown=iOBdown[iOBdown$CellType == celltype,]$GENE
tempiROSIup=iROSIup[iROSIup$CellType == celltype,]$GENE
tempiROSIdown=iROSIdown[iROSIdown$CellType == celltype,]$GENE

ListInput<-list("ob/ob_effect_UP"= tempiOBup,"Rosi_effect_DOWN"= tempiROSIdown,"ob/ob_effect_DOWN"=tempiOBdown, "Rosi_effect_UP"= tempiROSIup )
lt=fromList(ListInput)
m = make_comb_mat(lt)
pdf(paste0("/Volumes/data\ 1/mirian/newfiguretoshare/iWATDEG_",celltype,".pdf"), height = 3.5, width = 5)

UpSet(m,set_order=c("ob/ob_effect_UP","Rosi_effect_DOWN","ob/ob_effect_DOWN","Rosi_effect_UP"),pt_size = unit(5, "mm"),lwd = 3,
      comb_col = c("#C99f6e","#ffd28f")[comb_degree(m)])
dev.off()

ListInput<-list("ob/ob_effect_UP"= tempeOBup,"ob/ob_effect_DOWN"=tempeOBdown, "Rosi_effect_UP"= tempeROSIup,"Rosi_effect_DOWN"= tempeROSIdown )
lt=fromList(ListInput)
m = make_comb_mat(lt)
pdf(paste0("/Volumes/data\ 1/mirian/newfiguretoshare/eWATDEG_",celltype,".pdf"), height = 3.5, width = 5)

UpSet(m,set_order=c("ob/ob_effect_UP","Rosi_effect_DOWN","ob/ob_effect_DOWN","Rosi_effect_UP"),pt_size = unit(5, "mm"),lwd = 3,
      comb_col = c("#C99f6e","#ffd28f")[comb_degree(m)])
dev.off()

ListInput<-list("ob/ob_UP"= tempiOBup,"Rosi_DOWN"= tempiROSIdown,"ob/ob_DOWN"=tempiOBdown, "Rosi_UP"= tempiROSIup )
lt=fromList(ListInput)
m = make_comb_mat(lt)
pdf(paste0("/Volumes/data\ 1/mirian/newfiguretoshare/iWAT/iWATDEG_",celltype,".pdf"), height = 3.5, width = 5)

UpSet(m,set_order=c("ob/ob_UP","Rosi_DOWN","ob/ob_DOWN","Rosi_UP"),pt_size = unit(5, "mm"),lwd = 3,
      comb_col = c("#C99f6e","#ffd28f")[comb_degree(m)])
dev.off()

ListInput<-list("ob/ob_UP"= tempeOBup,"ob/ob_DOWN"=tempeOBdown, "Rosi_UP"= tempeROSIup,"Rosi_DOWN"= tempeROSIdown )
lt=fromList(ListInput)
m = make_comb_mat(lt)
pdf(paste0("/Volumes/data\ 1/mirian/newfiguretoshare/eWAT/eWATDEG_",celltype,".pdf"), height = 3.5, width = 5)

UpSet(m,set_order=c("ob/ob_UP","Rosi_DOWN","ob/ob_DOWN","Rosi_UP"),pt_size = unit(5, "mm"),lwd = 3,
      comb_col = c("#C99f6e","#ffd28f")[comb_degree(m)])
dev.off()


######################
#for(celltype in celltypes){
celltypes=union(unique(eOBup$`"CellType"`),unique(iOBup$`"CellType"`))
celltype=celltypes[2]
CT=celltype
print(celltype)

tempeOBup=eOBup[eOBup$`"CellType"` == celltype,]$`"GENE"`
tempeOBdown=eOBdown[eOBdown$`"CellType"` == celltype,]$`"GENE"`
tempeROSIup=eROSIup[eROSIup$`"CellType"` == celltype,]$`"GENE"`
tempeROSIdown=eROSIdown[eROSIdown$`"CellType"` == celltype,]$`"GENE"`

tempiOBup=iOBup[iOBup$`"CellType"` == celltype,]$`"GENE"`
tempiOBdown=iOBdown[iOBdown$`"CellType"` == celltype,]$`"GENE"`
tempiROSIup=iROSIup[iROSIup$`"CellType"` == celltype,]$`"GENE"`
tempiROSIdown=iROSIdown[iROSIdown$`"CellType"` == celltype,]$`"GENE"`

genese1=intersect(tempeOBup,tempeROSIdown)
genese2=intersect(tempeOBdown,tempeROSIup)
genesi1=intersect(tempiOBup,tempiROSIdown)
genesi2=intersect(tempiOBdown,tempiROSIup)
e1not=c(tempeOBdown,tempeROSIup,tempiOBup,tempiOBdown,tempiROSIup,tempiROSIdown)
agenese1=setdiff(genese1, e1not)
e2not=c(tempeOBup,tempeROSIdown,tempiOBup,tempiOBdown,tempiROSIup,tempiROSIdown)
agenese2=setdiff(genese2, e2not)
i1not=c(tempeOBup,tempeROSIdown,tempeOBdown,tempeROSIup,tempiOBdown,tempiROSIup)
agenesi1=setdiff(genesi1, i1not)
i2not=c(tempeOBup,tempeROSIdown,tempeOBdown,tempeROSIup,tempiOBup,tempiROSIdown)
agenesi2=setdiff(genesi2, i2not)
tempeOBup=eOBup[eOBup$`"CellType"` == celltype,]

####e
res=tempeOBup
genes=genese1
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="e_obUProsiDOWN"
e_obUProsiDOWN=selected

tempeOBdown=eOBdown[eOBdown$`"CellType"` == celltype,]
res=tempeOBdown
genes=genese2
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="e_obDOWNrosiUP"
e_obDOWNrosiUP=selected

##########i
tempiOBup=iOBup[iOBup$`"CellType"` == celltype,]
res=tempiOBup
genes=genesi1
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="i_obUProsiDOWN"
i_obUProsiDOWN=selected

tempiOBdown=iOBdown[iOBdown$`"CellType"` == celltype,]
res=tempiOBdown
genes=genesi2
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="i_obDOWNrosiUP"
i_obDOWNrosiUP=selected

combined = rbind(e_obUProsiDOWN,es_obUProsiDOWN,e_obDOWNrosiUP, es_obDOWNrosiUP,i_obUProsiDOWN,is_obUProsiDOWN,i_obDOWNrosiUP, is_obDOWNrosiUP)
combined$CT=CT
#write.csv(combined,"/Volumes/data/splittedDEG/preadipoctesplittedDEG.csv")

celltypes=union(unique(eOBup$`"CellType"`),unique(iOBup$`"CellType"`))
celltype=celltypes[11]
CT=celltype
print(celltype)

tempeOBup=eOBup[eOBup$`"CellType"` == celltype,]$`"GENE"`
tempeOBdown=eOBdown[eOBdown$`"CellType"` == celltype,]$`"GENE"`
tempeROSIup=eROSIup[eROSIup$`"CellType"` == celltype,]$`"GENE"`
tempeROSIdown=eROSIdown[eROSIdown$`"CellType"` == celltype,]$`"GENE"`

tempiOBup=iOBup[iOBup$`"CellType"` == celltype,]$`"GENE"`
tempiOBdown=iOBdown[iOBdown$`"CellType"` == celltype,]$`"GENE"`
tempiROSIup=iROSIup[iROSIup$`"CellType"` == celltype,]$`"GENE"`
tempiROSIdown=iROSIdown[iROSIdown$`"CellType"` == celltype,]$`"GENE"`

genese1=intersect(tempeOBup,tempeROSIdown)
genese2=intersect(tempeOBdown,tempeROSIup)
genesi1=intersect(tempiOBup,tempiROSIdown)
genesi2=intersect(tempiOBdown,tempiROSIup)
e1not=c(tempeOBdown,tempeROSIup,tempiOBup,tempiOBdown,tempiROSIup,tempiROSIdown)
agenese1=setdiff(genese1, e1not)
e2not=c(tempeOBup,tempeROSIdown,tempiOBup,tempiOBdown,tempiROSIup,tempiROSIdown)
agenese2=setdiff(genese2, e2not)
i1not=c(tempeOBup,tempeROSIdown,tempeOBdown,tempeROSIup,tempiOBdown,tempiROSIup)
agenesi1=setdiff(genesi1, i1not)
i2not=c(tempeOBup,tempeROSIdown,tempeOBdown,tempeROSIup,tempiOBup,tempiROSIdown)
agenesi2=setdiff(genesi2, i2not)

#这里读的是整列
tempeROSIdown=eROSIdown[eROSIdown$`"CellType"` == celltype,]

####e
res=tempeROSIdown
genes=genese1
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="e_obUProsiDOWN"
e_obUProsiDOWN=selected
tempeROSIup=eROSIup[eROSIup$`"CellType"` == celltype,]
res=tempeROSIup
genes=genese2
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="e_obDOWNrosiUP"
e_obDOWNrosiUP=selected

##########i
tempiROSIdown=iROSIdown[iROSIdown$`"CellType"` == celltype,]
res=tempiROSIdown
genes=genesi1
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="i_obUProsiDOWN"
i_obUProsiDOWN=selected

tempiROSIup=iROSIup[iROSIup$`"CellType"` == celltype,]
res=tempiROSIup
genes=genesi2
selected = data.frame(stringsAsFactors = FALSE)
for(gene in genes){
  temp = res[res$`"GENE"`==gene,]
  selected = rbind(selected, temp)
}
head(selected)
selected$"CellType"="i_obDOWNrosiUP"
i_obDOWNrosiUP=selected
combined = rbind(e_obUProsiDOWN,e_obDOWNrosiUP,i_obUProsiDOWN,i_obDOWNrosiUP)
combined$CT=CT
write.csv(combined,file=paste0("/Volumes/data/splittedDEG/",CT,"splittedDEG.csv")) 
  


###Ribo Dotplot
Fig5A_ribogenesinAPCs_fromGL <- read_csv("./raw_data/Fig5A_ribogenesinAPCs_fromGL.csv")
View(Fig5A_ribogenesinAPCs_fromGL)
Fig5A_ribogenesinAPCs_fromGL <- Fig5A_ribogenesinAPCs_fromGL[,-1]
colnames(Fig5A_ribogenesinAPCs_fromGL)[12] <- "Significant"
Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl <- Fig5A_ribogenesinAPCs_fromGL[grepl("^Rps|^Rpl", Fig5A_ribogenesinAPCs_fromGL$GENE),]
View(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl)
unique(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$GENE)
Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl <- Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl[ Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$GENE !="Rps6ka3",]
unique(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$GENE)
colnames(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl)
unique(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$Combine)
desired_order <- c(
  "eWAT_Progenitors_ob/ob_effect", "eWAT_Progenitors_Rosi_effect", 
  "iWAT_Progenitors_ob/ob_effect", "iWAT_Progenitors_Rosi_effect", 
  "eWAT_Preadipocytes_ob/ob_effect", "eWAT_Preadipocytes_Rosi_effect",
  "iWAT_Preadipocytes_ob/ob_effect", "iWAT_Preadipocytes_Rosi_effect"
)
Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$Combine <- factor(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$Combine, levels = desired_order)
saveRDS(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl,"intermediate_data/Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl_forggplot.rds")
p3 <- ggplot(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl, aes(x = Combine, 
                                                          y = GENE, 
                                                          size = `-Log[10]FDR`, 
                                                          color = avg_log2FC, 
                                                          alpha = Significant)) +
  geom_point() +  
  scale_color_gradient2(
    low = "#2D4693", 
    mid = "white", 
    high = "#EB362E", 
    midpoint = 0,
    limits = c(min(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$avg_log2FC), max(Fig5A_ribogenesinAPCs_fromGL_onlyRpsRpl$avg_log2FC)),
    oob = scales::rescale_none
  ) + 
  scale_alpha_manual(
    values = c(1, 0.2), 
    breaks = 1,          
    labels = "TRUE",      
    guide = guide_legend(title = "FDR<0.05 \n|avg.log2FC|>0.25")
  ) +
  labs(
    y = "Cytoplasmic ribosomal genes",
    color = expression("avg.log"[2]*"FC")  
  ) +
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    axis.text.y = element_text(face = "italic"),
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"), 
    axis.ticks.x = element_line(color = "black"),  
    axis.ticks.y = element_line(color = "black"), 
    axis.ticks.length = unit(0.05, "cm"), 
    panel.border = element_blank(), 
    axis.ticks.margin = unit(0.5, "cm") 
  ) +
  guides(
    size = guide_legend(title = expression("-log"[10]*"FDR")),
    color = guide_colorbar(
      title = expression("avg.log"[2]*"FC"),
      barwidth = 1,  
      barheight = 5  
    )
  ) 
pdf("output_plot/iWATeWAT_cytoribo(Fig5A)_DEGs(padj0.05_LFC0.25)_dotplot_usingGaoyansList_06042024.pdf", width = 4.5, height = 18)
plot(p3)
dev.off()





########Pathway cite Jessica utils here
all_DEG=combined
dir.create(file.path("/Volumes/data/",CT))
colnames(all_DEG)[colnames(all_DEG) == "\"GENE\""] <- 'GENE'
colnames(all_DEG)[colnames(all_DEG) == "\"p_val\""] <- 'p_val'
colnames(all_DEG)[colnames(all_DEG) == "\"avg_log2FC\""] <- 'avg_logFC'
colnames(all_DEG)[colnames(all_DEG) == "\"p_val_adj\""] <- 'p_val_adj'
types=unique(all_DEG$CellType)
types
head(all_DEG)
for(celltype_cluster in types){
  print(celltype_cluster)
  
  DEG_df=all_DEG[all_DEG$CellType == celltype_cluster,]
  rownames(DEG_df) <-DEG_df$GENE
  DEG_df$`Cell_type` = DEG_df$CellType
  DEG_df=addInfoAndTrim(DEG_df)
  identifier=celltype_cluster
  add_direction<-TRUE
  total=pathway_enrichment(deg_list=DEG_df,FDR_threshold=NULL, pval_threshold=0.01, celltype_cluster,identifier, pVal = TRUE)
  pathways=total
  pathways = pathways[pathways$nOverlap>3,]
  Direction=c()
  meanFC=c()
  if(add_direction){
    new_overlap = c()
    for(g in pathways$Overlap){
      genes = unlist(strsplit(g, split = ","))
      UP_genes = c("UP:")
      upcount=0
      DOWN_genes = c("DOWN:")
      downcount=0
      sumFC=0
      nFC=0
      for(gene in genes){
        
        FC = DEG_df$avg_logFC[DEG_df$HUMAN==gene]
        FC = FC[!is.na(FC)]
        if( length(FC)>0){
          if(length(FC)>1){
            FC=FC[1]
            #dup_genes = append(dup_genes, gene)
            #cat("dup\n")
          }
          
          #if( length(FC)>0){}
          if(FC>0){
            upcount=upcount+1
            UP_genes = append(UP_genes, paste0(gene,","))
            sumFC=sumFC+FC
            nFC=nFC+1
          }
          else{
            downcount=downcount+1
            DOWN_genes = append(DOWN_genes, paste0(gene,","))
            sumFC=sumFC+FC
            nFC=nFC+1
          }
        }
      }
      UP_genes = append(UP_genes, DOWN_genes)
      meanFC=c(meanFC,sumFC/nFC)
      detailed_overlap = concatenate(UP_genes, mysep = "")
      new_overlap = append(new_overlap, detailed_overlap)
      if(upcount>downcount){
        Direction=c(Direction,"UP")
      }
      else{
        Direction=c(Direction,"Down")
      }
    }
    pathways$Overlap = new_overlap
    pathways$Direction = Direction
    pathways$meanFC = meanFC
  }
  
  
  write.csv(pathways,file=paste0("/Volumes/data/",CT,"/",celltype_cluster,"path.csv")) 
}


directories = list.files(path = paste0("/Volumes/data/",CT,"/"), full.names = TRUE, recursive = TRUE)
directories = directories[grep("csv",directories)]
selected=data.frame()
for(sample in directories){
  print(sample)
  pathwaylist=(read.csv(sample)) 
  selected = rbind(selected,pathwaylist)
  rm(pathwaylist)
}
head(selected)
selected$CellType=CT
selected=selected[selected$FDR<0.05,]
write.csv(selected,file=paste0("/Volumes/data/",CT,"path.csv")) 

library(dplyr)
list_of_pathway_databases = list.files(path = "/Volumes/data//Resources/", pattern = "*.txt", full.names = TRUE)
pathB = read.table("/Volumes/data/Resources/Biocarta.txt", header=T, sep='\t', check.names=F, quote=NULL)
pathH= read.table("/Volumes/data/Resources/Hallmark.txt", header=T, sep='\t', check.names=F, quote=NULL)
pathK = read.table("/Volumes/data/Resources/Kegg.txt", header=T, sep='\t', check.names=F, quote=NULL)
pathR= read.table("/Volumes/data/Resources/Reactome.txt", header=T, sep='\t', check.names=F, quote=NULL)
path=rbind(pathB,pathH,pathK,pathR)

df <- path %>%
  group_by(module) %>%
  mutate(gene_count = n())

df_summary <- df %>%
  group_by(module) %>%
  summarize(gene_count = first(gene_count))


CB <- selected %>%
  left_join(df_summary, by = c("Pathway" = "module")) %>%
  rename(pathway_size = gene_count)


############ Pathway figure
df<-read.csv("/Volumes/data/totalpathwaycombine.csv")

df$`-Log[10]FDR` = -log10(df$FDR)
df$`FDR < 0.05` = df$FDR < 0.05
dfsaved=df
unique(df$Module)
df= df[df$Tissue == "i",]
df= df[df$Module == "obUProsiDOWN",]
#"obUProsiDOWN" "obDOWNrosiUP"
df <- df %>%
  mutate(FC = recode(FC, 'meanFC1' = 'obob', 'meanFC2' = 'Rosi'))
# Rename columns
df <- df %>%
  rename(
    Effect = FC,
    log2FC = meanFC
  )
df <- df %>% 
  group_by(Pathway) %>% 
  mutate(Pathway_count = n()) %>% 
  ungroup()
###For e-Rosi down
selected_pathways_list = c("Tnfa Signaling Via Nfkb","Hypoxia","Il2 Stat5 Signaling","Kras Signaling Up","Activation Of Chaperone Genes By Atf6 Alpha",
                           "Tgf Beta Signaling","Interferon Gamma Response","Mtorc1 Signaling",
                           "Mrna Splicing","Estrogen Response Early","Il6 Jak Stat3 Signaling","Androgen Response","Mapk Signaling Pathway","Ppara Pathway",
                           "Estrogen Response Late","Nfkb Pathway","Nthi Pathway","Pi3K Akt Mtor Signaling","Egf Pathway","Hivnef pathway","Il6 Pathway","Signaling By Notch",
                           "Chemokine Signaling Pathway","Complement","Jak Stat Signaling Pathway","Ccr5 Pathway","Growth Hormone Receptor Signaling","Pi3K Akt Activation","Eif4 Pathway","Il2 Pathway","Erk5 Pathway",
                           "Hdac Pathway","Carm1 Pathway","Il7 Pathway","Erbb Signaling Pathway","Leukocyte Transendothelial Migration","Tff Pathway","Gsk3 Pathway")

###For e-Rosi up
selected_pathways_list = c("Ribosome","Peptide Chain Elongation","Myc Targets V1","Complement","Apoptosis",
                           "Ppar Signaling Pathway","Platelet Activation Signaling And Aggregation","Cholesterol Homeostasis",
                           "Adipogenesis","Complement And Coagulation Cascades","Trafficking And Processing Of Endosomal Tlr","Formation Of Atp By Chemiosmotic Coupling",
                           "Respiratory Electron Transport","P53 Dependent G1 DNA Damage Response",
                           "Synthesis Of DNA","Proteasome Pathway","Allograft Rejection","Transferrin Endocytosis And Recycling","Cell Cycle Checkpoints",
                           "Il2 Stat5 Signaling","Metabolism Of Porphyrins",
                           "Glycolysis","Peroxisome","Activation Of Kainate Receptors Upon Glutamate Binding","Glycosaminoglycan Degradation","Egfr Downregulation","Coagulation","Kras Signaling Up")

###For i-Rosi up
selected_pathways_list = c("Ribosome","Peptide Chain Elongation","Myc Targets V1","Formation Of Atp By Chemiosmotic Coupling","Respiratory Electron Transport",
                           "Adipogenesis","P53 Dependent G1 DNA Damage Response","Synthesis Of DNA",
                           "Apoptosis","Etc Pathway","Eif Pathway","Fatty Acid Metabolism",
                           "Metabolism Of Amino Acids And Derivatives",
                           "Platelet Activation Signaling And Aggregation","Rho Pathway","Mta3 Pathway","Salmonella Pathway","Cdc42Rac Pathway",
                           "Ras Pathway","Hypoxia","Cholesterol Homeostasiss",
                           "Glycolysis Gluconeogenesis","Lysosome","Complement","Regulatory RNA Pathways","Citric Acid Cycle Tca Cycle","Ccr5 Pathway", "Tnfa Signaling Via Nfkb","Metabolism Of Non Coding RNA","Regulation Of Actin Cytoskeleton")

###For i Rosi down
selected_pathways_list = c("Epithelial Mesenchymal Transition","Coagulation","Apoptosis","Extracellular Matrix Organization","Complement",
                           "Hypoxia","Tnfa Signaling Via Nfkb","Mtorc1 Signaling",
                           "Interferon Gamma Response","Mitotic Spindle","Ncam1 Interactions","Signaling By Pdgf",
                           "Toll Receptor Cascades",
                           "Fatty Acid Metabolism")
selected_pathways = data.frame(stringsAsFactors = FALSE)
for(path in selected_pathways_list){
  temp = df[df$Pathway==path,]
  selected_pathways = rbind(selected_pathways, temp)
}
selected_pathways <- selected_pathways %>% 
  group_by(Pathway) %>% 
  mutate(Pathway_count = n()) %>% 
  ungroup()
write.csv(df,"/Volumes/data/iWAT_RosiDownregulated_pathway.csv")
write.csv(selected_pathways,"/Volumes/data/iWAT_RosiDownregulated_pathway_selected.csv")

my_plot <- ggplot(selected_pathways, aes(x = Pathway, y = CellType)) +
  geom_point(aes(size = Enrichment, colour = log2FC, shape = Effect), stroke = 0) +
  scale_shape_manual(values = c("obob_Up" = "\u25D6", "Rosi_Down" = "\u25D7")) +
  scale_colour_gradient2(
    low = "#263E85", mid = "white", high = "#D41111", 
    midpoint = 0, limits = c(-1, 1),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(3, 12), breaks = c(5, 10, 30)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial Unicode MS"))
pdf("/Volumes/data/celltypecombine/PathwayRosiUp.pdf", height=8,width = 15)
print(my_plot)
dev.off()




###Polysome/totalRNA speicifc and shared pathway plots
### --- 1. Get polysome/totalRNA-specific and shared DEGs in Fig6 --- ###
sig_TotalRNA <- read_excel("./raw_data/sig_TotalRNA.xlsx",  skip = 1)
sig_TotalRNA <- sig_TotalRNA[,-1]
sig_Polysome <- read_excel("./raw_data/sig_Polysome.xlsx", 
                           skip = 1)
sig_Polysome <- sig_Polysome[,-1]
## Total RNA specific 
TotalRNA_specific <- anti_join(sig_TotalRNA, sig_Polysome, by = "GENE")
## Polysome specific 
Polysome_specific <- anti_join(sig_Polysome, sig_TotalRNA, by = "GENE")
## Shared
TotalRNA_Polysome_shared <- inner_join(sig_TotalRNA, sig_Polysome, by = "GENE")
names(TotalRNA_Polysome_shared) <- sub("\\.x$", "_TotalRNA", names(TotalRNA_Polysome_shared))
names(TotalRNA_Polysome_shared) <- sub("\\.y$", "_Polysome", names(TotalRNA_Polysome_shared))
save(TotalRNA_specific,Polysome_specific,TotalRNA_Polysome_shared, file = "./intermediate_data/SharedORSpecific_degs.RData")

### --- 2. Pathway analysis on polysome/totalRNA-shared 85 DEGs but with totalRNA values --- ###
load("./intermediate_data/SharedORSpecific_degs.RData")
source("script/Jessica_pathway_GL_YZmodified.r")
shared_TotalRNAval <- TotalRNA_Polysome_shared
shared_TotalRNAval <- shared_TotalRNAval[,-c(8:13)]
shared_TotalRNAval$p_val <- shared_TotalRNAval$pvalue_TotalRNA
shared_TotalRNAval$p_val_adj <- shared_TotalRNAval$padj_TotalRNA
shared_TotalRNAval$Cell_type <- "iWAT"
shared_TotalRNAval$Comparison <- "Rosi.v.Veh"
shared_TotalRNAval$avg_logFC <- shared_TotalRNAval$log2FoldChange_TotalRNA
DEG_df <- shared_TotalRNAval
rownames(DEG_df) <-DEG_df$GENE
DEG_df=addInfoAndTrim(DEG_df, 
                      p_cutoff=0.01, 
                      multiple_comparisons=TRUE, 
                      convertToHuman=TRUE, 
                      lfc_cutoff=0.1) # Used all default params here.
View(DEG_df)
celltype_cluster <- "shared_degs_wTotalRNA.values"
identifier=celltype_cluster
add_direction <- TRUE
pathway_sharedDEGs_wTotalRNAValue = pathway_enrichment(deg_list=DEG_df,
                                                       FDR_threshold=NULL, 
                                                       pval_threshold=0.01, 
                                                       celltype_cluster,
                                                       identifier, 
                                                       pVal = TRUE)  # pVal this is our method of switching off between FDR and pval threshold
## Filter pathway results to only keep ones with nOverlap > 3
pathway_sharedDEGs_wTotalRNAValue = pathway_sharedDEGs_wTotalRNAValue[pathway_sharedDEGs_wTotalRNAValue$nOverlap>3,]
## Assign meanFC and direction to each pathway 
Direction=c()
meanFC=c() # meanFC here represents the mean of the log2FC of genes enriched in the pathway 
if(add_direction){
  new_overlap = c()
  for(g in pathway_sharedDEGs_wTotalRNAValue$Overlap){
    genes = unlist(strsplit(g, split = ","))
    UP_genes = c("UP:")
    upcount=0
    DOWN_genes = c("DOWN:")
    downcount=0
    sumFC=0
    nFC=0
    for(gene in genes){
      FC = DEG_df$avg_logFC[DEG_df$HUMAN==gene]
      FC = FC[!is.na(FC)]
      if( length(FC)>0){
        if(length(FC)>1){
          FC=FC[1]
          #dup_genes = append(dup_genes, gene)
          #cat("dup\n")
        }
        #if( length(FC)>0){}
        if(FC>0){
          upcount=upcount+1
          UP_genes = append(UP_genes, paste0(gene,","))
          sumFC=sumFC+FC
          nFC=nFC+1
        }
        else{
          downcount=downcount+1
          DOWN_genes = append(DOWN_genes, paste0(gene,","))
          sumFC=sumFC+FC
          nFC=nFC+1
        }
      }
    }
    # UP_genes = append(UP_genes, DOWN_genes)
    meanFC=c(meanFC,sumFC/nFC)
    # detailed_overlap = paste(UP_genes, collapse = "")
    # new_overlap = append(new_overlap, detailed_overlap)
    # if(upcount>downcount){
    #   Direction=c(Direction,"UP")
    # }
    # else{
    #   Direction=c(Direction,"Down")
    # }
  }
  # pathway_sharedDEGs_wTotalRNAValue$Overlap = new_overlap
  pathway_sharedDEGs_wTotalRNAValue$meanFC = meanFC
}
pathway_sharedDEGs_wTotalRNAValue$Direction <- ifelse(pathway_sharedDEGs_wTotalRNAValue$meanFC > 0, "UP", "DOWN")
pathway_sharedDEGs_wTotalRNAValue.sig <- pathway_sharedDEGs_wTotalRNAValue[pathway_sharedDEGs_wTotalRNAValue$FDR < 0.05,]
View(pathway_sharedDEGs_wTotalRNAValue.sig)
pathway_sharedDEGs_wTotalRNAValue.sig$Pathway <- sapply(pathway_sharedDEGs_wTotalRNAValue.sig$Pathway, transform_pathway_name)
write.table(pathway_sharedDEGs_wTotalRNAValue.sig,file="./output_data/shared85deg_enrichedpathways_wTotalRNAValues/pathway_sharedDEGs_wTotalRNAValue_sig(nOverlap.4_padj0.05)",quote = FALSE,sep='\t',row.names=F) 

### --- 3. Pathway analysis on polysome/totalRNA-shared 85 DEGs but with Polysome values --- ###
load("./intermediate_data/SharedORSpecific_degs.RData")
shared_Polysomeval <- TotalRNA_Polysome_shared
shared_Polysomeval <- shared_Polysomeval[,-c(2:7)]
shared_Polysomeval$p_val <- shared_Polysomeval$pvalue_Polysome
shared_Polysomeval$p_val_adj <- shared_Polysomeval$padj_Polysome
shared_Polysomeval$Cell_type <- "iWAT"
shared_Polysomeval$Comparison <- "Rosi.v.Veh"
shared_Polysomeval$avg_logFC <- shared_Polysomeval$log2FoldChange_Polysome
DEG_df <- shared_Polysomeval
rownames(DEG_df) <-DEG_df$GENE
DEG_df=addInfoAndTrim(DEG_df, 
                      p_cutoff=0.01, 
                      multiple_comparisons=TRUE, 
                      convertToHuman=TRUE, 
                      lfc_cutoff=0.1) # Used all default params here.
View(DEG_df)
celltype_cluster <- "shared_degs_wPolysome.values"
identifier=celltype_cluster
add_direction <- TRUE
pathway_sharedDEGs_wPolysomeValue = pathway_enrichment(deg_list=DEG_df,
                                                       FDR_threshold=NULL, 
                                                       pval_threshold=0.01, 
                                                       celltype_cluster,
                                                       identifier, 
                                                       pVal = TRUE)  # pVal this is our method of switching off between FDR and pval threshold
## Filter pathway results to only keep ones with nOverlap > 3
pathway_sharedDEGs_wPolysomeValue = pathway_sharedDEGs_wPolysomeValue[pathway_sharedDEGs_wPolysomeValue$nOverlap>3,]
## Assign meanFC and direction to each pathway 
meanFC=c() # meanFC here represents the mean of the log2FC of genes enriched in the pathway 
if(add_direction){
  new_overlap = c()
  for(g in pathway_sharedDEGs_wPolysomeValue$Overlap){
    genes = unlist(strsplit(g, split = ","))
    UP_genes = c("UP:")
    upcount=0
    DOWN_genes = c("DOWN:")
    downcount=0
    sumFC=0
    nFC=0
    for(gene in genes){
      FC = DEG_df$avg_logFC[DEG_df$HUMAN==gene]
      FC = FC[!is.na(FC)]
      if( length(FC)>0){
        if(length(FC)>1){
          FC=FC[1]
          #dup_genes = append(dup_genes, gene)
          #cat("dup\n")
        }
        #if( length(FC)>0){}
        if(FC>0){
          upcount=upcount+1
          UP_genes = append(UP_genes, paste0(gene,","))
          sumFC=sumFC+FC
          nFC=nFC+1
        }
        else{
          downcount=downcount+1
          DOWN_genes = append(DOWN_genes, paste0(gene,","))
          sumFC=sumFC+FC
          nFC=nFC+1
        }
      }
    }
    meanFC=c(meanFC,sumFC/nFC)
  }
  pathway_sharedDEGs_wPolysomeValue$meanFC = meanFC
}

pathway_sharedDEGs_wPolysomeValue$Direction <- ifelse(pathway_sharedDEGs_wPolysomeValue$meanFC > 0, "UP", "DOWN")
pathway_sharedDEGs_wPolysomeValue.sig <- pathway_sharedDEGs_wPolysomeValue[pathway_sharedDEGs_wPolysomeValue$FDR < 0.05,]
View(pathway_sharedDEGs_wPolysomeValue.sig)
pathway_sharedDEGs_wPolysomeValue.sig$Pathway <- sapply(pathway_sharedDEGs_wPolysomeValue.sig$Pathway, transform_pathway_name)
write.table(pathway_sharedDEGs_wPolysomeValue.sig,file="./output_data/shared85deg_enrichedpathways_wPolysomeValues/pathway_sharedDEGs_wPolysomeValue_sig(nOverlap.4_padj0.05)",quote = FALSE,sep='\t',row.names=F) 

### --- 4. Pathway visualization derived from section 2&3 --- ###
pathway_sharedDEGs_wTotalRNAValue.sig <- read.delim("./output_data/shared85deg_enrichedpathways_wTotalRNAValues/pathway_sharedDEGs_wTotalRNAValue_sig(nOverlap.4_padj0.05)")
pathway_sharedDEGs_wTotalRNAValue.sig$ModuleSize <- NULL
pathway_sharedDEGs_wPolysomeValue.sig <- read.delim("./output_data/shared85deg_enrichedpathways_wPolysomeValues/pathway_sharedDEGs_wPolysomeValue_sig(nOverlap.4_padj0.05)")
colnames(pathway_sharedDEGs_wPolysomeValue.sig)
pathway_sharedDEGs_wPolysomeValue.sig$ModuleSize <- NULL
# also load the original pathway results in Fig6F 
Polysome_TotalRNA_separate_pathways <- read_excel("./raw_data/Polysome_TotalRNA_separate_pathways.xlsx")
colnames(Polysome_TotalRNA_separate_pathways)
Polysome_TotalRNA_separate_pathways$Pathway <- sapply(Polysome_TotalRNA_separate_pathways$Pathway, transform_pathway_name)
# Combine all pathways from totalRNA, Polysome, shared by total rna value, and shared by polysome value
combined_pathways_for_plotting <- rbind(rbind(pathway_sharedDEGs_wTotalRNAValue.sig,pathway_sharedDEGs_wPolysomeValue.sig),Polysome_TotalRNA_separate_pathways)
View(combined_pathways_for_plotting)
## Apply Gaoyan's visualization script
combined_pathways_for_plotting$`neg.log10FDR` = -log10(combined_pathways_for_plotting$FDR)
combined_pathways_for_plotting$`FDR < 0.05` = combined_pathways_for_plotting$FDR < 0.05
range (combined_pathways_for_plotting$meanFC)
combined_pathways_for_plotting <-combined_pathways_for_plotting[order(combined_pathways_for_plotting$meanFC),]
order=unique(combined_pathways_for_plotting$`Pathway`)
combined_pathways_for_plotting$`Pathway` = factor(combined_pathways_for_plotting$`Pathway`, levels =order)
####control min and max color coding
num = dim(combined_pathways_for_plotting)[1]
combined_pathways_for_plotting1=combined_pathways_for_plotting
combined_pathways_for_plotting = data.frame(stringsAsFactors = FALSE)
for(i in 1:num){
  temp = combined_pathways_for_plotting1[i,]
  if (temp$meanFC<(-2)) {
    temp$meanFC=(-2)
  }
  combined_pathways_for_plotting = rbind(combined_pathways_for_plotting, temp)
}
range (combined_pathways_for_plotting$meanFC)
num=dim(combined_pathways_for_plotting)[1]
combined_pathways_for_plotting1=combined_pathways_for_plotting
combined_pathways_for_plotting = data.frame(stringsAsFactors = FALSE)
for(i in 1:num){
  temp = combined_pathways_for_plotting1[i,]
  if (temp$meanFC>2) {
    temp$meanFC=2
  }
  combined_pathways_for_plotting = rbind(combined_pathways_for_plotting, temp)
}
range (combined_pathways_for_plotting$meanFC)
combined_pathways_for_plotting_pathways=combined_pathways_for_plotting
combined_pathways_for_plotting_pathways$Module <- factor(x = combined_pathways_for_plotting_pathways$Module ,
                                                         levels = c("totalRNA","Polysome","shared_degs_wTotalRNA.values", "shared_degs_wPolysome.values"))
View(combined_pathways_for_plotting_pathways)
write.xlsx(combined_pathways_for_plotting_pathways,"output_data/NewFig6F_all_pathways.xlsx")
D <- ggplot(combined_pathways_for_plotting_pathways, aes(x=Pathway, y=`Module`, size=`Enrichment`, color=meanFC)) + theme_bw() +
  geom_point(aes(shape=`FDR < 0.05`))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        #axis.text.x = element_text(angle = 90, margin = margin(t=5), size=14, hjust = 1,vjust = 1),
        #axis.text.y = element_text(angle = 45,size = 14),
        axis.text.x = element_text(angle = 90, margin = margin(t=5), size=14, hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),         
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.position = "bottom") +
  labs(size=expression('Enrichment Score'), color=expression('avg_log2FC')) +
  scale_color_gradient2(low = "blue",  high = "red", mid = "gray") +
  scale_size(range = c(3, 12), breaks = c(0.5,10,25,50,75,100)) +
  scale_shape_manual(values=c(16,16)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) 
pdf("output_plot/pathway_dotplot_sharedDEG.pathways_added_horizontal.pdf", height=9,width = 20)
print(D)
dev.off()
## switch to vertical dotplot
D2 <- ggplot(combined_pathways_for_plotting_pathways, aes(y=Pathway, x=`Module`, size=`Enrichment`, color=meanFC)) + theme_bw() +
  geom_point(aes(shape=`FDR < 0.05`)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, size = 14, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y.right = element_text(angle = 0, size = 14, face="plain"),
        axis.text = element_text(face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.y.right = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.position = "right") +
  labs(size=expression('Enrichment Score'), color=expression('avg_log2FC')) +
  scale_color_gradient2(low = "blue",  high = "red", mid = "gray") +
  scale_size(range = c(3, 12), breaks = c(0.5,10,25,50,75,100)) +
  scale_shape_manual(values=c(16,16)) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  scale_y_discrete(position = "right")
pdf("output_plot/pathway_dotplot_sharedDEG.pathways_added_vertical.pdf", height=17, width = 10)
print(D2)
dev.off()





###Proteomics plot
### --- 1. Load in Proteomics data from Mirian --- ###
iWAT_Rosi_v_Veh_DEprotein_unfiltered_df <- read_excel("./raw_data/iWAT_Rosi.v.Veh_Proteomics_DE_results.xlsx", 
                                                      sheet = "All_Results")
View(iWAT_Rosi_v_Veh_DEprotein_unfiltered_df)
## remove unwanted columns 
iWAT_Rosi_v_Veh_DEprotein_unfiltered_df <- iWAT_Rosi_v_Veh_DEprotein_unfiltered_df[,-c(4,6:9)]
saveRDS(iWAT_Rosi_v_Veh_DEprotein_unfiltered_df,"intermediate_data/iWAT_Rosi_v_Veh_DEprotein_unfiltered_df(unwanted.cols.removed).rds")

### --- 2. Volcano plot for DE proteins --- ###
iWAT_Rosi_v_Veh_DEprotein_unfiltered_df <- readRDS("intermediate_data/iWAT_Rosi_v_Veh_DEprotein_unfiltered_df(unwanted.cols.removed).rds") 
View(iWAT_Rosi_v_Veh_DEprotein_unfiltered_df)
df <- iWAT_Rosi_v_Veh_DEprotein_unfiltered_df
df <- df[,-1]
df$Protein <- toupper(df$`Gene Name`)
colnames(df) <- c("Gene","log2FC","padj","Protein")
df <- df %>%
  mutate(Significance = ifelse(padj < 0.05 & abs(log2FC) > 2, "Significant", "Not Significant"))
df$Direction <- df$Significance
df[df$Significance == "Significant",]$Direction <- ifelse(df[df$Significance == "Significant",]$log2FC > 0, "UP", "DOWN")
nrow(df[df$Direction != "Not Significant",])
# Create the volcano plot
p <- ggplot(df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = Direction), size = 1) +
  scale_color_manual(values = c("UP" = "red", "Not Significant" = "grey", "DOWN" = "blue")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(p.adj)"), title = "iWAT Rosi vs Veh Proteins") +
  theme(legend.position = "none")
proteins_to_label <- c("ADIPOQ", "CPT1B", "CIDEC", "AQP7", "RBP7", "FASN", "PLIN1", "PLIN4", "ADIPOR2", "FABP4")
volcano_plot <- p +
  geom_text_repel(data = df %>% filter(Protein %in% proteins_to_label), 
                  aes(label = Protein),
                  size = 4, 
                  box.padding = unit(1, "lines"),
                  point.padding = unit(1, "lines"),
                  segment.color = 'grey50',
                  # nudge_y = 0.5,
                  # nudge_x = 0.5,
                  max.overlaps = Inf)
pdf("output_plot/Volcano_sigDEP(padj0.05_absLFC2)_colored.pdf")
volcano_plot
dev.off()

### --- 3. Dotplot heatmap for sig. cyto.ribosomal DE proteins x--- ###
iWAT_Rosi_v_Veh_DEprotein_unfiltered_df <- readRDS("intermediate_data/iWAT_Rosi_v_Veh_DEprotein_unfiltered_df(unwanted.cols.removed).rds") 
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_unfiltered_df <- iWAT_Rosi_v_Veh_DEprotein_unfiltered_df[grepl("^Rpl|^Rps",iWAT_Rosi_v_Veh_DEprotein_unfiltered_df$`Gene Name`),]
saveRDS(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_unfiltered_df,"intermediate_data/iWAT_Rosi_v_Veh_ribosomal_DEprotein_unfiltered_df(unwanted.cols.removed).rds")
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df <- iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_unfiltered_df[abs(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_unfiltered_df$`iWAT_Rosi_vs_iWAT_Veh_log2 fold change`) > 0.25 & iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_unfiltered_df$iWAT_Rosi_vs_iWAT_Veh_p.adj < 0.05,]
View(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df)
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$neg.log10padj <- -log10(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$iWAT_Rosi_vs_iWAT_Veh_p.adj)
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$Comparison <- "iWAT_Rosi.v.Veh"
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$ProteinName <- toupper(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$`Gene Name`)
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df <- iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df[-which( iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$`Gene Name` %in% c("Rps6ka3","Rps6ka4")),]
saveRDS(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df,"intermediate_data/iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig(padj0.05_LFC0.25)_df.rds")
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$log2FC <- iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$`iWAT_Rosi_vs_iWAT_Veh_log2 fold change`
iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$Significant <- iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$iWAT_Rosi_vs_iWAT_Veh_p.adj < 0.05
p <- ggplot(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df, aes(x = Comparison, 
                                                              y = ProteinName, 
                                                              size = neg.log10padj, 
                                                              color = log2FC, 
                                                              alpha = Significant)) +
  geom_point() +  
  scale_color_gradient2(
    low = "#2D4693", 
    mid = "white", 
    high = "#EB362E", 
    midpoint = 0,
    limits = c(min(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$log2FC, na.rm = TRUE), max(iWAT_Rosi_v_Veh_cyto.Ribo_Proteins_DE_sig_df$log2FC, na.rm = TRUE)),
    oob = scales::rescale_none  
  ) + 
  scale_alpha_manual(
    values = c(1, 0.2),  
    breaks = 1,         
    labels = "TRUE",     
    guide = guide_legend(title = "p.adj<0.05\n|log2FC|>0.25")
  ) +
  labs(
    y = "Cytosolic ribosomal proteins",
    color = expression("log"[2]*"FC")  
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "white"), 
    plot.background = element_rect(fill = "white", colour = "white"),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),  
    axis.ticks.x = element_line(color = "black"), 
    axis.ticks.y = element_line(color = "black"),  
    axis.ticks.length = unit(0.05, "cm"),  
    panel.border = element_blank(),  
    axis.ticks.margin = unit(0.5, "cm")  
  ) +
  guides(
    size = guide_legend(title = expression("-log"[10]*"p.adj")),
    color = guide_colorbar(
      title = expression("log"[2]*"FC"),
      barwidth = 1,  # Reduce bar width
      barheight = 5  # Standard height
    )
  ) 
pdf("output_plot/iWAT_Rosi.v.Veh_DE.cytoribo_sig(padj0.05_LFC0.25)_dotplot.pdf", width = 2.5, height = 10)
plot(p)
dev.off()

### --- 4. Dotplot heatmap for sig. mito.ribosomal DE proteins --- ###
iWAT_Rosi_v_Veh_DEprotein_unfiltered_df <- readRDS("intermediate_data/iWAT_Rosi_v_Veh_DEprotein_unfiltered_df(unwanted.cols.removed).rds") 
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_unfiltered_df <- iWAT_Rosi_v_Veh_DEprotein_unfiltered_df[grepl("^Mrpl|^Mrps",iWAT_Rosi_v_Veh_DEprotein_unfiltered_df$`Gene Name`),]
saveRDS(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_unfiltered_df,"intermediate_data/iWAT_Rosi_v_Veh_mitoribosomal_DEprotein_unfiltered_df(unwanted.cols.removed).rds")
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df <- iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_unfiltered_df[abs(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_unfiltered_df$`iWAT_Rosi_vs_iWAT_Veh_log2 fold change`) > 0.25 & iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_unfiltered_df$iWAT_Rosi_vs_iWAT_Veh_p.adj < 0.05,]
View(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df)
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$neg.log10padj <- -log10(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$iWAT_Rosi_vs_iWAT_Veh_p.adj)
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$Comparison <- "iWAT_Rosi.v.Veh"
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$ProteinName <- toupper(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$`Gene Name`)
saveRDS(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df,"intermediate_data/iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig(padj0.05_LFC0.25)_df.rds")
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$log2FC <- iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$`iWAT_Rosi_vs_iWAT_Veh_log2 fold change`
iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$Significant <- iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$iWAT_Rosi_vs_iWAT_Veh_p.adj < 0.05
p <- ggplot(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df, aes(x = Comparison, 
                                                              y = ProteinName, 
                                                              size = neg.log10padj, 
                                                              color = log2FC, 
                                                              alpha = Significant)) +
  geom_point() +  
  scale_color_gradient2(
    low = "#2D4693", 
    mid = "white", 
    high = "#EB362E", 
    midpoint = 0,
    limits = c(min(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$log2FC, na.rm = TRUE), max(iWAT_Rosi_v_Veh_mito.Ribo_Proteins_DE_sig_df$log2FC, na.rm = TRUE)),
    oob = scales::rescale_none 
  ) +  
  scale_alpha_manual(
    values = c(1, 0.2), 
    breaks = 1,           
    labels = "TRUE",      
    guide = guide_legend(title = "p.adj<0.05\n|log2FC|>0.25")
  ) +
  labs(
    y = "Mitochondrial ribosomal proteins",
    color = expression("log"[2]*"FC") )
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better visibility
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),  
    axis.ticks.x = element_line(color = "black"),  
    axis.ticks.y = element_line(color = "black"), 
    axis.ticks.length = unit(0.05, "cm"),  
    panel.border = element_blank(),  
    axis.ticks.margin = unit(0.5, "cm")  
  ) +
  guides(
    size = guide_legend(title = expression("-log"[10]*"p.adj")),
    color = guide_colorbar(
      title = expression("log"[2]*"FC"),
      barwidth = 1,  
      barheight = 5 
    )
  ) 
pdf("output_plot/iWAT_Rosi.v.Veh_DE.mitoribo_sig(padj0.05_LFC0.25)_dotplot.pdf", width = 2.5, height = 15)
plot(p)
dev.off()





