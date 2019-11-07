library(Seurat, lib = "/frazer01/home/mdonovan/alt_r_library/")
library(dplyr)
library(biomaRt)
library(stringr)

meta = read.csv("/projects/PPC/analysis/ppc_pilot/tables/ppc_pilot_metadata.csv")

setwd("/projects/PPC/analysis/ppc_pilot/data/")
dir.create("qc", showWarnings = FALSE)
dir.create("robjs", showWarnings = FALSE)
dir.create("UMAP", showWarnings = FALSE)
dir.create("differential_expression", showWarnings = FALSE)
dir.create("gene_expression", showWarnings = FALSE)
dir.create("integrated", showWarnings = FALSE)
dir.create("prediction", showWarnings = FALSE)

marker_genes <- c("INS", "NKX6-1", "NPTX2", "ITGA1", "PCSK1",
                         "GCG", "DPP4",   "ARX",   "IRX2",  "ETV1",
                         "TPH1", "SLC18A1", "LMX1A", "ADRA2A", "MME",
                         "CHGA", "FEV", "SYP", "CPE", "CDKN1C",
                         "SOX9", "PTF1A", "PDX1",
                         "TOP2A", "AURKB", "ISL1", "PAX4", "NEUROG3")

#### load melvin data #####
load("/projects/PPC/analysis/ppc_pilot/data/melvin_data/Stages_3_to_6.x1_S3c.Robj")
S3c = tiss
load("/projects/PPC/analysis/ppc_pilot/data/melvin_data/Stages_3_to_6.x1_S4c.Robj")
S4c = tiss
load("/projects/PPC/analysis/ppc_pilot/data/melvin_data/Stages_3_to_6.x1_S5c.Robj")
S5c = tiss
load("/projects/PPC/analysis/ppc_pilot/data/melvin_data/Stages_3_to_6.x1_S6c.Robj")
S6c = tiss
load("/projects/PPC/analysis/ppc_pilot/data/melvin_data/ESB_only_integrated.Robj")
scb_integrated = integrated

S3c@meta.data$Celltype <- S3c@meta.data$Assigned_cluster
S4c@meta.data$Celltype <- S4c@meta.data$Assigned_cluster
S5c@meta.data$Celltype <- S5c@meta.data$Assigned_cluster
S6c@meta.data$Celltype <- S6c@meta.data$Assigned_cluster
S3c@meta.data$set    <- S3c@meta.data$Lib_prep_batch
S4c@meta.data$set    <- S4c@meta.data$Lib_prep_batch
S5c@meta.data$set    <- S5c@meta.data$Lib_prep_batch
S6c@meta.data$set    <- S6c@meta.data$Lib_prep_batch
###########################


for (m in seq(1, nrow(meta))){
    
    sample_name = meta[m, "name"]
    path        = meta[m, "path"]
    

    tiss <- Read10X(data.dir = path)
    tiss <- CreateSeuratObject(counts = tiss, project = sample_name)
    
    tiss[["percent.mt"]] <- PercentageFeatureSet(tiss, pattern = "^MT")
    
    png(file = paste("./qc/", sample_name, "_umi.vs.mt.png",sep = ""), width = 5, height = 5, units = "in", res = 300)
    plot(tiss@meta.data$nCount_RNA, tiss@meta.data$percent.mt, pch = 20, cex = .5, xlab = "nUMI", ylab = "%MT", las = 1) 
    abline(h = 10, col = "red")
    dev.off()
    
    png(file = paste("./qc/", sample_name, "_umi.vs.genes.png",sep = ""), width = 5, height = 5, units = "in", res = 300)
    plot(tiss@meta.data$nCount_RNA, tiss@meta.data$nFeature_RNA, pch = 20, cex = .5, xlab = "nUMI", ylab = "n Unique Genes", las = 1) 
    abline(h = 500, col = "red")
    abline(h = 6000, col = "red")
    dev.off()
    
    save(tiss, paste("./robjs/", sample_name, "_seurat_unfiltered.robj",sep = ""))
    
    tiss <- subset(tiss, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
    tiss <- NormalizeData(tiss)
    tiss <- FindVariableFeatures(tiss, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(tiss)
    tiss <- ScaleData(tiss, features = all.genes)
    
    tiss <- RunPCA(tiss, features = VariableFeatures(object = tiss))
    
    pdf(paste("./qc/", sample_name, "_elbow_plot.pdf",sep = ""), height = 5, width = 10)
    ElbowPlot(tiss, ndims = 50)
    dev.off()
    
    n.pcs = 20
    res.used <- .5

    tiss <- FindNeighbors(tiss, dims = 1:n.pcs)
    tiss <- FindClusters( tiss, resolution = res.used)
    tiss <- RunUMAP(tiss, dims = 1:n.pcs, min_dist = .8, seed = 10)
    
    png(file = paste("./UMAP/", sample_name, "_cluster_UMAP.png",sep = ""), width = 10, height = 5, units = "in", res = 300)
    DimPlot(tiss, reduction.use = "umap", label = TRUE, pt.size = 1)
    dev.off()
    
    pdf(paste("./gene_expression/", sample_name, "_heatmap.pdf",sep = ""), height = 5, width = 10)
    DoHeatmap(tiss, features = marker_genes) + NoLegend()
    dev.off()
    
    pdf(paste("./gene_expression/", sample_name, "_heatmap.pdf",sep = ""), height = 5, width = 10)
    FeaturePlot(object = tiss, features = toupper(unique(marker_genes)))
    dev.off()
    
    clusters = seq(1, length(unique(tiss@meta.data$RNA_snn_res.0.5)))
    for( i in clusters){
        clusterX.markers <- FindMarkers(tiss, ident.1 = i, min.pct = 0.25)
        write.csv(clusterX.markers, file = paste("./differential_expression/", sample_name, "_cluster_", i, ".csv", sep = ""))
    }
    
    save(tiss, paste("./robjs/", sample_name, "_seurat_filtered.robj",sep = ""))
    
    
    ### integrate
    tiss@meta.data$set <- sample_name
    tiss@meta.data$Celltype <- tiss@meta.data$RNA_snn_res.0.5
    
    data2integrate = list()
    data2integrate[["tiss"]] = tiss
    data2integrate[["S3c"]]    = S3c
    data2integrate[["S4c"]]    = S4c
    data2integrate[["S5c"]]    = S5c
    data2integrate[["S6c"]]    = S6c
    
    anchors <- FindIntegrationAnchors(object.list = data2integrate, dims = 1:30)
    integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated, verbose = FALSE)
    integrated <- RunPCA(   integrated, npcs = 30, verbose = FALSE)
    integrated <- RunUMAP(  integrated, reduction = "pca", dims = 1:30, seed = 10, min.dist = 0.1, n.neighbors = 30)
    
    save(integrated, paste("./robjs/", sample_name, "_seurat_SCBintegrated.robj",sep = ""))   
    
    png(file = paste("./UMAP/", sample_name, "_SCBintegrated_byset_UMAP.png",sep = ""), width = 10, height = 5, units = "in", res = 300)
    DimPlot(integrated, group.by = "set", label = T, pt.size = .15)
    dev.off()
    
    png(file = paste("./UMAP/", sample_name, "_SCBintegrated_bycelltype_UMAP.png",sep = ""), width = 10, height = 5, units = "in", res = 300)
    DimPlot(integrated, group.by = "Celltype", label = T, pt.size = .15)
    dev.off()
    
    
    ######project########
    query = tiss
    anchors <- FindTransferAnchors(reference = scb_integrated, query = query, dims = 1:30)
    predictions <- TransferData(anchorset = anchors, refdata = scb_integrated$Standard_celltype, dims = 1:30)
    query <- AddMetaData(query, metadata = predictions)
    
    save(query, paste("./robjs/", sample_name, "_seurat_SCBprojection.robj",sep = ""))    
    
    png(file = paste("./UMAP/", sample_name, "_SCB_projected_UMAP.png",sep = ""), width = 10, height = 5, units = "in", res = 300)
    DimPlot(query, group.by = "predicted.id", label = T)
    dev.off()
    
    ############
    
}