library(Seurat)
library(tidyverse)
library(harmony)
library(xlsx)
library(patchwork)

anno=read.table("npc_anno.txt",header = T,sep = "\t",stringsAsFactors = F)
anno$CB=str_replace(anno$CB,"-","_")
count=readRDS("all_sparse_mat.rds")

### 标准流程 ##########################################################
allseu = CreateSeuratObject(counts = count)
allseu[["percent.mt"]] <- PercentageFeatureSet(allseu, pattern = "^MT-")
allseu <- NormalizeData(allseu, normalization.method = "LogNormalize", scale.factor = 10000)
allseu <- FindVariableFeatures(allseu, selection.method = "vst", nfeatures = 2000)
allseu <- ScaleData(allseu)
allseu <- RunPCA(allseu, npcs = 50, verbose = FALSE)

### 添加注释 ##########################################################
allseu@meta.data$CB=rownames(allseu@meta.data)
allseu@meta.data = inner_join(allseu@meta.data,anno,by="CB")
rownames(allseu@meta.data)=allseu@meta.data$CB

### 不去批次，直接降维&聚类 ###########################################
allseu <- FindNeighbors(allseu, dims = 1:20)
allseu <- FindClusters(allseu, resolution = 0.5)
allseu <- RunUMAP(allseu, dims = 1:20)
allseu <- RunTSNE(allseu, dims = 1:20)

### 简单可视化 ########################################################
DimPlot(allseu, reduction = "tsne",group.by = "cluster",pt.size = 1)+
  DimPlot(allseu, reduction = "tsne",group.by = "sample",pt.size = 1)
# saveRDS(allseu,"allseu.withbe.1112.rds")
# allseu=readRDS("allseu.withbe.1112.rds")

# ### harmony去批次
# allseu=allseu %>% RunHarmony("sample", plot_convergence = TRUE)
# 
# ### 降维&聚类
# allseu <- allseu %>%
#   RunUMAP(reduction = "harmony", dims = 1:20) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# allseu <- allseu %>%
#   RunTSNE(reduction = "harmony", dims = 1:20)
# 
# ### 简单可视化
# DimPlot(allseu, reduction = "tsne",group.by = "cluster",pt.size = 0.5)+
#   DimPlot(allseu, reduction = "tsne",group.by = "sample",pt.size = 0.5)
# ggsave("withoutbe.pdf",width = 40,height = 20,units = "cm")

### 设置配色 ###########################################################
#Snipaste
library(RColorBrewer)
library(scales)

color_cluster=c("#e61a2e","#1f6cb9","#3fb04b","#873f98","#f27024")
names(color_cluster)=c("B cell","CAFs","Epithelial cell","Myeloid","T cell")
show_col(color_cluster)

color_sample=c("#e5192c","#3a77b7","#3cac4c","#813c93","#f36c24",
               "#37b8c3","#a54922","#6b7627","#28996b",
               "#965b6a","#e9148f","#595b5e","#76c3ad",
               "#80d08a","#d29099","#f2e010")
names(color_sample)=sort(unique(allseu@meta.data$sample))
show_col(color_sample)

### figure1 b ##########################################################
### 基因集打分
marker.gene.set=read.xlsx("gene_set.xlsx",sheetIndex = 1,header = F)
colnames(marker.gene.set)=c("gene","set")
marker.gene.set=marker.gene.set[marker.gene.set$gene %in% rownames(allseu),]

for (i in unique(marker.gene.set$set)) {
  marker.gene.set_small=marker.gene.set%>%filter(set==i)
  genes.for.scoring <- list(marker.gene.set_small$gene)
  allseu <- AddModuleScore(object = allseu, features = genes.for.scoring, name = i)
}

### 降维图
DimPlot(
  allseu,reduction = "tsne",group.by = "cluster",pt.size = 1,
  label = T,repel = T,label.size = 6,#label.color = color_cluster[unique(allseu@meta.data$cluster)]
)+scale_color_manual(values = color_cluster)+
  theme(plot.title = element_blank())
ggsave("fig1/tsne_cluster.pdf",width = 20,height = 16,units = "cm")

### 降维图 DIY 1
Idents(allseu)="cluster"
tsne_data=DimPlot(allseu,reduction = "tsne")
tsne_data=tsne_data$data
# head(tsne_data)

tsne_dim=as.data.frame(Embeddings(allseu,"tsne"))
tsne_dim$CB=rownames(tsne_dim)
allseu@meta.data=allseu@meta.data%>%inner_join(tsne_dim,by = "CB")
rownames(allseu@meta.data)=allseu@meta.data$CB
# tsne_data=allseu@meta.data[,c("cluster","tSNE_1","tSNE_2")]
# head(tsne_data)

p1=ggplot()+geom_point(data = tsne_data,mapping = aes(tSNE_1,tSNE_2,color=ident),size=0.5)+
  scale_color_manual(values = color_cluster)+
  theme_classic()

label.df=data.frame(x=c(24,-12,35,50,-40),
                    y=c(31,-27,-30,-3,22),
                    label=c("B cell","CAFs","Epithelial cell","Myeloid","T cell"))

p2=p1+geom_text(data = label.df,mapping = aes(x=x,y=y,label=label,color=label),size=10)+
  scale_x_continuous(limits = c(-49,60),expand = c(0,0))+
  scale_y_continuous(limits = c(-44,44),expand = c(0,0))+
  theme(
    legend.position = "none"
  )

segment.df=data.frame(x=c(-49,-49),
                      xend=c(-24,-49),
                      y=c(-44,-44),
                      yend=c(-44,-19))

p2+geom_segment(data = segment.df,mapping = aes(x=x,xend=xend,y=y,yend=yend),arrow = arrow(length=unit(0.3, "cm")),size=1)+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x.bottom = element_text(hjust = 0.1,size = 18),
    axis.title.y.left = element_text(hjust = 0.1,size = 18)
  )+
  scale_x_continuous(limits = c(-50,60),expand = c(0,0))+
  scale_y_continuous(limits = c(-45,44),expand = c(0,0))

ggsave("./fig1/tsne.cluster.2.pdf",width = 22,height = 18,units = "cm")

### 降维图 DIY 2
ggplot(tsne_data, aes(x=tSNE_1,y=tSNE_2))+
  geom_density_2d_filled(bins=10)+
  scale_fill_manual(values = colorRampPalette(c("#ffffff","#75aadb"))(10))+
  geom_density_2d(colour = "#7b4d4d",alpha=0.8)+
  geom_point(aes(color=ident),alpha=0.5,shape=16)+
  scale_color_manual(values = color_cluster)+
  scale_x_continuous(expand = c(0.01,0))+
  scale_y_continuous(expand = c(0.01,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
ggsave("./fig1/tsne.cluster.3.pdf",width = 32,height = 30,units = "cm")

### 投影图
#知识点：点的顺序有助于显示；在Seurat系列绘图函数中用&设置多图的属性
FeaturePlot(allseu,features = c("Epithelial_cells1","Myeloid_cells1","T_cells1","B_cells1","CAFs1"),
            reduction = "tsne",pt.size = 0.5,ncol = 3,cols = c("#ccccca", "#e61a2e"),order = T)&
  theme_bw()&
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 20)
  )
ggsave("./fig1/marker.gene.set.pdf",width = 21,height = 16,units = "cm")

### 改写函数
source("syFeaturePlot.R")
syFeaturePlot(seu_obj = allseu,
              gene_plot = c("Epithelial_cells1","Myeloid_cells1","T_cells1","B_cells1","CAFs1"),
              group = "cluster",
              highlight = c("Epithelial cell","Myeloid","T cell","B cell","CAFs"),
              trick = "under",
              redim = "tsne",
              color_low = "#ccccca",
              color_high = "#e61a2e",
              filename = "./fig1/geneset1"
)
syFeaturePlot(seu_obj = allseu,
              gene_plot = c("Epithelial_cells1","Myeloid_cells1","T_cells1","B_cells1","CAFs1"),
              group = "cluster",
              highlight = c("Epithelial cell","Myeloid","T cell","B cell","CAFs"),
              trick = "random",
              redim = "tsne",
              color_low = "#ccccca",
              color_high = "#e61a2e",
              filename = "./fig1/geneset2"
)

### 额外的fig 1可能用到的图 ###################################
###cell type 找差异基因
allseu@meta.data$cluster=factor(allseu@meta.data$cluster,levels = sort(unique(allseu@meta.data$cluster)))
Idents(allseu)="cluster"
markers2 <- FindAllMarkers(allseu, logfc.threshold = 0.25, min.pct = 0.1,
                           only.pos = T, test.use = "wilcox")
markers_celltype=markers2
markers_celltype=markers_celltype%>%filter(p_val_adj < 0.01)
markers_celltype=markers_celltype%>%arrange(cluster,desc(avg_log2FC))
markers_celltype=markers_celltype %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
markers_celltype=as.data.frame(markers_celltype)
write.xlsx(markers_celltype,file = "./fig1/allseu.markers_5celltype_log2fc0.25_padj0.01_top200.xlsx",row.names = F,col.names = T)

### 热图展示top20基因
top20 = as.data.frame(markers_celltype %>% group_by(cluster) %>% top_n(20, avg_log2FC))
top20$cluster = factor(top20$cluster,levels = sort(as.character(unique(top20$cluster))))
top20 = top20 %>% arrange(cluster,desc(avg_log2FC))
top20 = as.data.frame(top20)

allseu=ScaleData(allseu,features = union(VariableFeatures(allseu),top20$gene))
DoHeatmap(allseu,features = top20$gene,
          group.by = "cluster",group.bar = TRUE,group.bar.height = 0.03,group.colors = color_cluster,
          label = T,size=4,
          slot = "scale.data", disp.min = -2.5, disp.max = 2.5)
ggsave(filename = "./fig1/5celltypes_top20_heatmap.png",width = 80,height = 45,units = "cm")

### 气泡图展示经典marker基因
library(RColorBrewer)
library(scales)

bubble_gene=read.table("./fig1/bubble_gene.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(bubble_gene)=c("cluster","gene")
bubble_gene$cluster=factor(bubble_gene$cluster,levels = c("B cell","CAFs","Epithelial cell","Myeloid","T cell"))
bubble_gene=bubble_gene%>%arrange(cluster,gene)

DotPlot(allseu,features = bubble_gene$gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  coord_flip()+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
ggsave("./fig1/marker.bubble.pdf",device = "pdf",width = 9.5,height = 11,units = "cm")

### 堆叠小提琴图展示经典marker基因
VlnPlot(allseu,features = bubble_gene$gene,pt.size = 0,stack = T)&
  scale_fill_manual(values = c(brewer.pal(11,"Set3"),brewer.pal(3,"Dark2")))&
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(angle = 45,size = 10,hjust = 0,vjust = 0),
    legend.position = "none"
  )
ggsave("./fig1/marker.Vln.pdf",device = "pdf",width = 18,height = 10,units = "cm")

### 分样本看比例
library(patchwork)

bar.df=allseu@meta.data
text.df=as.data.frame(table(bar.df$sample))
p1=bar.df%>%ggplot(aes(x=sample))+geom_bar(aes(fill=cluster))+
  scale_x_discrete("")+
  scale_y_continuous("cell number",expand = c(0.02,0))+
  scale_fill_manual(values = color_cluster)+
  geom_text(data = text.df,aes(x=Var1,y=Freq,label=Freq),size=4)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )

p2=bar.df%>%ggplot(aes(x=sample))+geom_bar(aes(fill=cluster),position = "fill")+
  scale_x_discrete("")+
  scale_y_continuous("cell ratio",expand = c(0.02,0),labels = scales::label_percent())+
  scale_fill_manual(values = color_cluster)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  )

p1 / p2
ggsave("./fig1/cell_number_and_ratio.bar.pdf",width = 21,height = 18,units = "cm")

### QC结果
VlnPlot(allseu,features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),pt.size = 0,group.by = "cluster")&
  scale_fill_manual(values = color_cluster)&
  theme(
    axis.title.x.bottom = element_blank()
  )
ggsave("./fig1/QC.pdf",width = 20,height = 9,units = "cm")

### 可视化结束