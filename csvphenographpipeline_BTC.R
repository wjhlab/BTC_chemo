rm(list = ls())

#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

metaDataFile = paste0(work,"/Config/metadata.xlsx")
panelDataFile = paste0(work,"/Config/cleanpanel.xlsx")
dataDirectory = paste0(work,"/Data")
PDACcompositionfile = paste0(work,"/Config/PDAC_composition.xlsx")

#libraries 
require(scales);require(readxl);require(plyr);require(dplyr);require(DataEditR); require(promises); require(qgraph)
require(xfun); require(Rphenograph);require(Hmisc); require(ComplexHeatmap); require(pals); require(matrixStats);
require(reshape2); require(ggplot2); require(ggpubr); require(tidyr); require(stringr); require(spatstat); require(gridExtra);
require(rstatix); require(ggrepel); require(tidyverse); require(magrittr)

### load data & functions=====
source(paste0(work,'/functions/clusterfcs.R'))
source(paste0(work,'/functions/do_CI_quantification.R'))
source(paste0(work,'/functions/do_umap.R'))
source(paste0(work,'/functions/plotUmap.R'))
source(paste0(work,'/functions/plot_clustering_heatmap_wrapper2.R'))

## START HERE IF PRIOR DATA WAS SAVED THEN SKIP TO LINE 196====
output<-readRDS("global_data.RDS")
data_full <- data.frame(output[1])
data <- data.matrix(output[2])
data01 <- output[3]
csv_full <- output[4]

## Read-in metadata and clean
ifelse(grepl(metaDataFile,pattern='.xlsx'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
md$file_name <- factor(md$file_name)
md$file_order <- factor(md$file_order)
md$response <- factor(md$response)
md$case <- factor(md$case)
md$sample_id <- factor(md$sample_id)

## Read-in Clean Panel
cleanpanel <- read_xlsx(panelDataFile)
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
tcellmarkers <- cleanpanel$clean_names[cleanpanel$tcell == 1]
myelmarkers <- cleanpanel$clean_names[cleanpanel$myeloid == 1]
stromamarkers <- cleanpanel$clean_names[cleanpanel$stroma == 1]
tumormarkers <- cleanpanel$clean_names[cleanpanel$tumor == 1]

## Set up levels
samplevels=md$sample_id

responselevels=c("S", "P")

caselevels=c("P01", "P03", "P09", "P13")

clusterlevels=c("B",
                "Endoth",
                "Myeloid",
                "Stroma",
                "T",
                "Tumor",
                "UA")

clusterdelevels=unique(data_full$cluster_detailed)

nontumorcells = c("Monocyte",
                  "CAF_I",
                  "CAF_II",
                  "ThM_Cyt",
                  "Endothelial",
                  "CAF_III",
                  "Mac_I",
                  "DC_prolif",
                  "ThM",
                  "Tc_Cyt",
                  "Tc_Eff",
                  "Treg",
                  "Mac_II",
                  "NK",
                  "Tmem",
                  "Monocyte_prolif",
                  "B",
                  "ThM_prolif",
                  "Mac_I_prolif")

immunecells = c("Monocyte",
                "ThM_Cyt",
                "Mac_I",
                "DC_prolif",
                "ThM",
                "Tc_Cyt",
                "Tc_Eff",
                "Treg",
                "Mac_II",
                "NK",
                "Tmem",
                "Monocyte_prolif",
                "B",
                "ThM_prolif",
                "Mac_I_prolif")

stromalcells = c("CAF_I",
                 "CAF_III",
                 "Endothelial",
                 "CAF_II")

CAFs = c("CAF_I",
         "CAF_II",
         "CAF_III")

tumorcells = c("Cancer_NOS",
               "Classical",
               "Classical_prolif",
               "Mixed",
               "Mixed_prolif",
               "Basal")

####START HERE IF NO PRIOR DATA WAS SAVED=====
## Read-in metadata and clean
ifelse(grepl(metaDataFile,pattern='.xlsx'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
md$file_name <- factor(md$file_name)
md$file_order <- factor(md$file_order)
md$response <- factor(md$response)
md$case <- factor(md$case)
md$sample_id <- factor(md$sample_id)

## input image id into metadata
image_id<-c()
for (i in 1:length(md$file_name)){
  tempfile <- read.csv(paste0(dataDirectory,"/",md$file_name[i]))
  df<- as.data.frame(cbind(paste0(md$file_name[i]), unique(tempfile$ImageId)))
  image_id<-rbind(image_id,df)
}
md$ImageId <- image_id$V2[match(image_id$V1,md$file_name)]


## Make sure all files in metadata present in datadirectory
if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])){
  print(paste('ERR: not all filenames in metadata present in data folder - missing',
              md$file_name[!which(data$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),
                                                                                   pattern = '.csv')])],'Subsetting...'))
  md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])),]
}

## Set up levels
samplevels=md$sample_id

responselevels=c("S", "P")

caselevels=c("P01", "P03", "P09", "P13")

## Read csv into csv_raw
csv_raw <- lapply(paste0(dataDirectory,"/",md$file_name),read.csv)
csv_raw_full <- plyr::ldply(csv_raw, rbind)
csv_raw_full$ImageId <- md$sample_id[match(csv_raw_full$ImageId,md$ImageId)]

## clean csv_raw_full dataframe to csv_full containing analysis markers only
cleanpanel <- read_xlsx(panelDataFile)
colnames(csv_raw_full) <- cleanpanel$clean_names
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]
csv_full <- csv_raw_full[,colnames(csv_raw_full) %in% panel]
data_full <- csv_full
data_full$image_id <- md$ImageId[match(csv_full$sample_id,md$sample_id)]
data_full <- data_full[, c(1, ncol(data_full), 2:(ncol(data_full) - 1))]


## sort panels into different categories
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
tcellmarkers <- cleanpanel$clean_names[cleanpanel$tcell == 1]
myelmarkers <- cleanpanel$clean_names[cleanpanel$myeloid == 1]
stromamarkers <- cleanpanel$clean_names[cleanpanel$stroma == 1]
tumormarkers <- cleanpanel$clean_names[cleanpanel$tumor == 1]

## Cluster heatmap for unannotated clusters
data <- data.matrix(csv_full[,-1])
data <- asinh(data[, union(subtype_markers,functional_markers)] / 0.8)

## phenograph clustering of data
rng <- colQuantiles(data, probs = c(0.01, 0.99))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0; data01[data01 > 1] <- 1;data01 <-data01[,union(subtype_markers,functional_markers)]

set.seed(1234)
phenographout<-Rphenograph(data01)
data_full$cluster<-factor(membership(phenographout[[2]]))

#save as RDS file
global_data <- list(data_full, data, data01, csv_full)
saveRDS(global_data, "global_data.RDS")
        
#Start here if you saved data_full with the unannotated clusters====
cluster_mean <- data.frame(data01, cluster = data_full$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_mat<-as.matrix(cluster_mean[,union(subtype_markers,functional_markers)])

rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)

cluster_scaled<-t(scale(t(cluster_mean_mat)))

rownames(cluster_scaled)<-1:nrow(cluster_scaled)


## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_full$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))


pdf("unannotated_clusterheatmap.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="Phenograph Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off()

#cluster heatmap for merged annotations
clusterMergeFile = paste0(work,"/Config/merge.xlsx")
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("B",
                "Endothelial",
                "Myeloid",
                "Stroma",
                "T",
                "Tumor",
                "UA")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_full$cluster, cluster_merging$original_cluster)
data_full$cluster1m <- cluster_merging$new_cluster[mm1]

cluster_mean_merged <- data.frame(data01, cluster = data_full$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])

cluster_scaled_merged<-t(scale(t(cluster_mean_merged_mat)))

rownames(cluster_scaled_merged)<-1:nrow(cluster_scaled_merged)

## Annotation for the merged clusters

if(!is.null(clusterMergeFile)){
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
  annotation_row$Merged <- cluster_merging$new_cluster
  color_clusters2 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Merged))
  names(color_clusters2) <- levels(cluster_merging$new_cluster)
  annotation_colors$Merged <- color_clusters2
}

## Colors for the heatmap

legend_breaks = seq(from = 0, to = 1, by = 0.2)

clusternames<-clusterlevels

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))

names(colorassigned)<-clusternames

rownames(cluster_scaled_merged)<-cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

color_list_byoriginal = colorassigned[match(unique(cluster_merging$new_cluster),names(colorassigned))]

cp<-rowAnnotation(col=color_list,
                  gp = gpar(col = "white", lwd = .5),
                  counts= anno_barplot(
                    as.vector(table(data_full$cluster1m)),
                    gp = gpar(fill=colorassigned),
                    border = F,
                    bar_width = 0.75, 
                    width = unit(2,"cm")))

pdf("clusterheatmap_merged.pdf",width=10,height=4)
Heatmap(cluster_scaled_merged,
        column_title="Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()

## save RDS file again with cluster1m
global_data <- list(data_full, data, data01, csv_full)
saveRDS(global_data, "global_data.RDS")

#Start here if you saved data_full with the annotated clusters====
## Sub-clustering Global Clusters====
## T cell sub clusters
data_immune <- data_full[data_full$cluster1m %in% c("T"),]
data_immune_expr <- data.matrix(asinh(data_immune[,tcellmarkers]/0.8))

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,tcellmarkers]


set.seed(1234)
phenographout_immune<-Rphenograph(data_immune_expr01)
data_immune$cluster<-factor(membership(phenographout_immune[[2]]))
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tcellmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
tcell_data <- list(data_immune, data_immune_expr, data_immune_expr01)
#saveRDS(tcell_data, "tcell_data.RDS")
#readRDS("tcell_data.RDS") -> tcell_data; tcell_data[[1]] -> data_immune

data_immune -> data_immune_t

#Load output of tcell Data
t_data <- readRDS("tcell_data.RDS")
data_immune <- t_data[[1]]
data_immune_t <- data_immune
data_immune_expr <- t_data[[2]]
data_immune_expr01 <- t_data[[3]]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tcellmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCtcells.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph T Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(1000)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

##Myel sub clusters
data_immune <- data_full[data_full$cluster1m %in% c("Myeloid"),]
data_immune_expr <- data.matrix(asinh(data_immune[,myelmarkers]/0.8))

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,myelmarkers]


set.seed(1234)
phenographout_immune<-Rphenograph(data_immune_expr01)
data_immune$cluster<-factor(membership(phenographout_immune[[2]]))
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,myelmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
myelcell_data <- list(data_immune, data_immune_expr, data_immune_expr01)
#saveRDS(myelcell_data, "myel_data.RDS")
#readRDS("myel_data.RDS")->myelcell_data; myelcell_data[[1]]->data_immune

data_immune -> data_immune_myel

#Load output of myel data
myel_data <- readRDS("myel_data.RDS")
data_immune <- myel_data[[1]]
data_immune_myel <- data_immune
data_immune_expr <- myel_data[[2]]
data_immune_expr01 <- myel_data[[3]]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,myelmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCMyel.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Myeloid Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(1000)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

## Tumor sub clusters
data_immune <- data_full[data_full$cluster1m %in% c("Tumor"),]
data_immune_expr <- data.matrix(asinh(data_immune[,tumormarkers]/0.8))

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,tumormarkers]


set.seed(1234)
phenographout_immune<-Rphenograph(data_immune_expr01)
data_immune$cluster<-factor(membership(phenographout_immune[[2]]))
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tumormarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
Tumorcell_data <- list(data_immune, data_immune_expr, data_immune_expr01)
#saveRDS(Tumorcell_data, "Tumor_data.RDS")
#readRDS("Tumor_data.RDS")->Tumorcell_data; Tumorcell_data[[1]]->data_immune

data_immune -> data_immune_Tumor

#Load output of Tumor data
Tumor_data <- readRDS("Tumor_data.RDS")
data_immune <- Tumor_data[[1]]
data_immune_tumor <- data_immune
data_immune_expr <- Tumor_data[[2]]
data_immune_expr01 <- Tumor_data[[3]]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tumormarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCTumor.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Tumor Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(1000)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

## stroma sub clusters
data_immune <- data_full[data_full$cluster1m %in% c("Stroma"),]
data_immune_expr <- data.matrix(asinh(data_immune[,stromamarkers]/0.8))

rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,stromamarkers]


set.seed(1234)
phenographout_immune<-Rphenograph(data_immune_expr01)
data_immune$cluster<-factor(membership(phenographout_immune[[2]]))
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,stromamarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#save as RDS file
stromyelell_data <- list(data_immune, data_immune_expr, data_immune_expr01)
#saveRDS(stromyelell_data, "stroma_data.RDS")
#readRDS("stroma_data.RDS")->stromyelell_data; stromyelell_data[[1]]->data_immune

data_immune -> data_immune_stroma

#Load output of stroma data
stroma_data <- readRDS("stroma_data.RDS")
data_immune <- stroma_data[[1]]
data_immune_stroma <- data_immune
data_immune_expr <- stroma_data[[2]]
data_immune_expr01 <- stroma_data[[3]]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,stromamarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)


## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCstroma.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Stroma Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(1000)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

### annotated t cell heatmap
clusterMergeFileT = paste0(work,"/Config/tcellannotations.xlsx") #create dummy merger numbers prior to annotation
cluster_mergingT <- read_excel(clusterMergeFileT)

clusterlevelsT=c("NK",
                 "ThM_prolif",
                 "ThM_Cyt",
                 "ThM",
                 "Treg",
                 "Tmem",
                 "Tc_Cyt",
                 "Tc_Eff")

colorassignedT<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mergingT$new_cluster)))
clusternamesT<-clusterlevelsT
names(colorassignedT)<-clusternamesT
mm1 <- match(data_immune_t$cluster, cluster_mergingT$original_cluster)
data_immune_t$cluster1m <- cluster_mergingT$new_cluster[mm1]

## matrix for the heatmap
data_immune_expr <- data.matrix(asinh(data_immune_t[,tcellmarkers]/0.8))
rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,tcellmarkers]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune_t$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tcellmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]
rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune_t$cluster1m)),
                       gp = gpar(fill=color_list_byoriginal),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCtcells_merged.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph T Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

### annotated myel heatmap

clusterMergeFileM = paste0(work,"/Config/myelannotations.xlsx") #create dummy merger numbers prior to annotation
cluster_mergingM <- read_excel(clusterMergeFileM)

clusterlevelsM=c("Monocyte",
                 "Mac_I",
                 "DC_prolif",
                 "Mac_II",
                 "Monocyte_prolif",
                 "Mac_I_prolif")

colorassignedM<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mergingM$new_cluster)))
clusternamesM<-clusterlevelsM
names(colorassignedM)<-clusternamesM
mm1 <- match(data_immune_myel$cluster, cluster_mergingM$original_cluster)
data_immune_myel$cluster1m <- cluster_mergingM$new_cluster[mm1]

## matrix for the heatmap
data_immune_expr <- data.matrix(asinh(data_immune_myel[,myelmarkers]/0.8))
rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,myelmarkers]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune_myel$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,myelmarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]
rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune_myel$cluster1m)),
                       gp = gpar(fill=color_list_byoriginal),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCmyel_merged.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Myel Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off()

### annotated tumor heatmap

clusterMergeFileE = paste0(work,"/Config/tumorannotations.xlsx") #create dummy merger numbers prior to annotation
cluster_mergingE <- read_excel(clusterMergeFileE)

clusterlevelsE=c("Basal",
                 "Classical_prolif",
                 "Classical",
                 "Mixed",
                 "Cancer_NOS",
                 "Mixed_prolif")

colorassignedE<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mergingE$new_cluster)))
clusternamesE<-clusterlevelsE
names(colorassignedE)<-clusternamesE
mm1 <- match(data_immune_tumor$cluster, cluster_mergingE$original_cluster)
data_immune_tumor$cluster1m <- cluster_mergingE$new_cluster[mm1]

## matrix for the heatmap
data_immune_expr <- data.matrix(asinh(data_immune_tumor[,tumormarkers]/0.8))
rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,tumormarkers]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune_tumor$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,tumormarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]
rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune_tumor$cluster1m)),
                       gp = gpar(fill=color_list_byoriginal),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCtumor_merged.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Tumor Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

### annotated stroma heatmap

clusterMergeFileS = paste0(work,"/Config/stromaannotations.xlsx") #create dummy merger numbers prior to annotation
cluster_mergingS <- read_excel(clusterMergeFileS)

clusterlevelsS=c("CAF_I",
                 "CAF_II",
                 "CAF_III",
                 "Endothelial")

colorassignedS<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mergingS$new_cluster)))
clusternamesS<-clusterlevelsS
names(colorassignedS)<-clusternamesS
mm1 <- match(data_immune_stroma$cluster, cluster_mergingS$original_cluster)
data_immune_stroma$cluster1m <- cluster_mergingS$new_cluster[mm1]

## matrix for the heatmap
data_immune_expr <- data.matrix(asinh(data_immune_stroma[,stromamarkers]/0.8))
rng <- colQuantiles(data_immune_expr, probs = c(0.05, 0.95))
data_immune_expr01 <- t((t(data_immune_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_immune_expr01[data_immune_expr01 < 0] <- 0; data_immune_expr01[data_immune_expr01 > 1] <- 1;data_immune_expr01 <-data_immune_expr01[,stromamarkers]
cluster_mean <- data.frame(data_immune_expr01, cluster = data_immune_stroma$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,stromamarkers])
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-1:nrow(cluster_scaled)
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassignedS)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassignedS))]
rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_immune_stroma$cluster1m)),
                       gp = gpar(fill=color_list_byoriginal),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCstroma_merged.pdf",width=10,height=8)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Stroma Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

### merging into one dataframe

data_immune_t$full_id <- paste(data_immune_t$sample_id, data_immune_t$CellId,sep="_")
data_immune_myel$full_id <- paste(data_immune_myel$sample_id, data_immune_myel$CellId,sep="_")
data_immune_stroma$full_id <- paste(data_immune_stroma$sample_id, data_immune_stroma$CellId,sep="_")
data_immune_tumor$full_id <- paste(data_immune_tumor$sample_id, data_immune_tumor$CellId,sep="_")
data_full$full_id <- paste(data_full$sample_id, data_full$CellId,sep="_")

data_full$cluster_detailed <- data_full$cluster1m

data_full$cluster_detailed[match(data_immune_t$full_id, data_full$full_id)] <- data_immune_t$cluster1m
data_full$cluster_detailed[match(data_immune_myel$full_id, data_full$full_id)] <- data_immune_myel$cluster1m
data_full$cluster_detailed[match(data_immune_tumor$full_id, data_full$full_id)] <- data_immune_tumor$cluster1m
data_full$cluster_detailed[match(data_immune_stroma$full_id, data_full$full_id)] <- data_immune_stroma$cluster1m

#Cluster Heatmap by 
data <- data.matrix(data_full[,union(subtype_markers,functional_markers)])
data <- asinh(data / 0.8)

#rescale
rng <- colQuantiles(data, probs = c(0.05, 0.95))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0; data01[data01 > 1] <- 1;data01 <-data01[,union(subtype_markers,functional_markers)]
cluster_mean <- data.frame(data01, cluster = data_full$cluster_detailed, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-cluster_mean$cluster

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_full$cluster_detailed)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("clusterheatmap_BTCfull.pdf",width=10,height=10)
Heatmap(cluster_scaled,
        column_title="BTC Phenograph Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

#Custom order
marker_order <- c("Collagen", "VIM", "SMA", "CD45", "CD20", "CD3", "CD4", "CD8", "CD57", "CD45RO", "CD45RA", "CXCL12",
                  "CD16", "S100A4", "CD68", "HLADR", "CD74", "CD163", "PDGFRA", "FAP", "NCAD", "PDPN",
                  "CD31", "FOXP3", "CD105", "IL8", "IL6", "PDL1", "PD1", "GZMB", "CK", "ECAD", "TFFI",
                  "GATA6", "KRT17", "S100A2", "DCSIGN")

cluster_order <- c("B", "Tc_Cyt", "Tc_Eff", "ThM", "ThM_Cyt", "ThM_prolif", "Tmem", "Treg", "NK", "DC_prolif",
                   "Mac_I", "Mac_I_prolif", "Mac_II", "Monocyte", "Monocyte_prolif", "Endothelial", "CAF_I",
                   "CAF_II", "CAF_III", "UA", "Epith", "Classical", "Classical_prolif", "Mixed", "Mixed_prolif",
                   "Basal")

#Prepare scaled mean matrix for clusters and markers
data <- data.matrix(data_full[, union(subtype_markers, functional_markers)])
data <- asinh(data / 0.8)

#Rescale
rng <- colQuantiles(data, probs = c(0.05, 0.95))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0; data01[data01 > 1] <- 1;data01 <-data01[,union(subtype_markers,functional_markers)]
cluster_mean <- data.frame(data01, cluster = data_full$cluster_detailed, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-as.matrix(cluster_mean[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)
cluster_scaled<-t(scale(t(cluster_mean_mat)))
rownames(cluster_scaled)<-cluster_mean$cluster

#Reorder matrix by your markers (columns) and clusters (rows)
present_markers <- intersect(marker_order, colnames(cluster_scaled))
present_clusters <- intersect(cluster_order, rownames(cluster_scaled))
cluster_scaled <- cluster_scaled[present_clusters, present_markers, drop = FALSE]

#Annotation and colors
annotation_row <- data.frame(Cluster = factor(present_clusters, levels = cluster_order))
rownames(annotation_row) <- present_clusters

color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(length(present_clusters))
names(color_clusters1) <- present_clusters
annotation_colors <- list(Cluster = color_clusters1)
colorassigned <- color_clusters1

rAbar <- rowAnnotation(
  clusters = factor(present_clusters, levels = cluster_order),
  col = list(clusters = colorassigned),
  gp = gpar(col = "white", lwd = .5),
  counts = anno_barplot(
    as.vector(table(factor(data_full$cluster_detailed, levels = cluster_order))),
    gp = gpar(fill = colorassigned),
    border = F,
    bar_width = 0.75, 
    width = unit(2, "cm"))
)

#Plot to PDF (markers=x axis, clusters=y axis)
pdf("clusterheatmap_BTCfull_ordered.pdf", width = 10, height = 10)
Heatmap(
  cluster_scaled,
  column_title = "BTC Phenograph Clusters",
  name = "scaled",
  col = rev(brewer.rdbu(100)),
  cluster_columns = FALSE,           # Keep marker order
  cluster_rows = FALSE,              # Keep cluster order
  border = NA,
  rect_gp = gpar(col = "white", lwd = .5),
  right_annotation = rAbar,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    at = seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))
  ),
  width = ncol(cluster_scaled) * unit(3, "mm"), 
  height = nrow(cluster_scaled) * unit(4, "mm")
)
dev.off()


####Density and Abundance plots======
#totalcounts
totcounts_table<-table(data_full$cluster1m, data_full$sample_id)
totprops_table <- t(t(totcounts_table) / colSums(totcounts_table)) * 100
totcounts <- as.data.frame.matrix(totcounts_table)
totprops <-as.data.frame.matrix(totprops_table)
totcell_table <- table(data_full$sample_id)

#subcounts
counts_table<-table(data_full$cluster_detailed, data_full$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <-as.data.frame.matrix(props_table)
cell_table <- table(data_full$sample_id)

ggdfccounts <- melt(
  data.frame(cluster = rownames(counts),
             counts,
             check.names = FALSE),
  id.vars = "cluster", value.name = "counts", variable.name = "sample_id"
)

ggdfccounts$sample_id <- factor(ggdfccounts$sample_id, levels=samplevels)
ggdfccounts$response <- factor(md$response[match(ggdfccounts$sample_id,md$sample_id)], levels=responselevels)
ggdfccounts$case <- factor(md$case[match(ggdfccounts$sample_id,md$sample_id)], levels=caselevels)

caf_clusters <- c("CAF_I", "CAF_II", "CAF_III")
caf_counts_summary <- subset(ggdfccounts, cluster %in% caf_clusters & !is.na(response))
caf_counts_summary <- aggregate(counts ~ cluster + response,
                                  data = caf_counts_summary,
                                  FUN = sum)

caf_counts_by_sample <- subset(ggdfccounts,
                                 cluster %in% caf_clusters & !is.na(response))
caf_counts_by_sample <- aggregate(counts ~ cluster + sample_id + response,
                                    data = caf_counts_by_sample, FUN = sum)

caf_counts_by_case <- subset(ggdfccounts,
                               cluster %in% caf_clusters & !is.na(response))
caf_counts_by_case <- aggregate(counts ~ cluster + case + response,
                                  data = caf_counts_by_case, FUN = sum)


write.csv(caf_counts_by_sample,'caf_counts_by_sample.csv')
write.csv(caf_counts_by_case,'caf_counts_by_case.csv')

#Densities
areas <- read_xlsx(paste0(work,'/Config/area.xlsx'))
densities <- t(t(counts)/areas$TotalArea)

write.csv(counts,'counts.csv')
write.csv(props,'props.csv')
write.csv(densities, 'densities.csv')

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities), densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$response <- factor(md$response[match(ggdfd$sample_id,md$sample_id)], levels=responselevels)
ggdfd$case <- factor(md$case[match(ggdfd$sample_id,md$sample_id)], levels=caselevels)

#Abundance dataframe
ggdfp <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
              id.vars = "cluster", value.name = "proportion", 
              variable.name = "sample_id")
#Non-ordered
ggdfp$sample_id <- factor(ggdfp$sample_id, levels=samplevels)
#Ordered
ggdfp$sample_id <- factor(ggdfp$sample_id, levels= c( "H1_R1", "H1_R2", "H1_N1", "H1_N2", 
                                                      "H3_R1", "H3_R2", "H3_N1", "H3_N2",
                                                      "H9_R1", "H9_R2", "H9_N1", "H9_N2",
                                                      "H13_R1", "H13_R2", "H13_N1", "H13_N2"))
ggdfp$response <- factor(md$response[match(ggdfp$sample_id,md$sample_id)], levels=responselevels)
ggdfp$case <- factor(md$case[match(ggdfp$sample_id,md$sample_id)], levels=caselevels)

# Define statistical tests with their corresponding labels
statistical_tests <- list("t.test" = "T-Test", "wilcox.test" = "Wilcoxon-Test")

# Start PDF output
pdf('Density_and_Abundance_box_response.pdf', width = 9, height = 12)

# Loop through each statistical test and create a plot for each (Density first, then Abundance)
for (test in names(statistical_tests)) {
  
  # Density Plot by response
  ggp_density <- ggplot(ggdfd, aes(x = response, y = densities, fill = response)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.5) +
    geom_jitter(width = 0.05) +
    scale_shape_manual(values = c(1:10, 1:18, 1:8)) +
    facet_wrap(~cluster, ncol = 6, scales = "free") +
    ylab("Cells per mm²") +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.x = element_blank(),
      axis.line.x = element_line(linewidth = 0.25, color = "black"),
      axis.line.y = element_line(linewidth = 0.25, color = "black"),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 8, color = "black"),
      strip.background = element_rect(fill = NA),
      strip.text = element_text(size = 8, color = "black"),
      panel.background = element_rect(fill = "white"),
      legend.title = element_blank(),
      legend.key.size = unit(2, 'lines'),
      legend.text = element_text(size = 8),
      legend.key = element_rect(fill = "white"),
      strip.text.x = element_text(size = 8),
      panel.spacing = unit(1, "lines")
    ) +
    stat_summary(
      fun = "mean",
      geom = "crossbar",
      width = 0.4,
      linewidth = 0.1,
      position = position_dodge(0.9),
      colour = "black"
    ) +
    stat_compare_means(
      aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
      method = test,
      label.x = 1.4,
      size = 3,
      color = "black",
      hide.ns = TRUE
    ) +
    ggtitle(paste("Cells by response -", statistical_tests[[test]], "Test (Density)"))
  
  # Print Density plot to PDF
  print(ggp_density)
  
  # Abundance Plot by response
  ggp_abundance <- ggplot(ggdfp, aes(x = response, y = proportion, fill = response)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.5) +
    geom_jitter(width = 0.05) +
    scale_shape_manual(values = c(1:10, 1:18, 1:8)) +
    facet_wrap(~cluster, ncol = 6, scales = "free") +
    ylab("% of Cells") +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.x = element_blank(),
      axis.line.x = element_line(linewidth = 0.25, color = "black"),
      axis.line.y = element_line(linewidth = 0.25, color = "black"),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 8, color = "black"),
      strip.background = element_rect(fill = NA),
      strip.text = element_text(size = 8, color = "black"),
      panel.background = element_rect(fill = "white"),
      legend.title = element_blank(),
      legend.key.size = unit(2, 'lines'),
      legend.text = element_text(size = 8),
      legend.key = element_rect(fill = "white"),
      strip.text.x = element_text(size = 8),
      panel.spacing = unit(1, "lines")
    ) +
    stat_summary(
      fun = "mean",
      geom = "crossbar",
      width = 0.4,
      linewidth = 0.1,
      position = position_dodge(0.9),
      colour = "black"
    ) +
    stat_compare_means(
      aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
      method = test,
      label.x = 1.4,
      size = 3,
      color = "black",
      hide.ns = TRUE
    ) +
    ggtitle(paste("Cells by response -", statistical_tests[[test]], "Test (Abundance)"))
  
  # Print Abundance plot to PDF
  print(ggp_abundance)
}

# Close the PDF device
dev.off()

#Functional Plots====
supmarkers <- c("Collagen", "HLADR", "CD105", "IL6", "VIM")

## Functional Suppressive markers P vs S
exprtco <- data.frame(data_full[, supmarkers], sample_id = data_full$sample_id, cluster = data_full$cluster_detailed)
exprtcoe <- exprtco %>% filter(cluster %in% c(unique(ggdfd$cluster)))
exprtcoe$response <- factor(md$response[match(exprtcoe$sample_id,md$sample_id)], levels=responselevels)
ggdf2<-melt(exprtcoe, id.var=c("cluster","response","sample_id"))
ggdf2$cluster <- factor(ggdf2$cluster, levels=unique(ggdfd$cluster))
ggdf2$sample_id <- factor(ggdf2$sample_id, levels=samplevels)
ggdf2$response <- factor(md$response[match(ggdf2$sample_id,md$sample_id)], levels=responselevels)
ggdf2$case <- factor(md$case[match(ggdf2$sample_id,md$sample_id)], levels=caselevels)

pdf("plot_functional_vio_suppressive_markers_t_test.pdf",width=6,height=7)
for(i in 1:length(supmarkers)){
  ggp <- ggplot(ggdf2[ggdf2$variable==supmarkers[i],], aes(x=response, y=value, fill=response))+
    ggtitle(supmarkers[i])+
    facet_wrap(~cluster, ncol=6, scales='free')+
    geom_violin()+
    stat_summary(fun = "mean", geom = "crossbar", width = 0.4, size = 0.1, position = position_dodge(0.9), colour = "black") +
    stat_compare_means(method = "t.test", 
                       aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))), 
                       label.x = 1.4,  # Adjust x position if needed
                       label.y.npc = .95,  # Adjust y position if needed
                       size = 2,
                       color = "black",
                       hide.ns = FALSE) +
    #geom_boxplot(fill = NA , outlier.shape = NA)+
    theme(axis.text.x = element_text(size=8, color="black", angle=45, hjust=1),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=0.25),
          axis.line = element_line(size=0.25),
          axis.text = element_text(color="black"),
          axis.title.y = element_text(size=8, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=8),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(1, 'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white"),
          plot.title = element_text(size=8))
  print(ggp)}
dev.off()

## Stacked Bar ====
#list out the cell types and create a legends data frame
allcelltypes<-clusterdelevels
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"Cancer_NOS")]<-"Cancer_NOS"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Classical")]<-"Classical"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Classical_prolif")]<-"Classical"
legendctype$maintype[str_detect(legendctype$allcelltypes,"B")]<-"Lymphoid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Basal")]<-"Basal"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mixed")]<-"Mixed"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mixed_prolif")]<-"Mixed"
legendctype$maintype[str_detect(legendctype$allcelltypes,"ThM_Cyt")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"ThM_prolif")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"ThM")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tc_Cyt")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tc_Eff")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Treg")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tmem")]<-"T"
legendctype$maintype[str_detect(legendctype$allcelltypes,"NK")]<-"Lymphoid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Endothelial")]<-"Endothelial"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CAF_I")]<-"CAF_I"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CAF_II")]<-"CAF_II"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CAF_III")]<-"CAF_III"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DC_prolif")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Monocyte")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Monocyte_prolif")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac_I")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac_I_prolif")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac_II")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"UA")]<-"Other"

#Create Custom Cluster Order
ggdfp$cluster <- factor(ggdfp$cluster, levels =  c("UA", "B", "Tc_Cyt", "Tc_Eff", "ThM", "ThM_Cyt", "ThM_prolif", "Tmem", "Treg", "NK", "DC_prolif",
                                                   "Mac_I", "Mac_I_prolif", "Mac_II", "Monocyte", "Monocyte_prolif", "Endothelial", "CAF_I",
                                                   "CAF_II", "CAF_III", "Cancer_NOS", "Classical", "Classical_prolif", "Mixed", "Mixed_prolif",
                                                   "Basal"))

colorassigned <- c("#0030F5",
                   "#1857C6",
                   "#1C719A",
                   "#338571",
                   "#409646",
                   "#50A51F",
                   "#73B112",
                   "#9BBA16",
                   "#BEC31D",
                   "#E2CA22",
                   "#F7C423",
                   "#F9B120",
                   "#F79C1B",
                   "#F48517",
                   "#F16E16",
                   "#F05B23",
                   "#F65C4D",
                   "#FD6B84",
                   "#FF7EBD",
                   "#FD92FA",
                   "#E5E5E5",
                   "#CCCCCC",
                   "#B2B2B2",
                   "#999999",
                   "#5C5C5C",
                   "#212121")

names(colorassigned) <-  levels(ggdfp$cluster)

pdf("stacked_bar_custom_cluster_order.pdf", width = 5, height = 2)

ggplot(ggdfp, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.4),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.25),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 6, 22, 6)
  ) +
  ylab("Proportion") +
  scale_y_continuous(
    expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = colorassigned) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size=4))) +
  ggtitle("Cell Composition by ROI")

dev.off()

# --- Annotation strip data ---
annotdf <- unique(ggdfp[, c("sample_id", "response", "case")])

# Color palettes
response_palette <- c("#55FF55", "#FF5555")
response_colors <- setNames(response_palette[1:length(responselevels)], responselevels)

case_palette <- RColorBrewer::brewer.pal(max(3, length(caselevels)), "Set2")
case_colors <- setNames(case_palette[1:length(caselevels)], caselevels)

ggp <- ggplot(ggdfp, aes(x = sample_id, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.4),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.25),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 6, 22, 6)
  ) +
  ylab("Proportion") +
  scale_y_continuous(
    expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = colorassigned) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size=4))) +
  ggtitle("Cell Composition by ROI")

# --- Response annotation strip ---
response_tile <- ggplot(annotdf, aes(x = sample_id, y = "Response", fill = response)) +
  geom_tile() +
  scale_fill_manual(values = response_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(-6, 0, -6, 0),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 7)
  )

# --- Case annotation strip ---
case_tile <- ggplot(annotdf, aes(x = sample_id, y = "Case", fill = case)) +
  geom_tile() +
  scale_fill_manual(values = case_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(-6, 0, -6, 0),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 7)
  )

# --- Response legend plot ---
response_legend <- ggplot(data.frame(response = responselevels), aes(x = response, fill = response)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = response_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(title = "Response", override.aes = list(size = 5)))

# --- Case legend plot ---
case_legend <- ggplot(data.frame(case = caselevels), aes(x = case, fill = case)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = case_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(title = "Case", override.aes = list(size = 5)))

# --- Combine and plot ---
pdf("stacked_bar_custom_cluster_order_response_case.pdf", width = 5, height = 4)
(ggp +
    response_tile +
    case_tile +
    response_legend +
    case_legend +
    plot_layout(ncol = 1, heights = c(8, 0.47, 0.47, 0.5, 0.5)))
dev.off()

#Non-epithelial
#Which clusters are epithelial?
epithelial_types <- c("Cancer_NOS", "Classical", "Classical_prolif", "Mixed", "Mixed_prolif", "Basal")
ggdfp_noepi <- ggdfp[!ggdfp$cluster %in% epithelial_types, ]

ggdfp_noepi <- ggdfp_noepi %>%
  group_by(sample_id) %>%
  mutate(prop = 100 * proportion / sum(proportion)) %>%
  ungroup()

#ROI/sample order
samplevels <- c(
  "H1_R1", "H1_R2", "H1_N1", "H1_N2",
  "H3_R1", "H3_R2", "H3_N1", "H3_N2",
  "H9_R1", "H9_R2", "H9_N1", "H9_N2",
  "H13_R1", "H13_R2", "H13_N1", "H13_N2"
)


# 9. Color palette
colorassigned <- c(
  UA              = "#0030F5",
  B               = "#1857C6",
  Tc_Cyt          = "#1C719A",
  Tc_Eff          = "#338571",
  ThM             = "#409646",
  ThM_Cyt         = "#50A51F",
  ThM_prolif      = "#73B112",
  Tmem            = "#9BBA16",
  Treg            = "#BEC31D",
  NK              = "#E2CA22",
  DC_prolif       = "#F7C423",
  Mac_I           = "#F9B120",
  Mac_I_prolif    = "#F79C1B",
  Mac_II          = "#F48517",
  Monocyte        = "#F16E16",
  Monocyte_prolif = "#F05B23",
  Endothelial     = "#F65C4D",
  CAF_I           = "#FD6B84",
  CAF_II          = "#FF7EBD",
  CAF_III         = "#FD92FA"
)

custom_cluster_order <- c("UA", "B", "Tc_Cyt", "Tc_Eff", "ThM", "ThM_Cyt", "ThM_prolif", "Tmem", "Treg", "NK", "DC_prolif",
                          "Mac_I", "Mac_I_prolif", "Mac_II", "Monocyte", "Monocyte_prolif", "Endothelial", "CAF_I",
                          "CAF_II", "CAF_III")

names(colorassigned) <- custom_cluster_order

# 10. Plot
pdf("stacked_bar_custom_cluster_order_nonepi.pdf", width = 5, height = 2)

ggplot(ggdfp_noepi, aes(x = sample_id, y = prop, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.4),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.25),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 6, 22, 6)
  ) +
  ylab("Proportion (non-epithelial)") +
  scale_y_continuous(
    expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size=4))) +
  scale_fill_manual(values = colorassigned) +
  ggtitle("Non-Epithelial Cell Composition by ROI")

dev.off()


#Only Epithelial
#Which clusters are epithelial?
epithelial_types <- c("Cancer_NOS", "Classical", "Classical_prolif", "Mixed", "Mixed_prolif", "Basal")

nonepithelial_types <- c("UA", "B", "Tc_Cyt", "Tc_Eff", "ThM", "ThM_Cyt", "ThM_prolif", "Tmem", "Treg", "NK", "DC_prolif",
                         "Mac_I", "Mac_I_prolif", "Mac_II", "Monocyte", "Monocyte_prolif", "Endothelial", "CAF_I",
                         "CAF_II", "CAF_III", "Epith")

#ROI/sample order
samplevels <- c(
  "H1_R1", "H1_R2", "H1_N1", "H1_N2",
  "H3_R1", "H3_R2", "H3_N1", "H3_N2",
  "H9_R1", "H9_R2", "H9_N1", "H9_N2",
  "H13_R1", "H13_R2", "H13_N1", "H13_N2"
)

ggdfp_epi <- ggdfp[!ggdfp$cluster %in% nonepithelial_types,]

ggdfp_epi <- ggdfp_epi %>%
  group_by(sample_id) %>%
  mutate(prop = 100 * proportion / sum(proportion)) %>%
  ungroup()

ggdfp_epi$cluster <- factor(ggdfp_epi$cluster, levels=custom_cluster_order)
ggdfp_epi$sample_id <- factor(ggdfp_epi$sample_id, levels=samplevels)



# 9. Color palette
colorassigned <- c(
  Cancer_NOS        = "#E5E5E5",
  Classical         = "#CCCCCC",
  Classical_prolif  = "#B2B2B2",
  Mixed             = "#999999",
  Mixed_prolif      = "#5C5C5C",
  Basal             = "#212121"
)

custom_cluster_order <- c("Cancer_NOS", "Classical", "Classical_prolif", "Mixed", "Mixed_prolif", "Basal")

names(colorassigned) <- custom_cluster_order

# 10. Plot
pdf("stacked_bar_custom_cluster_order_epi.pdf", width = 5, height = 2)

ggplot(ggdfp_epi, aes(x = sample_id, y = prop, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 0.4),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.ticks.y = element_line(size = 0.25),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 6, 22, 6)
  ) +
  ylab("Proportion (Epithelial)") +
  scale_y_continuous(
    expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size=4))) +
  scale_fill_manual(values = colorassigned) +
  ggtitle("Epithelial Cell Composition by ROI")

dev.off()

#### Mixed Cluster - Persistent PDAC====
## Read-in PDAC composition and clean
ifelse(grepl(PDACcompositionfile,pattern='.xlsx'),comp <- read_excel(PDACcompositionfile),comp <- read.csv(PDACcompositionfile,header = TRUE))#must be in xl format or csv

#Make the composition a data frame to pull information from
compdf <- as.data.frame(comp)

#Density dataframe
ggdfd <- melt(data.frame(cluster = rownames(densities), densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$response <- factor(md$response[match(ggdfd$sample_id,md$sample_id)], levels=responselevels)
ggdfd$case <- factor(md$case[match(ggdfd$sample_id,md$sample_id)], levels=caselevels)
ggdfd$sensitive_PDAC <- factor(compdf$sensitive_PDAC[match(ggdfd$sample_id,compdf$sample_id)])
ggdfd$persistent_PDAC <- factor(compdf$persistent_PDAC[match(ggdfd$sample_id,compdf$sample_id)])

# Filter the data to include only the "Mixed" cluster and combine Mixed and Mixed_prolif clusters
ggdfdmixedonly <- ggdfd %>%
  mutate(
    cluster_combined = case_when(
      cluster %in% c("Mixed", "Mixed_prolif") ~ "Mixed_Combined",  # Combine Mixed and Mixed_prolif
      TRUE ~ as.character(cluster)  # Keep other clusters as they are
    )
  ) %>%
  group_by(sample_id, case, cluster_combined) %>%  # Group by sample_id and the new combined cluster
  summarise(
    densities = sum(densities, na.rm = TRUE),  # Sum densities for each sample_id
    persistent_PDAC = first(persistent_PDAC),  # Keep the first (same) value for responsive_PDAC in each group
    .groups = "drop"  # Ungroup after summarizing
  ) %>%
  filter(cluster_combined == "Mixed_Combined") %>%
  mutate(
    case = factor(case, levels = caselevels),
    sample_id = factor(sample_id, levels = samplevels),  # Ensure sample_id follows samplevels order
    persistent_PDAC = as.numeric(as.character(persistent_PDAC)),  # Convert to numeric
    densities = as.numeric(as.character(densities))  # Convert to numeric
  )

# Check if the columns are numeric
if (!is.numeric(ggdfdmixedonly$persistent_PDAC) || !is.numeric(ggdfdmixedonly$densities)) {
  stop("Both persistent_PDAC and densities must be numeric.")
}

# Correlation Test for Persistent PDAC in Mixed cluster
cor_test_mixed <- cor.test(ggdfdmixedonly$persistent_PDAC, ggdfdmixedonly$densities, method = "pearson")
cor_value_mixed <- round(cor_test_mixed$estimate, 2)
p_value_mixed <- cor_test_mixed$p.value
p_value_text_mixed <- ifelse(p_value_mixed < 0.001, "< 0.001", paste0("= ", signif(p_value_mixed, 2)))

# Make sure to define casecolors as before
casecolors <- c(
  "P01" = "#1b9e77",
  "P03" = "#d95f02",
  "P09" = "#7570b3",
  "P13" = "#e7298a"
)

ggp_dot_mixed <- ggplot(
  ggdfdmixedonly,
  aes(x = persistent_PDAC, y = densities, color = case)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.5) +
  xlab("% of Persistent PDAC") +
  ylab("# of Mixed Cells per mm²") +
  annotate(
    "text",
    x = max(ggdfdmixedonly$persistent_PDAC, na.rm = TRUE) * 0.6,
    y = max(ggdfdmixedonly$densities, na.rm = TRUE) * 0.9,
    label = paste0("R = ", cor_value_mixed, "\np ", p_value_text_mixed),
    size = 5,
    color = "black"
  ) +
  scale_color_manual(
    values = casecolors,
    breaks = names(casecolors),
    labels = names(casecolors)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(ggdfdmixedonly$persistent_PDAC, na.rm = TRUE), by = 5),
    limits = c(0, max(ggdfdmixedonly$persistent_PDAC, na.rm = TRUE) * 1.1)
  ) +
  scale_y_continuous(
    breaks = seq(0, max(ggdfdmixedonly$densities, na.rm = TRUE), by = 100),
    limits = c(0, max(ggdfdmixedonly$densities, na.rm = TRUE) * 1.1)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),           # REMOVE ALL GRIDLINES!
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(size = 0.3, color = "black"),
    axis.text.x = element_text(size = 8, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, hjust = 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("# Mixed Cells vs % Persistent PDAC")

# Save to file
pdf("Scaled_Mixed_Combined_cells_nondotplot_with_pval.pdf", width = 6, height = 5)
print(ggp_dot_mixed)
dev.off()

#Average Spatial Influence of each cell type====
#Create combined spatial influence dataframe

## Gaussian Kernel Spatial Analysis
dfgk<-data.frame(data_full[,c("X_coord","Y_coord")],
                 #cell_cluster=output$cell_clustering,
                 celltype=data_full$cluster_detailed,
                 sample_id=data_full$sample_id)

colnames(dfgk) <- c("X_position", "Y_position", "celltype", "sample_id")

dfgk$response <- factor(md$response[match(dfgk$sample_id,md$sample_id)], levels=responselevels)

dfgk$sample_id<-factor(data_full$sample_id, levels=unique(data_full$sample_id))

dim(dfgk) # 125091 cells x 5 (x,y coord, celltype, sample_id, response)


#saveRDS(dfgk, "dfgk_data.RDS")

## Spatial weight calculations for each core 
#run if you have no prior data
sample_id <- unique(dfgk$sample_id)
spatwt <- do_CI_quantification(expr = dfgk,  # should contain sample_id, X, Y position, celltype in the df
                               sample_id = sample_id,
                               kernels ="gaussian",
                               clusterdelevels = clusterdelevels,
                               sigma = 10)
spatwt$sample_id <- factor(spatwt$sample_id, levels=samplevels)
spatwt$response <- factor(md$response[match(spatwt$sample_id,md$sample_id)], levels=responselevels)
spatwt_P <- as.data.frame(spatwt)  %>% filter(response != "S")
spatwt_S <- as.data.frame(spatwt) %>% filter(response != "P")

# Save spatwt 
#saveRDS(spatwt, "spatwt.rds")

# Load spatwt 
#Start Here if you have prior data
spatwt <- readRDS("spatwt.RDS")

#TME analysis
spatwt_df <- as.data.frame(spatwt[1:(length(spatwt)-1)]) # remove sample_id
spatwt_df <- spatwt[, clusterdelevels]


# TME of interest analysis
# Basal cells
spatwt_df_basal<- spatwt_df[grepl("Basal", rownames(spatwt_df)), ]
somClustering <- clusterfcs(fcs=as.matrix(spatwt_df_basal), 
                            cluster_by=colnames(spatwt_df), 
                            numclusters = 10, 
                            scaleoption = F)

# Extract sample_id while ensuring all are retained
spatwt_basal <- spatwt[grepl("Basal", rownames(spatwt)), , drop = FALSE]

# Classical cells
spatwt_df_classical<- spatwt_df[grepl("Classical", rownames(spatwt_df)), ]
somClustering <- clusterfcs(fcs=as.matrix(spatwt_df_classical), 
                            cluster_by=colnames(spatwt_df), 
                            numclusters = 10, 
                            scaleoption = F)

# Extract sample_id while ensuring all are retained
spatwt_classical <- spatwt[grepl("Classical", rownames(spatwt)), , drop = FALSE]

# Mixed cells
spatwt_df_mixed<- spatwt_df[grepl("Mixed", rownames(spatwt_df)), ]
somClustering <- clusterfcs(fcs=as.matrix(spatwt_df_mixed), 
                            cluster_by=colnames(spatwt_df), 
                            numclusters = 10, 
                            scaleoption = F)

# Extract sample_id while ensuring all are retained
spatwt_mixed <- spatwt[grepl("Mixed", rownames(spatwt)), , drop = FALSE]

# Cancer_NOS cells
spatwt_df_cancer_NOS<- spatwt_df[grepl("Cancer_NOS", rownames(spatwt_df)), ]
somClustering <- clusterfcs(fcs=as.matrix(spatwt_df_cancer_NOS), 
                            cluster_by=colnames(spatwt_df), 
                            numclusters = 10, 
                            scaleoption = F)

# Extract sample_id while ensuring all are retained
spatwt_cancer_NOS <- spatwt[grepl("Cancer_NOS", rownames(spatwt)), , drop = FALSE]

#Create dataframes for barplots
spatwt_basal_df <- spatwt_basal[, nontumorcells]
spatwt_classical_df <- spatwt_classical[, nontumorcells]
spatwt_mixed_df <- spatwt_mixed[, nontumorcells]
spatwt_cancer_NOS_df <- spatwt_cancer_NOS[, nontumorcells]

spatwt_basal_df$Epithelial <- "Basal"
spatwt_classical_df$Epithelial <- "Classical"
spatwt_mixed_df$Epithelial <- "Mixed"
spatwt_cancer_NOS_df$Epithelial <- "Cancer_NOS"

spatwt_all <- bind_rows(spatwt_basal_df, spatwt_classical_df, spatwt_mixed_df, spatwt_df_cancer_NOS)

#Get mean Spatwt Influence 
long_all <- spatwt_all %>%
  pivot_longer(cols = -Epithelial, names_to = "CellType", values_to = "SpatialInfluence")

mean_influence <- long_all %>%
  filter(!is.na(CellType)) %>%        # drop NA CellType rows
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    mean_spatial = mean(SpatialInfluence, na.rm = TRUE),
    n = dplyr::n(),
    n_nonzero = sum(SpatialInfluence > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_spatial))

# Drop tumor cells and re-order
mean_influence_no_tumor <- mean_influence %>%
  dplyr::filter(!is.na(CellType), !CellType %in% tumorcells) %>%
  dplyr::arrange(desc(mean_spatial))
mean_influence_no_tumor$CellType <- factor(mean_influence_no_tumor$CellType,
                                           levels = mean_influence_no_tumor$CellType)

print(mean_influence_no_tumor, n = Inf)

mean_influence <- mean_influence_no_tumor %>%
  arrange(desc(mean_spatial))
mean_influence$CellType <- factor(mean_influence$CellType, levels = mean_influence$CellType)

#Barplot for avergae spatial influence
ggp <- ggplot(mean_influence, aes(x = CellType, y = mean_spatial, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits = c(0, .06),
    expand = c(0, 0)
  ) +
  scale_fill_hue() +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = .2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Average Spatial Influence on Epithelial Cells by Cell Type",
    x = "Cell Type",
    y = "Average Spatial Influence",
    fill = "Cell Types"
  )

pdf("spatial_influence_barplot.pdf", width = 7, height = 5)
print(ggp)
dev.off()
