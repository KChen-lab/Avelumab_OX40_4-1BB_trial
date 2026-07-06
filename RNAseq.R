## This script includes RNAseq data analysis of manuscript titled as "Impact of avelumab and immune-stimulating agent combinations on the tumor microenvironment in patients with advanced solid tumors"

# RNAseq data 

set.seed(2024)

library(ggplot2)
library(data.table)
library(DESeq2)
library(dplyr)
library(tximport)

#read in rsem

rsem.data = tximport(f,type = "rsem", txIn = F, txOut = F)

## Remove genes that has length ==0
subset_txi = function(tx,x){
  #this is specifically written for rsem import
  var = lapply(tx[1:3], function(i){
    i[,x]
  })
  var[["countsFromAbundance"]] = "no"
  co = var$length 
  co[co == 0] =1
  var$length = co
  return(var)
}

coldata = read.csv("Pfizer_RNA_metadata_paired_major_disease_ARM.csv", row.names =1)

rsem.data = subset_txi(rsem.data,rownames(coldata))

## Analysis with DESEQ2 using tximport

coldata$dis_arm_time = paste(coldata$Disease, coldata$Arm, coldata$TimePoint, sep = "_")

dds = DESeqDataSetFromTximport(txi = rsem.data,colData = coldata,design = ~  Batch + dis_arm_time) 

smallestGroupSize <- 50
keep <- rowSums(counts(dds) >= 100) >= smallestGroupSize
dds <- dds[keep,]

dds = estimateSizeFactors(dds)
dds = DESeq(dds)
resultsNames(dds)



## Fig 3 Pathway analysis

library(fgsea)

## Fig 3A Cholangiocarcinoma

# Fig 3A.1
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Bile Duct|Cholangocarcinoma_A_C1D15","Bile Duct|Cholangocarcinoma_A_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3A_1.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Bile Duct ARM A") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

## Fig 3A_2 Cholangiocarcinoma

resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Bile Duct|Cholangocarcinoma_C_C1D15","Bile Duct|Cholangocarcinoma_C_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3A_2.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Bile Duct ARM C") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

## Fig 3A_3 Cholangiocarcinoma

resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Bile Duct|Cholangocarcinoma_C_C3D15","Bile Duct|Cholangocarcinoma_C_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3A_3.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C3D15 and baseline of Bile Duct ARM C") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

###########

## Fig 3B Pancreatic

# Fig 3B.1
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Pancreatic_B_C1D15","Pancreatic_B_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3B_1.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Pancreatic ARM B") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

# Fig 3B.2
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Pancreatic_B_C3D15","Pancreatic_B_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3B_2.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C3D15 and baseline of Pancreatic ARM B") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()


# Fig 3B.3
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Pancreatic_C_C1D15","Pancreatic_C_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3B_3.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Pancreatic ARM C") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

# Fig 3B.4
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Pancreatic_C_C3D15","Pancreatic_C_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3B_4.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C3D15 and baseline of Pancreatic ARM C") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

#################################

# Fig 3C
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Ovarian_B_C1D15","Ovarian_B_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3C.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Ovarian ARM B") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

#########
# Fig 3D_1 Endometrial 
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Endometrial_B_C1D15","Endometrial_B_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3D_1.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Endometrial ARM B") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

# Fig 3D_2 Endometrial 
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Endometrial_B_C3D15","Endometrial_B_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3D_2.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C3D15 and baseline of Endometrial ARM B") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()


# Fig 3E Colorectal
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Colorectal_C_C1D15","Colorectal_C_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3E.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Colorectal ARM C") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

##########

# Fig 3F_1 Cervical
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Cervical_A_C1D15","Cervical_A_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3F_1.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Cervical ARM A") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()


# Fig 3F_2 Cervical
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Cervical_A_C3D15","Cervical_A_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3F_2.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C3D15 and baseline of Cervical ARM A") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

# Fig 3F_4 Cervical
resultsNames(dds)
res <- results(dds,
               contrast = c("dis_arm_time","Cervical_D_C1D15","Cervical_D_Baseline"))
summary(res)
res<- data.frame(res)
res1<- res[complete.cases(res), ]

## Replace ensembl id with symbols
genes <- read.delim("human_gene_name.txt")
res1=merge(genes,as.data.frame(res1),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
res1<- res1 %>% distinct(gene_name, .keep_all = TRUE)

rownames(res1) <- res1[,2]
res1[,2] <- NULL
#remove gene_ID column 
res1[1] <- NULL
res1 <- data.frame(res1)
res1<- res1[complete.cases(res1), ]
res1$symbol <- rownames(res1)

res2 <- res1 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
head(exampleRanks)

barplot(sort(ranks, decreasing = T))

pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark%>% 
  head() %>% 
  lapply(head)
fgseaResP<- fgsea(pathways=pathways.hallmark, stats=ranks, minSize =15, maxSize = 500)

fgseaResTidy <- fgseaResP %>%
  as_tibble() %>%
  arrange(desc(NES))


#Significant Pathways
fgseaRes_Tidy_1 <- subset(fgseaResTidy, padj<0.05)

fgseaRes_Tidy_1 <- apply(fgseaRes_Tidy_1,2,as.character)

pdf("Fig_3F_4.pdf")
Bubble <- ggplot(fgseaRes_Tidy_1, aes(reorder(pathway, NES), NES, size=size, color=padj))+
  geom_point(alpha=0.8)+coord_flip() +
  labs(x="Significant Hallmark Pathways", y="Normalized Enrichment Score",
       title="Significant Hallmark pathways NES from GSEA between C1D15 and baseline of Cervical ARM D") + 
  theme_minimal()+ scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(4, 8), name="Coverage, [%]")
Bubble
dev.off()

#####################################################

## Fig 4

#Normalized data
norm <- assay(vsd)
##### Batch correction
norm_batch<- limma::removeBatchEffect(norm, batch = coldata$Batch)

#gene names
norm_batch = read.csv("Pfizer_Trial_RNAseq_normalized_batch_corrected_major_disease_data_keep100_by50.csv", row.names=1)
genes <- read.delim("human_gene_name.txt")
namevsd=merge(genes,as.data.frame(norm_batch),by.x='ensembl_ID',by.y=0)

# Remove duplicated rows based on Gene_ID
namevsd<- namevsd %>% distinct(gene_name, .keep_all = TRUE)

rownames(namevsd) <- namevsd[,2]
namevsd[,2] <- NULL
#remove gene_ID column 
namevsd[1] <- NULL

## CibersortX with absolute mode was run in Cibersort website using batch corrected normalized data.

library(tidyverse)
library(ggpubr)
#Immunedeconvolution
cibersort <- read.csv("Cibersort.csv")
colnames(cibersort)
cibersort_samples <- select(cibersort,-`P.value_CIBERSORT`,-`Correlation_CIBERSORT`,
                            -`RMSE_CIBERSORT`,-`Absolute_score_.sig_score._CIBERSORT`)


long_cibersort_samples <- gather(cibersort_samples, cells, values, B_cells_naive_CIBERSORT:Neutrophils_CIBERSORT, factor_key=T)

## CD4 T cells (Combination of 
#################1. T_cells_CD4_naive_CIBERSORT,
#################2. T_cells_CD4_memory_resting_CIBERSORT,
#################3. T_cells_CD4_memory_activated_CIBERSORT,
#################4. T_cells_regulatory_(Tregs)_CIBERSORT)

T_CD4 <- long_cibersort_samples[long_cibersort_samples$cells== "T_cells_CD4",]

p1 <-ggplot(T_CD4,aes(x= TimePoint,y = values, color=Arm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Arm),size = 1,position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.7))+
  #geom_line(aes(group = MRN), colour = "black")+
  theme_bw()+ ggtitle("Cibersortx Scores - CD4 T cells")+
  ylab("Scores")+ xlab("")+
  facet_wrap(~Disease, scales = "free", ncol=3)+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=18),
        title = element_text(size=18),
        strip.text= element_text(size=18))

pdf(file="Fig_4_A.pdf", height = 8, width = 11)

p1
dev.off()

## CD8 T cells
T_CD8 <- long_cibersort_samples[long_cibersort_samples$cells== "T_cells_CD8_CIBERSORT",]

p1 <-ggplot(T_CD8,aes(x= TimePoint,y = values, color=Arm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Arm),size = 1,position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.7))+
  #geom_line(aes(group = MRN), colour = "black")+
  theme_bw()+ ggtitle("Cibersortx Scores - CD4 T cells")+
  ylab("Scores")+ xlab("")+
  facet_wrap(~Disease, scales = "free", ncol=3)+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=18),
        title = element_text(size=18),
        strip.text= element_text(size=18))

pdf(file="Fig_4_B.pdf", height = 8, width = 11)

p1
dev.off()

## Fig 4C - Macrophages_M2
M2 <- long_cibersort_samples[long_cibersort_samples$cells== "Macrophages_M2_CIBERSORT",]

p1 <-ggplot(M2,aes(x= TimePoint,y = values, color=Arm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Arm),size = 1,position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.7))+
  #geom_line(aes(group = MRN), colour = "black")+
  theme_bw()+
  ylab("Scores")+ xlab("")+
  facet_wrap(~Disease, scales = "free", ncol=3)+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=18),
        title = element_text(size=18),
        strip.text= element_text(size=18))

pdf(file="Fig_4_C.pdf", height = 8, width = 11)

p1
dev.off()

