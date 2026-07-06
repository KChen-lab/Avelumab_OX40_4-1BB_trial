
## This script includes WES data analysis of manuscript titled as "Impact of avelumab and immune-stimulating agent combinations on the tumor microenvironment in patients with advanced solid tumors"

# WES data 

## Create maf file from GATK-mutect pipeline to visualize somatic mutations.

set.seed(2024)
library(maftools)

mut_pin <- read.delim("Pfizer_Trial_mutect_pindel_specific_columns_major_disease_ARM.maf")

WEX_meta = read.csv(file = "Pfizer_WEX_metadata_major_disease_ARM.csv")
samples <- WEX_meta[c(WEX_meta$TimePoint != "PBMC Sample"),]
maf = read.maf(maf = mut_pin,removeDuplicatedVariants =T,useAll = T, clinicalData = samples)

maf@variants.per.sample$Tumor_Sample_Barcode

## Calculate TMB

tmb <- tmb(maf, captureSize = 51, logScale = F)
tmb_d <- select(tmb,Tumor_Sample_Barcode,total_perMB )
tmb_d <- as.data.frame(tmb_d)

## Fig S1

library(ggplot2)
library(tidyverse)
library(ggpubr)
pdf(file="Fig_S1.pdf", h=5, w=8)
ggplot(tmb_d, aes(x =Disease, y = total_perMB)) +
  theme_bw() + 
  geom_boxplot(aes(fill=TimePoint)) +
  geom_point(position=position_dodge(width=0.75),aes(group=TimePoint))+
  ggtitle("Pfizer Trial Tumor Mutation Burden (TMB)") +
  ylab("TMB / MB") + xlab("")+
  scale_color_manual(values=c("#f8766c","#7bad00","#00bec4","#c67bff"))
dev.off()
  
## Fig 2.A

pdf(file="Fig_2A.pdf")
oncoplot(maf = maf,showTumorSampleBarcodes =FALSE,sortByMutation = T,drawColBar = T,topBarData = tmb_d,
         sampleOrder = c("222453",'222454','222455','222463','222464','226425','226426','226446','226447','226435',
                         '226436','226422','226423','226424','226432','226433','226434','222477','222478','222479',
                         '222482','222483','222484','222486','222487','222491','222492',"222427",'222428','222429',
                         '222456','222457','222430','222431','222432','222433','222434','222461','222462','222424',
                         '222465','222480','222481','222469','222470','222474','222475','222471','222472','222443',
                         '222444','222445','222435','222436','222437','222440','222441','222442','222438','222439',
                         '222448','222449','222450','222467','222468','222488','222489',
                         '226459','226460','226463','226464','226453','226454','226455','226461','226462','226448',
                         '226449','226450','226440','226441','226442','226451','226452','226430','226431','226427',
                         '226428','226429','222412','222413','222409','222410','222411','222417','222418','222419',
                         '222414','222415','222416','222420','222421','226457','226458','222497','222498','222495',
                         '222496','222493','222494'),
         clinicalFeatures = c("Disease","TimePoint","Arm","Patient","Clinical.benefits"),bgCol="gray",gene_mar = 7,sortByAnnotation = TRUE,
         removeNonMutated = F, top=40)
dev.off()

## Fig 2.B
pws = pathways(maf = maf, plotType = 'treemap')
pdf(file="Fig_2B.pdf")

plotPathways(maf = maf, pathlist = pws,sampleOrder = c("222453",'222454','222455','222463','222464','226425','226426','226446','226447','226435',
                                                       '226436','226422','226423','226424','226432','226433','226434','222477','222478','222479',
                                                       '222482','222483','222484','222486','222487','222491','222492',"222427",'222428','222429',
                                                       '222456','222457','222430','222431','222432','222433','222434','222461','222462','222424',
                                                       '222465','222480','222481','222469','222470','222474','222475','222471','222472','222443',
                                                       '222444','222445','222435','222436','222437','222440','222441','222442','222438','222439',
                                                       '222448','222449','222450','222467','222468','222488','222489',
                                                       '226459','226460','226463','226464','226453','226454','226455','226461','226462','226448',
                                                       '226449','226450','226440','226441','226442','226451','226452','226430','226431','226427',
                                                       '226428','226429','222412','222413','222409','222410','222411','222417','222418','222419',
                                                       '222414','222415','222416','222420','222421','226457','226458','222497','222498','222495',
                                                       '222496','222493','222494'),
             showTumorSampleBarcodes =FALSE)

dev.off()

## Fig 2.C
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('NMF')
maf.tnm = trinucleotideMatrix(maf = maf, prefix = 'chr', add = TRUE,useSyn = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

maf.sign = estimateSignatures(mat = maf.tnm, nTry = 6)

maf.sig = extractSignatures(mat = maf.tnm, n = 3)

pdf(file="Fig_2C.pdf")
maftools::plotSignatures(nmfRes = maf.sig, title_size = 1.2, sig_db = "SBS",contributions = TRUE,show_barcodes = FALSE,
                         patient_order =  c("222453",'222454','222455','222463','222464','226425','226426','226446','226447','226435',
                                                       '226436','226422','226423','226424','226432','226433','226434','222477','222478','222479',
                                                       '222482','222483','222484','222486','222487','222491','222492',"222427",'222428','222429',
                                                       '222456','222457','222430','222431','222432','222433','222434','222461','222462','222424',
                                                       '222465','222480','222481','222469','222470','222474','222475','222471','222472','222443',
                                                       '222444','222445','222435','222436','222437','222440','222441','222442','222438','222439',
                                                       '222448','222449','222450','222467','222468','222488','222489',
                                                       '226459','226460','226463','226464','226453','226454','226455','226461','226462','226448',
                                                       '226449','226450','226440','226441','226442','226451','226452','226430','226431','226427',
                                                       '226428','226429','222412','222413','222409','222410','222411','222417','222418','222419',
                                                       '222414','222415','222416','222420','222421','226457','226458','222497','222498','222495',
                                                       '222496','222493','222494'))

dev.off()

## Fig 2.D and 2.E

pdf(file="Fig_2D_E.pdf")
plotApobecDiff(tnm = maf.tnm, maf = maf)
dev.off()

## Fig 2.F

CNV <- read.csv("Pfizer_trial_WES_exomecn_gene_annotated.csv")

## Select samples only have baseline and pbmc and a timepoint
samples <- WEX_meta[c(WEX_meta$TimePoint != "PBMC Sample"),]
samples$SpecimenID <- paste("X", samples$SpecimenID, sep="")
select <-samples$SpecimenID 

# Remove duplicated
# Remove duplicates based on Sepal.Width columns
CNV <- CNV[!duplicated(CNV$genename), ]
# Remove unwanted columns for heatmap matrix

rownames(CNV) <- CNV$genename

CNV[1:5] <- NULL

CNV<- CNV[,select]


library(ComplexHeatmap)
library(circlize)


CNV_rep_gene <- CNV[c("PTEN","TP53","CDKN2A","PIK3CA","MDM2","CDK4","BRAF","FGFR1","ALK","FGFR3",
           "MET","EGFR","FLT3","HRAS","KRAS","TSC1","BAP1","BRCA1","BRCA2",'CTNNB1',
           'RAD51C',"CCNE1","PARP","EMSY"),]

CNV_rep_gene <- na.omit(CNV_rep_gene)
col_fun = c("D"="blue","N"="#e5e1e1","A"="red")

col_fun = colorRamp2(c(-0.8,0,0.8), c("blue", "white", "red"))

ann <- data.frame(
  TimePoint =samples$TimePoint,
  Arm =samples$Arm,
  Disease = samples$Disease,
  stringsAsFactors = FALSE)

# create the colour mapping
colours <- list(Disease =c("Pancreatic"="brown", "Cervical"= "yellow4","Ovarian"= "green","Colorectal"= "blue",
                          "Endometrial"= "magenta","Bile Duct|Cholangocarcinoma"="purple"),
                Arm = c("A"= 'yellow',"B"="blue","C"="#ea1b82","D"="black"),
  TimePoint = c('Baseline'='#f8766d', 'TP-2'='#7cae00','TP-3'="#00bfc4","TP-4"="#0e1d9e"))


# now create the ComplexHeatmap annotation object
# as most of these parameters are self-explanatory, comments will only appear where needed
###############

###############################################3

colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'black', # default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    TimePoint= list(
      nrow =5, # number of rows across which the legend will be arranged
      title = 'TimePoint',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12)),
    Arm = list(
      nrow = 4,
      title = 'Arm',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12)),
    Disease= list(
      nrow = 6, # number of rows across which the legend will be arranged
      title = 'Disease',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12))))


ht_list = Heatmap(CNV_rep_gene,name="log2", col=col_fun,cluster_columns = F,cluster_rows = T,
                  column_order=c("X222453",'X222454','X222455','X222463','X222464','X226425','X226426','X226446','X226447','X226435',
                                 'X226436','X226422','X226423','X226424','X226432','X226433','X226434','X222477','X222478','X222479',
                                 'X222482','X222483','X222484','X222486','X222487','X222491','X222492',"X222427",'X222428','X222429',
                                 'X222456','X222457','X222430','X222431','X222432','X222433','X222434','X222461','X222462','X222424',
                                 'X222465','X222480','X222481','X222469','X222470','X222474','X222475','X222471','X222472','X222443',
                                 'X222444','X222445','X222435','X222436','X222437','X222440','X222441','X222442','X222438','X222439',
                                 'X222448','X222449','X222450','X222467','X222468','X222488','X222489',
                                 'X226459','X226460','X226463','X226464','X226453','X226454','X226455','X226461','X226462','X226448',
                                 'X226449','X226450','X226440','X226441','X226442','X226451','X226452','X226430','X226431','X226427',
                                 'X226428','X226429','X222412','X222413','X222409','X222410','X222411','X222417','X222418','X222419',
                                 'X222414','X222415','X222416','X222420','X222421','X226457','X226458','X222497','X222498','X222495',
                                 'X222496','X222493','X222494'),
                  heatmap_legend_param = list(
                    color_bar = 'continuous',
                    title_position = 'topleft',
                    title_gp=gpar(fontsize = 12, fontface = 'bold'),
                    labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                  bottom_annotation = colAnn, 
                  #clustering_distance_rows = "spearman",clustering_distance_columns = "spearman",
                  row_names_side = "left",#rect_gp = gpar(col = "black", lwd = 0.5),
                  show_row_names = T, show_column_names = F,
                  row_names_gp = gpar(fontsize=8),
                  column_title = "CNV")

heatmap <- draw(ht_list,
                heatmap_legend_side = 'right',
                annotation_legend_side = 'right',
                row_sub_title_side = 'left')


pdf(file = "Fig_2F.pdf",w=12, h=8)
print(heatmap)
dev.off()


