#Script to analyze LiP-MS data
#use Markdown provided by Malinovska et al., 2023 (Nature Protocols)

rm(list=ls())

#after Markdown "Analysis of a LiP-MS data set"

knitr::opts_chunk$set(include = FALSE)
# library(devtools)
# devtools::install_github("Vitek-Lab/MSstatsConvert")
# devtools::install_github("Vitek-Lab/MSstats")
# devtools::install_github("Vitek-Lab/MSstatsTMT")
# devtools::install_github("Vitek-Lab/MSstatsPTM")
# devtools::install_github("Vitek-Lab/MSstatsLiP")


#first install R packages specified in paper
#needs to be done only once 

# install.packages("devtools")
# install.packages("checkmate")
# install.packages("factoextra")
# install.packages("gghighlight")
# install.packages("gridExtra")
# install.packages("goeveg")
# install.packages("magrittr")
# install.packages('dendextend')
# install.packages('metafolio')
# install.packages('corrplot')
# install.packages('FactoMineR')
# install.packages("GGally")
# install.packages("usethis")
# library(devtools)
# BiocManager::install("MSstatsPTM")
# 
# install.packages("rmarkdown")
# devtools::install_github("Vitek-Lab/MSstatsLiP", build_vignettes = TRUE)
# install_github("vqv/ggbiplot")

library(MSstatsLiP)
library(gghighlight)
library(grid)
library(gridExtra)
library(magrittr)
library(dendextend)
library(GGally)
library(RColorBrewer)
library(ggbiplot)
library(dplyr)
library(data.table)
library(ggrepel)
library(tidyverse)
library(MSstatsPTM)
library(MSstats)
library(utils)

#set working directory
# setwd(r"[Z:\Users\Jasmin\_PhD\Research\Collaborations\RhoA (Dani, Scheffner group)\XL_experiment5_092023\202310_LiP\Spectronaut\R_script]")

#load file from Spectronaut (directDIA search; without trypsin; default MSStats Report)
#choose.files is an interactive function that opens a window dialog in which you can choose the corresponding .csv file

my_raw_lip <- read_delim(
  # file=choose.files(caption="Choose LiP dataset"), #choose "20240215_153149_F2310_JJ_RhoA_LiP_unspecific_newCont_Norm_noImp_Report_excerpt.csv"
  file="20240215_153149_F2310_JJ_RhoA_LiP_unspecific_newCont_Norm_noImp_Report_excerpt.csv",  #DB: direct path
  delim=",", escape_double = FALSE, trim_ws = TRUE) 

my_raw_prot <- read_delim(
  # file=choose.files(caption="Choose TrP dataset"), #select "20240215_181727_F2310_JJ_RhoA_Trp_specific_newCont_Norm_noImp_Report_excerpt.csv"
  file="20240215_181727_F2310_JJ_RhoA_Trp_specific_newCont_Norm_noImp_Report_excerpt.csv", #DB: direct path
  delim=",", escape_double = FALSE, trim_ws = TRUE) 

#load example data from MSStatsLiP Protocol paper

raw_lip_exmp <- read_delim(
  # file=choose.files(caption="Choose LiP dataset"), #select "2019-04-29_Correction_Protein_LiP_MSStats_Report_excerpt.csv"
  file="2019-04-29_Correction_Protein_LiP_MSStats_Report_excerpt.csv", # DB: Direct path
  delim=",", escape_double = FALSE, trim_ws = TRUE) 

raw_prot_exmp <- read_delim(
  # file=choose.files(caption="Choose TrP dataset"), #select "2019-04-29_Correction_Protein_Trp_MSStats_Report_excerpt.csv"
  file="2019-04-29_Correction_Protein_Trp_MSStats_Report_excerpt.csv", # DB: Direct path
  delim=",", escape_double = FALSE, trim_ws = TRUE) 

# change all protein names to Protein Accession Numbers

#raw_lip2[raw_lip2 == 'His-E6AP'] <- 'Q05086-2'
#raw_lip2[raw_lip2 == 'GST-16E6'] <- 'P03126'
#raw_lip2[raw_lip2 == 'RhoA_G14V'] <- 'P61586'

#raw_prot2[raw_prot2 == 'His-E6AP'] <- 'Q05086-2'
#raw_prot2[raw_prot2 == 'GST-16E6'] <- 'P03126'
#raw_prot2[raw_prot2 == 'RhoA_G14V'] <- 'P61586'

#remove irt rows from both dataframes

#my_raw_lip <- my_raw_lip[!grepl("iRT_Peptides_Fusion", my_raw_lip$PG.ProteinAccessions),]
#my_raw_prot <- my_raw_prot[!grepl("iRT_Peptides_Fusion", my_raw_prot$PG.ProteinAccessions),]

#raw_lip_exmp <- raw_lip_exmp[!grepl("iRT-Kit_WR_fusion", raw_lip_exmp$PG.ProteinAccessions),]
#raw_prot_exmp <- raw_prot_exmp[!grepl("iRT-Kit_WR_fusion", raw_prot_exmp$PG.ProteinAccessions),]

#load fasta file that was used
#important that Uniprot identifiers are used!!!
#change in dataframe and fasta 

# my_fasta_file=choose.files(caption = "Choose FASTA file") #select "RhoA_G14V_E6AP_E6_Exp2301_UniprotDB_inclFrankenfield_Contaminants.fasta" 
my_fasta_file <- "RhoA_G14V_E6AP_E6_Exp2301_UniprotDB_inclFrankenfield_Contaminants.fasta" # DB: Direct path

# fasta_file_exmp=choose.files(caption = "Choose FASTA file") #select "yeast_UP000002311_2024_01_07_inclSNCA.fasta"
fasta_file_exmp <- "yeast_UP000002311_2024_01_07_inclSNCA.fasta"


#convert files to MSStatsLib format (normally includes the trp dataset as well)
#by default, iRT peptides are removed removeiRT = TRUE. Also, peptides containing modifications are filtered by default, removeModifications.
#peptides with mutliple protein annotations are filtered as well, removeNonUniqueProteins = FALSE --> this should not be an issue for my dataset

#debug(SpectronauttoMSstatsLiPFormat)
my_msstats_data <- SpectronauttoMSstatsLiPFormat(my_raw_lip, my_fasta_file, Trp.data = my_raw_prot)

msstats_data_exmp <- SpectronauttoMSstatsLiPFormat(raw_lip_exmp, fasta_file_exmp, Trp.data = raw_prot_exmp) 

#my_msstats_data only contains accessions of my three main proteins
#lip_uni <- unique(my_raw_lip$PG.ProteinAccessions)
#trp_uni <- unique(my_raw_prot$PG.ProteinAccessions)
#intersect(lip_uni, trp_uni)

#Summarise data
#log2 transformation and normalization by median equalization


my_MSstatsLiP_Summarized <- dataSummarizationLiP(my_msstats_data, normalization.LiP = "equalizeMedians")

MSstatsLiP_Summarized_exmp <- dataSummarizationLiP(msstats_data_exmp, normalization.LiP = "equalizeMedians")

##inspect new dataframe
#names(MSstatsLiP_Summarized[["LiP"]])

#head(MSstatsLiP_Summarized[["LiP"]]$FeatureLevelData)

#head(MSstatsLiP_Summarized[["LiP"]]$ProteinLevelData)

#unique(MSstatsLiP_Summarized$LiP$ProteinLevelData$GROUP)

#unique(msstats_data[["LiP"]]$Condition)%in%unique(msstats_data[["TrP"]]$Condition)

##Modelling
#output is altered peptide levels

#debug(groupComparisonLiP)
#debug(groupComparisonPlotsPTM)
my_MSstatsLiP_model = groupComparisonLiP(my_MSstatsLiP_Summarized, fasta=my_fasta_file)
##output error for function groupComparison: "Error in `[.data.frame`(as.data.frame(comparisons), , cols) : undefined columns selected"

MSstatsLiP_model_exmp = groupComparisonLiP(MSstatsLiP_Summarized_exmp, fasta=fasta_file_exmp)
##no output error




#Intensity distribution
dataProcessPlotsLiP(MSstatsLiP_Summarized,
                    type = 'QCPLOT',
                    address = FALSE,
                    which.Peptide = "allonly",
                    lip.title = "All LiP Peptides",
                    protein.title = "All Proteins")

#Coefficient of variation for samples
MSstatsLiP_Summarized$LiP$FeatureLevelData%>%
  group_by(FEATURE, GROUP) %>% 
  summarize(cv = sd(INTENSITY) / mean(INTENSITY)) %>% 
  ggplot(aes(x = GROUP, y = cv, fill = GROUP)) + 
  geom_violin() + 
  theme_bw()+
  labs(title = "Coefficient of Variation between Condtions",
       subtitle = "LiP samples",
       y = "Coefficient of Variation", 
       x = "Conditon")+    
  scale_fill_brewer(palette = "Paired")+
  stat_summary(fun.data=function(x){return(data.frame(y=0, label=round(median(x, na.rm = TRUE),4)))}, geom="text", hjust=0.5, vjust=1)


#Trypticity of samples
trypticHistogramLiP(MSstatsLiP_Summarized, fasta_file,
                    color_scale = "grey",
                    address = FALSE)


#Number of peptide and protein identifications
grid.arrange(
  MSstatsLiP_Summarized$LiP$FeatureLevelData%>%
    filter(!is.na(ABUNDANCE))%>%
    group_by(GROUP,SUBJECT)%>%
    summarize(distinct_peptides=n_distinct(PEPTIDE))%>%
    ggplot(aes(x=SUBJECT, y=distinct_peptides))+
    geom_col(fill="white",color="black")+
    facet_grid(~GROUP, scales = "free")+
    theme_bw()+
    geom_text(aes(label=distinct_peptides), vjust=0, hjust=0.5)+
    labs(x="Replicate", 
         y="Number of identified peptides", 
         title="Number of identified peptides",
         subtitle="LiP sample"),
  
  MSstatsLiP_Summarized$LiP$FeatureLevelData%>%
    filter(!is.na(ABUNDANCE))%>%
    group_by(GROUP,SUBJECT)%>%
    summarize(distinct_proteins=n_distinct(PROTEIN))%>%
    ggplot(aes(x=SUBJECT, y=distinct_proteins))+
    geom_col(fill="white",color="black")+
    facet_grid(~GROUP, scales = "free")+
    theme_bw()+
    geom_text(aes(label=distinct_proteins), vjust=0, hjust=0.5)+
    labs(x="Replicate", 
         y="Number of identified proteins", 
         title="Number of identified proteins",
         subtitle="LiP sample"),
  
  ncol=2)


#Correlation of data
correlationPlotLiP(MSstatsLiP_Summarized, address = FALSE)

color_fn <- function(data, mapping, method="p", use="complete.obs", min_val=0.5, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method=method, use=use)
  colFn <- colorRampPalette(c("white", "steelblue"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(min_val, 1, length=100))]
  ggally_cor(data = data, mapping = mapping, color="black",method=method, use=use,...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

ggpairs(MSstatsLiP_Summarized$LiP$FeatureLevelData%>%
          ungroup()%>%
          select(INTENSITY, SUBJECT, FEATURE)%>%
          spread(., SUBJECT, INTENSITY)%>%
          select(-FEATURE), 
        upper = list(continuous = wrap(color_fn, min_val=0.97)),
        lower = list(continuous = wrap("points", alpha=0.3)),
        diag = list(continuous = "densityDiag"))+
  ggtitle("LiP samples")


#Protein abundance

adj.pvalue.cutoff=0.01
log2FC.cutoff=2

#non-adjusted
grid.arrange(
  MSstatsLiP_model$LiP.Model%>%
    filter(is.na(issue)&abs(log2FC)>log2FC.cutoff&adj.pvalue<adj.pvalue.cutoff)%>%
    group_by(Label)%>%
    summarize(distinct_peptides=n_distinct(PeptideSequence))%>%
    ggplot(aes(x=Label, y=distinct_peptides))+
    geom_col(fill="white",color="black")+
    theme_bw()+
    geom_text(aes(label=distinct_peptides), vjust=0, hjust=0.5)+
    labs(x="Comparison", 
         y="Number of peptides", 
         title="Structurally altered peptides",
         subtitle="LiP sample, non-adjusted"),
  
  MSstatsLiP_model$LiP.Model%>%
    filter(is.na(issue)&abs(log2FC)>log2FC.cutoff&adj.pvalue<adj.pvalue.cutoff)%>%
    group_by(Label)%>%
    summarize(distinct_proteins=n_distinct(ProteinName))%>%
    ggplot(aes(x=Label, y=distinct_proteins))+
    geom_col(fill="white",color="black")+
    theme_bw()+
    geom_text(aes(label=distinct_proteins), vjust=0, hjust=0.5)+
    labs(x="Comparison", 
         y="Number of proteins", 
         title="Structurally altered proteins",
         subtitle="LiP sample, non-adjusted"),
  nrow=1)

#adjusted
grid.arrange(
  MSstatsLiP_model$Adjusted.LiP.Model%>%
    filter(is.na(issue)&abs(log2FC)>log2FC.cutoff&adj.pvalue<adj.pvalue.cutoff)%>%
    group_by(Label)%>%
    summarize(distinct_peptides=n_distinct(PeptideSequence))%>%
    ggplot(aes(x=Label, y=distinct_peptides))+
    geom_col(fill="white",color="black")+
    theme_bw()+
    geom_text(aes(label=distinct_peptides), vjust=0, hjust=0.5)+
    labs(x="Comparison", 
         y="Number of peptides", 
         title="Structurally altered peptides",
         subtitle="LiP sample, adjusted"),
  
  MSstatsLiP_model$Adjusted.LiP.Model%>%
    filter(is.na(issue)&abs(log2FC)>log2FC.cutoff&adj.pvalue<adj.pvalue.cutoff)%>%
    group_by(Label)%>%
    summarize(distinct_proteins=n_distinct(ProteinName))%>%
    ggplot(aes(x=Label, y=distinct_proteins))+
    geom_col(fill="white",color="black")+
    theme_bw()+
    geom_text(aes(label=distinct_proteins), vjust=0, hjust=0.5)+
    labs(x="Comparison", 
         y="Number of proteins", 
         title="Structurally altered proteins",
         subtitle="LiP sample, adjusted"),
  nrow=1)

#Try Barcode plots
#Careful function from script doesn't exist anymore
#following command creates plots for all possible comparisons
#default settings are FC.cutoff 0, adj.pvalue.cutoff is 0.05
StructuralBarcodePlotLiP(MSstatsLiP_model, fasta_file,
                         model_type = "Unadjusted",
                         address=FALSE)

#add parameters for individualized plot
StructuralBarcodePlotLiP(MSstatsLiP_model, fasta_file,
                         model_type = "Unadjusted",
                         FC.cutoff = 2, adj.pvalue.cutoff = 0.01,
                         which.comp = "Wt vs Wt + OF232",
                         height = 3,
                         width = 10,
                         address=FALSE)

#plot is not very informative since peptides with log2FC Inf values can't be used for filtering properly
#try to filter out all rows that contain (-)Inf values and check if result is better

MSstatsLiP_model$LiP.Model <- filter(MSstatsLiP_model$LiP.Model, log2FC != 'Inf')
MSstatsLiP_model$LiP.Model <- filter(MSstatsLiP_model$LiP.Model, log2FC != '-Inf')

#use only negative FC ratios 

MSstatsLiP_model$LiP.Model <- filter(MSstatsLiP_model$LiP.Model, log2FC > 0)

#test Barcode again

StructuralBarcodePlotLiP(MSstatsLiP_model, fasta_file,
                         model_type = "Unadjusted",
                         FC.cutoff = 2, adj.pvalue.cutoff = 0.01,
                         which.comp = "wt vs wt + E6",
                         height = 3,
                         width = 10,
                         address=FALSE)

#better results than before but dataset still includes peptides yielded from this algorithm with low p-values and high FC
#that Spectronaut yields as far less significant and also look inconsistently quantified when checking profile plots

#try to feed own dataframe into algorithm for figure
#necessary columns are

#export MSstatsLiP_model, fasta file for github issue

write.csv(MSstatsLiP_model$LiP.Model, "E6APwt_unspecific_imputed_RIGHTfile_afterMSstatsLiP.csv", row.names = FALSE)
                         