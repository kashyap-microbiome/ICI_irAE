#for Borenstein metabolomics
require(readr)

#for curatedMetagenomicData
library(dplyr)
library(DT)
require(curatedMetagenomicData)
library(stringr)
library(mia)
library(scater)
library(vegan)
library(TreeSummarizedExperiment)
library(lmerTest)
require(RColorBrewer)
library(ggplot2)
library(patchwork)


#----------------------------------
#setwd to where this file is
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file_dir) 

source("helper functions/MannWhitney_function.R")


#----------------------------------
#downloaded Borenstein meta-analysis Repo
#there are 3 studies on IBD with controls, Jacobs, iHMP, Franzosa
#orotic acid is annotated. Get these metabolite values and associate with IBD subtypes


#download the meta-analysis repo to the following folder: 
setwd("../data/microbiome-metabolome-curated-data-main/")

sel_colors <- brewer.pal(8, "Dark2")
plot(1:8, 1:8, pch=16, cex=15, col=sel_colors)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Jacobs

mx_jacobs <- as.data.frame(read_tsv("../data/processed_data/JACOBS_IBD_FAMILIES_2016/mtb.tsv"))
meta_jacobs <- as.data.frame(read_tsv('../data/processed_data/JACOBS_IBD_FAMILIES_2016/metadata.tsv'))

dim(mx_jacobs); dim(meta_jacobs)

rownames(mx_jacobs) <- mx_jacobs$Sample
rownames(meta_jacobs) <- meta_jacobs$Sample
mx_jacobs <- mx_jacobs[rownames(meta_jacobs),]
cbind(rownames(mx_jacobs), rownames(meta_jacobs))

#from corresponding mtb.map.tsv file
#Compound	Compound.Name	ESI.mode	m.z	Retention.time	Putative.ID.that.did.not.validate	HMDB.identification	KEGG.identification	LIPIDMAPS.identification	BioCyc.identification	HMDB	KEGG	High.Confidence.Annotation
#Negative_155.0084_0.9167	NA	Negative	155.0084	0.9167	NA	1/1,[Orotic acid]	1/1,[C00295: Orotate; Orotic acid; Uracil-6-carboxylic acid ]	0	4,[12-PROPANEDIOL-1-PHOSPHATE: 1,2-propanediol 1-phosphate],[CPD-629: uracil 5-carboxylate],[CPD0-1458: 1,3-propanediol-P],[OROTATE: orotate]	NA	NA	TRUE

oro_col <- grep("Negative_155.0084_0.9167", colnames(mx_jacobs))
oro_vec_jacobs <- mx_jacobs[,oro_col]

table(meta_jacobs$Study.Group)
#CD Normal     UC 
#26     54     10 


#only 10 UC patients; very small cohort
#two datasets with more than 100 samples; ihmp and Franzosa 

boxplot(oro_vec_jacobs ~ meta_jacobs$Study.Group, las=1)
summary(lm(oro_vec_jacobs ~ meta_jacobs$Study.Group))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#IHMP study

mx_ihmp <- as.data.frame(read_tsv("../data/processed_data/iHMP_IBDMDB_2019/mtb.tsv"))
meta_ihmp <- as.data.frame(read_tsv('../data/processed_data/iHMP_IBDMDB_2019/metadata.tsv'))


dim(mx_ihmp); dim(meta_ihmp)

rownames(mx_ihmp) <- mx_ihmp$Sample
rownames(meta_ihmp) <- meta_ihmp$Sample
mx_ihmp <- mx_ihmp[rownames(meta_ihmp),]
cbind(rownames(mx_ihmp), rownames(meta_ihmp))


#inspecting
colnames(mx_ihmp)[grep("oro", tolower(colnames(mx_ihmp)))]
oro_col <- grep("HILn_QI93__orotate", colnames(mx_ihmp))
oro_vec_ihmp <- mx_ihmp[,oro_col]


#"4-Hydroxybenzoic.acid" "Alpha-aminobutyric.acid" "Cytosine" "Indolelactic.acid" "Orotic.acid" "Xanthine"  

#ILA (not annotated), 
#colnames(mx_ihmp)[grep("indo", tolower(colnames(mx_ihmp)))] #ILA not annotated
#colnames(mx_ihmp)[grep("benz", tolower(colnames(mx_ihmp)))] #"4-Hydroxybenzoic.acid" not annotated
colnames(mx_ihmp)[grep("butyr", tolower(colnames(mx_ihmp)))] #"Alpha-aminobutyric.acid" annotated as "HILn_QI3__2-aminobutyrate" 
colnames(mx_ihmp)[grep("cyto", tolower(colnames(mx_ihmp)))] #"Cytosine" annotated as "HILp_QI473__cytosine"
colnames(mx_ihmp)[grep("xanth", tolower(colnames(mx_ihmp)))] #"Xanthine" annotated as "HILn_QI130__xanthine"

hmb_col <- grep("HILn_QI3__2-aminobutyrate", colnames(mx_ihmp))
ab_vec_ihmp <- mx_ihmp[,hmb_col]

hmb_col <- grep("HILp_QI473__cytosine", colnames(mx_ihmp))
cyt_vec_ihmp <- mx_ihmp[,hmb_col]

hmb_col <- grep("HILn_QI130__xanthine", colnames(mx_ihmp))
xan_vec_ihmp <- mx_ihmp[,hmb_col]



#----------------------------------
#longitudinal? Yes, correct for Subject column
head(meta_ihmp)
table(meta_ihmp$Subject)

table(meta_ihmp$Study.Group)
#CD    nonIBD   UC 
#177    104    101 

#how many subjects and number of samples per subjects
## 1) Number of unique subjects per group
meta_ihmp %>%
  group_by(Study.Group) %>%
  summarise(n_subjects = n_distinct(Subject))
# Study.Group        n_subjects
#   1 CD                  n=49 (d=177, median=4)
#   2 UC                  n=30 (d=101, median=3)
#   3 nonIBD              n=26 (d=104, median=4)

meta_ihmp %>%
  dplyr::count(Study.Group, Subject) %>%        # samples per subject within each group
  group_by(Study.Group) %>%
  summarise(median_samples_per_subject = median(n))
# Study.Group             median_samples_per_subject
#   1 CD                                   4
#   2 UC                                   3
#   3 nonIBD                               4
#----------------------------------


boxplot(log10(oro_vec_ihmp) ~ meta_ihmp$Study.Group, las=1)

#not corrected for subject
summary(lm(oro_vec_ihmp ~ meta_ihmp$Study.Group))


#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                   169267      16708  10.131  < 2e-16 ***
#  meta_ihmp$Study.GroupnonIBD    13684      27464   0.498  0.61859    
#  meta_ihmp$Study.GroupUC       -87513      27720  -3.157  0.00172 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 222300 on 379 degrees of freedom
#Multiple R-squared:  0.0335,	Adjusted R-squared:  0.0284 
#F-statistic: 6.568 on 2 and 379 DF,  p-value: 0.001571


#----------------------------------
#correct for subject

meta_ihmp$oro_vec_ihmp <- log10(oro_vec_ihmp)
meta_ihmp$Study.Group <- relevel(factor(meta_ihmp$Study.Group), ref = "nonIBD")

model <- lmer(oro_vec_ihmp ~ Study.Group + (1 | Subject), data = meta_ihmp)
summary(model)
#                   Estimate Std. Error       df t value Pr(>|t|)    
#   (Intercept)    5.07037    0.05813 86.15166  87.227  < 2e-16 ***
#   Study.GroupCD -0.11481    0.07244 88.65637  -1.585    0.117    
#   Study.GroupUC -0.37747    0.08108 92.20932  -4.655 1.08e-05 ***

p.vals_ihmp <- "<0.0001"


#----------------------------------

hmb_col <- grep("HILn_QI3__2-aminobutyrate", colnames(mx_ihmp))
ab_vec_ihmp <- mx_ihmp[,hmb_col]

hmb_col <- grep("HILp_QI473__cytosine", colnames(mx_ihmp))
cyt_vec_ihmp <- mx_ihmp[,hmb_col]

hmb_col <- grep("HILn_QI130__xanthine", colnames(mx_ihmp))
xan_vec_ihmp <- mx_ihmp[,hmb_col]



#aminobutyrate; higher in UC; not relevant to plot
meta_ihmp$ab_vec_ihmp <- log10(ab_vec_ihmp)
model <- lmer(ab_vec_ihmp ~ Study.Group + (1 | Subject), data = meta_ihmp)
summary(model) #sign in CD p 0.00651.
#               Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)    4.94160    0.07670 94.70546  64.429  < 2e-16 ***
#Study.GroupCD  0.26515    0.09532 96.41482   2.782  0.00651 ** 
#Study.GroupUC  0.16809    0.10627 99.33654   1.582  0.11688    
boxplot(ab_vec_ihmp ~ Study.Group, data = meta_ihmp)


#cytosine
meta_ihmp$cyt_vec_ihmp <- log10(cyt_vec_ihmp)
model <- lmer(cyt_vec_ihmp ~ Study.Group + (1 | Subject), data = meta_ihmp)
summary(model) #not sign
#not mentioned in manuscript
#               Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)    5.24724    0.09366 94.00457  56.023   <2e-16 ***
#Study.GroupCD -0.15606    0.11572 95.98840  -1.349    0.181    
#Study.GroupUC  0.01923    0.12872 98.92874   0.149    0.882    
boxplot(cyt_vec_ihmp ~ Study.Group, data = meta_ihmp)



#xanthine
meta_ihmp$xan_vec_ihmp <- log10(xan_vec_ihmp)
model <- lmer(xan_vec_ihmp ~ Study.Group + (1 | Subject), data = meta_ihmp)
summary(model) #sign in UC, p lower 0.0265
#Xanthine is mentioned in Extended Data Fig. 4 associated to dysbiosis; 
#"Top 10 differentially abundant metabolites during dysbiosis not shown elsewhere in this manuscript"

#               Estimate    Std. Error  df    t value Pr(>|t|)    
#(Intercept)    6.30958    0.06378 87.06292  98.921   <2e-16 ***
#Study.GroupCD -0.05670    0.07947 89.50593  -0.714   0.4774    
#Study.GroupUC -0.20049    0.08892 93.02647  -2.255   0.0265 *  
boxplot(xan_vec_ihmp ~ Study.Group, data = meta_ihmp)

ihmp_p_xanth_cd <- 0.4774
ihmp_p_xanth_uc <- 0.0265



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#as the sampling is longitudinal it would be ideal to associate with symptom severity; do this only for UC

full_meta_ihmp <- as.data.frame(read.csv('../data/processed_data/iHMP_IBDMDB_2019/hmp2_metadata_2018-08-20.csv'))

full_meta_ihmp_metabolomics <- full_meta_ihmp[full_meta_ihmp$data_type == "metabolomics" & full_meta_ihmp$diagnosis == "UC",]
dim(full_meta_ihmp_metabolomics)

#number of timepoints per subject
table(full_meta_ihmp_metabolomics$Participant.ID)
#C3003 C3004 C3005 C3006 C3011 C3013 C3015 C3029 C3032 C3034 C3037 E5004 H4010 H4019 H4027 H4035 H4040 H4042 H4044 M2026 M2064 M2069 M2071 M2083 
#5     5     5     5     5     5     5     3     8     5     8     5     7     6     5     7     4     4     2     5     5     6     5     5 
#M2103 P6012 P6013 P6025 P6035 P6038 
#3     3     4     2     4     5 

cols_of_int <- c("Was.SCCAI.completed.",
                 "Bowel.frequency.during.the.day",
                 "Bowel.frequency.during.the.night",
                 "Urgency.of.defecation",
                 "Blood.in.the.stool",
                 "General.well.being.over.the.past.24.hours",
                 "sccai")

head(full_meta_ihmp_metabolomics[,cols_of_int])


#compare with the metadata from the reprocessed data; looks like External.ID is the match
#first subset meta_ihmp for only UC
meta_ihmp_UC <- meta_ihmp[meta_ihmp$Study.Group == "UC",]
dim(full_meta_ihmp_metabolomics); dim(meta_ihmp_UC) 
#[1] 146 490 
#[1] 101  25 why these 45 missing samples from reprocessed samples?

full_meta_ihmp_metabolomics$External.ID %in% rownames(meta_ihmp_UC)
#these are the missing ones
full_meta_ihmp_metabolomics$External.ID[!full_meta_ihmp_metabolomics$External.ID %in% rownames(meta_ihmp_UC)]


full_meta_ihmp_metabolomics_sub <- full_meta_ihmp_metabolomics[full_meta_ihmp_metabolomics$External.ID %in% rownames(meta_ihmp_UC),]

#couple of NAs, remove these as well
full_meta_ihmp_metabolomics_sub <- full_meta_ihmp_metabolomics_sub[!is.na(full_meta_ihmp_metabolomics_sub$sccai),]

hist(full_meta_ihmp_metabolomics_sub$sccai)

rownames(full_meta_ihmp_metabolomics_sub) <- full_meta_ihmp_metabolomics_sub$External.ID

oro_vec_ihmp_UC <- mx_ihmp[rownames(full_meta_ihmp_metabolomics_sub), oro_col]
full_meta_ihmp_metabolomics_sub$oro <- oro_vec_ihmp_UC

plot(full_meta_ihmp_metabolomics_sub$sccai, full_meta_ihmp_metabolomics_sub$oro)
cor.test(full_meta_ihmp_metabolomics_sub$sccai, full_meta_ihmp_metabolomics_sub$oro)

#which subjects have multiple samples and a high sccai score in at least one sample?
# 0 to 19 is the range, only 1 has 12 and next highest is 7, 7, and 6.

#look for those subjects
high_symptom_subjects <- unique(full_meta_ihmp_metabolomics_sub$Participant.ID[full_meta_ihmp_metabolomics_sub$sccai >=6])
#"C3013" "H4042" "M2069" "P6012"
table(full_meta_ihmp_metabolomics_sub$Participant.ID)
#M2069 has 6 samples, P6012 has 3 samples, the other ones have only 1.

meta_M2069 <- full_meta_ihmp_metabolomics_sub[full_meta_ihmp_metabolomics_sub$Participant.ID == "M2069",]
meta_P6012 <- full_meta_ihmp_metabolomics_sub[full_meta_ihmp_metabolomics_sub$Participant.ID == "P6012",]

plot(meta_M2069$sccai, meta_M2069$oro) #has the 12 sccai value
plot(meta_P6012$sccai, meta_P6012$oro) #had the 7 vsccai alue


#not sufficient number of subjects to assess relationship between orotic acid and symptoms


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Franzosa

mx_franz <- as.data.frame(read_tsv("../data/processed_data/FRANZOSA_IBD_2019/mtb.tsv"))
meta_franz <- as.data.frame(read_tsv('../data/processed_data/FRANZOSA_IBD_2019/metadata.tsv'))

dim(mx_franz); dim(meta_franz)

rownames(mx_franz) <- mx_franz$Sample
rownames(meta_franz) <- meta_franz$Sample
mx_franz <- mx_franz[rownames(meta_franz),]
#cbind(rownames(mx_franz), rownames(meta_franz))

colnames(mx_franz)[grep("oro", tolower(colnames(mx_franz)))]
oro_col <- grep("HILIC.neg_Cluster_0203..orotate", colnames(mx_franz))
oro_vec_franz <- mx_franz[,oro_col]


#"4-Hydroxybenzoic.acid" "Alpha-aminobutyric.acid" "Cytosine" "Indolelactic.acid" "Orotic.acid" "Xanthine"  
colnames(mx_franz)[grep("but", tolower(colnames(mx_franz)))] #"Alpha-aminobutyric.acid" "HILIC-neg_Cluster_0032: 2-aminobutyric acid" 
colnames(mx_franz)[grep("ind", tolower(colnames(mx_franz)))]
colnames(mx_franz)[grep("benz", tolower(colnames(mx_franz)))]
colnames(mx_franz)[grep("cyt", tolower(colnames(mx_franz)))] #"Cytosine" annotated as "HILIC-pos_Cluster_0032: cytosine"
colnames(mx_franz)[grep("xant", tolower(colnames(mx_franz)))] #"Xanthine" annotated as "HILIC-neg_Cluster_0187: xanthine"


hmb_col <- grep("HILIC-neg_Cluster_0032: 2-aminobutyric acid" , colnames(mx_franz))
ab_vec_franz <- mx_franz[,hmb_col]

hmb_col <- grep("HILIC-pos_Cluster_0032: cytosine", colnames(mx_franz))
cyt_vec_franz <- mx_franz[,hmb_col]

hmb_col <- grep("HILIC-neg_Cluster_0187: xanthine", colnames(mx_franz))
xan_vec_franz <- mx_franz[,hmb_col]



table(meta_franz$Study.Group)
#CD Control      UC 
#88      56      76 


#repeated measures? No
#table(meta_franz$Subject)

meta_franz$oro_vec_franz <- log10(oro_vec_franz)
meta_franz$ab_vec_franz <- log10(ab_vec_franz)
meta_franz$cyt_vec_franz <- log10(cyt_vec_franz)
meta_franz$xan_vec_franz <- log10(xan_vec_franz)


#orotic acid
#comparing to control is appropriate
meta_franz$Study.Group <- relevel(factor(meta_franz$Study.Group), ref = "Control")
summary(lm(oro_vec_franz ~ meta_franz$Study.Group))

#                         Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                2764.9      262.2  10.545  < 2e-16 ***
#meta_franz$Study.GroupCD   -638.7      335.4  -1.904   0.0582 .  
#meta_franz$Study.GroupUC  -1600.2      345.5  -4.631 6.27e-06 ***
boxplot(log10(oro_vec_franz) ~ meta_franz$Study.Group, las=1)

p.vals_franz<- "<0.0001"


#alpha-aminobutyric acid; higher in UC
summary(lm(ab_vec_franz ~ meta_franz$Study.Group)) 
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                351.46     164.12   2.141 0.033350 *  
#meta_franz$Study.GroupCD   754.08     209.94   3.592 0.000406 ***
#meta_franz$Study.GroupUC    89.18     216.29   0.412 0.680502    
boxplot(ab_vec_franz ~ Study.Group, data = meta_franz)


#cytosine
summary(lm(cyt_vec_franz ~ meta_franz$Study.Group)) #trending but not significant
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                387.54      70.27   5.515 9.87e-08 ***
#meta_franz$Study.GroupCD  -161.59      89.89  -1.798   0.0736 .  
#meta_franz$Study.GroupUC  -171.48      92.61  -1.852   0.0654 . 
boxplot(cyt_vec_franz ~ Study.Group, data = meta_franz)


#xanthine
summary(lm(xan_vec_franz ~ meta_franz$Study.Group)) #very significant with both CD and UC
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)               11220.8      753.8  14.887  < 2e-16 ***
#meta_franz$Study.GroupCD  -2842.7      964.2  -2.948 0.003546 ** 
#meta_franz$Study.GroupUC  -3708.8      993.4  -3.734 0.000241 ***
boxplot(xan_vec_franz ~ Study.Group, data = meta_franz)

franz_p_xanth_cd <- 0.003546
franz_p_xanth_uc <- 0.000241



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#combine the two plots as boxplots for main figure

par(mfrow=c(1,2))
boxplot(log10(oro_vec_ihmp) ~ meta_ihmp$Study.Group, las=1)
boxplot(log10(oro_vec_franz) ~ meta_franz$Study.Group, las=1)


sel_colors2 <- sel_colors[8:6]


#to do only use PRISM for Franzosa? There is an issue with sample ids in here that was 
#corrected in the online version of the manuscript. It does not really matter unless matching the pathways
#1:1 with the metabolite abundances. That is not possible with the sample names from microbiome-metabolome-curated-data
#only regarding UMGC samples e.g. Validation.UMCG1180653 


plotname <- "orotic acid metaanalysis manuscript"

png(paste("../../figures/", plotname, ".png", sep=""), width=6.5*250, height=5*250, res=300)    

set.seed(42)
par(mfrow=c(1,2))
par(mar=(c(2.5,3.5,2,0.5)))

x <- meta_ihmp$oro_vec_ihmp
groups <- meta_ihmp$Study.Group
boxplot(x ~ groups, varwidth=F, las=1, xlab="", ylab="", outline=F, ylim=c(3, 7), border=sel_colors2, 
        horizontal=F, col="white")
mtext("log10(abundance)", 2, line=2)
mtext("Lloyd-Price et al, iHMP\nHILn_QI93__orotate", side=3, line=0, adj=0, cex=0.8, font.main=1)
stripchart(x ~ groups, vertical = T, method = "jitter", jitter= 0.25, add = TRUE, pch=16, col=sel_colors2, cex=0.8)

scale_temp <- 1.1
parusr <- par('usr')
segments(1, parusr[4]/scale_temp, 3, parusr[4]/scale_temp) 
text(2, parusr[4]/(scale_temp-0.02), paste("p", p.vals_ihmp[1], sep=""), cex=0.8, font=1) 


par(mar=(c(2.5,3.5,2,0.5)))


x <- meta_franz$oro_vec_franz
groups <- meta_franz$Study.Group
boxplot(x ~ groups, varwidth=F, las=1, xlab="", ylab="", outline=F, ylim=c(0.5,5.5), border=sel_colors2, 
        horizontal=F, col="white")
mtext("log10(abundance)", 2, line=2)
mtext("Franzosa et al, PRISM+validation\nHILIC.neg_Cluster_0203..orotate", side=3, line=0, adj=0, cex=0.8, font.main=1)
stripchart(x ~ groups, vertical = T, method = "jitter", jitter= 0.25, add = TRUE, pch=16, col=sel_colors2, cex=0.8)

scale_temp <- 1.3
parusr <- par('usr')
segments(1, parusr[4]/scale_temp, 3, parusr[4]/scale_temp) 
text(2, parusr[4]/(scale_temp-0.04), paste("p", p.vals_franz[1], sep=""), cex=0.8, font=1) 


dev.off()



#-------------------------------------------------------------------------------
#make the plot for xanthine for supplementary data


plotname <- "metaanalysis manuscript xanthine"

png(paste("../../figures/", plotname, ".png", sep=""), width=7*250, height=5*250, res=300)    

set.seed(42)
par(mfrow=c(1,2))
par(mar=(c(2.5,3.5,2,1)))

x <- meta_ihmp$xan_vec_ihmp
groups <- meta_ihmp$Study.Group
boxplot(x ~ groups, varwidth=F, las=1, xlab="", ylab="", outline=F, ylim=c(3, 8), border=sel_colors2, 
        horizontal=F, col="white")
mtext("log10(abundance)", 2, line=2)
mtext("Lloyd-Price et al, iHMP\nHILn_QI130__xanthine", side=3, line=0, adj=0, cex=0.8, font.main=1)
stripchart(x ~ groups, vertical = T, method = "jitter", jitter= 0.25, add = TRUE, pch=16, col=sel_colors2, cex=0.8)

scale_temp <- 1.1
parusr <- par('usr')
segments(1, parusr[4]/scale_temp, 3, parusr[4]/scale_temp) 
text(2, parusr[4]/(scale_temp-0.02), paste("p=", round(ihmp_p_xanth_uc, digits=3), sep=""), cex=0.8, font=1) 

scale_temp <- 1.05
parusr <- par('usr')
segments(1, parusr[4]/scale_temp, 2, parusr[4]/scale_temp) 
text(1.5, parusr[4]/(scale_temp-0.02), paste("p=", round(ihmp_p_xanth_cd, digits=3), sep=""), cex=0.8, font=1) 



par(mar=(c(2.5,3.5,2,0.5)))
x <- meta_franz$xan_vec_franz
groups <- meta_franz$Study.Group
boxplot(x ~ groups, varwidth=F, las=1, xlab="", ylab="", outline=F, ylim=c(2,5), border=sel_colors2, 
        horizontal=F, col="white")
mtext("log10(abundance)", 2, line=2.5)
mtext("Franzosa et al, PRISM+validation\nHILIC-neg_Cluster_0187: xanthine", side=3, line=0, adj=0, cex=0.8, font.main=1)
stripchart(x ~ groups, vertical = T, method = "jitter", jitter= 0.25, add = TRUE, pch=16, col=sel_colors2, cex=0.8)

scale_temp <- 1.1
parusr <- par('usr')
segments(1, parusr[4]/scale_temp, 3, parusr[4]/scale_temp) 
text(2, parusr[4]/(scale_temp-0.02), paste("p=", round(franz_p_xanth_uc, digits=5), sep=""), cex=0.8, font=1) 

scale_temp <- 1.05
parusr <- par('usr')
segments(1, parusr[4]/scale_temp, 2, parusr[4]/scale_temp) 
text(1.5, parusr[4]/(scale_temp-0.02), paste("p=", round(franz_p_xanth_cd, digits=3), sep=""), cex=0.8, font=1) 


dev.off()




#-------------------------------------------------------------------------------

# Create data frame
data <- data.frame(
  group = c("non IBD / control", "non IBD / control",
            "Crohn's disease (CD)", "Crohn's disease (CD)",
            "Ulcerative Colitis (UC)", "Ulcerative Colitis (UC)"),
  study = c("Lloyd-Price et al, iHMP", "Franzosa et al, PRISM+validation",
            "Lloyd-Price et al, iHMP", "Franzosa et al, PRISM+validation",
            "Lloyd-Price et al, iHMP", "Franzosa et al, PRISM+validation"),
  size = c(26, 56, 49, 88, 30, 76),
  d = c(104, NA, 177, NA, 101, NA)
)

# Create labels for bars
data <- data %>%
  mutate(label = ifelse(!is.na(d), 
                        paste0("n=", size, "\n(d=", d, ")"),
                        paste0("n=", size)))

# Define color mapping for disease groups
disease_colors <- c(
  "non IBD / control" = "#666666",
  "Crohn's disease (CD)" = "#A6761D",
  "Ulcerative Colitis (UC)" = "#E6AB02"
)

# Create plot for Lloyd-Price et al, iHMP (left panel)
plot_lloyd <- data %>%
  filter(study == "Lloyd-Price et al, iHMP") %>%
  ggplot(aes(x = group, y = size, fill = group)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = label), 
            vjust = -0.2,
            size = 3.5) +
  scale_fill_manual(values = disease_colors) +
  labs(title = "Lloyd-Price et al, iHMP",
       x = "",
       y = "Cohort Size (n)") +
  ylim(0, max(data$size) * 1.15) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black")
  )

# Create plot for Franzosa et al, PRISM+validation (right panel)
plot_franzosa <- data %>%
  filter(study == "Franzosa et al, PRISM+validation") %>%
  ggplot(aes(x = group, y = size, fill = group)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = label), 
            vjust = -0.5,
            size = 3.5) +
  scale_fill_manual(values = disease_colors) +
  labs(title = "Franzosa et al, PRISM+validation",
       x = "",
       y = "Cohort Size (n)") +
  ylim(0, max(data$size) * 1.15) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black")
  )

# Combine plots side by side
combined_plot <- plot_lloyd + plot_franzosa +
  plot_layout(ncol = 2) &
  plot_annotation(
    title = "Comparison of Cohort Sizes by Disease Group",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

print(combined_plot)


ggsave(
  filename = "../../figures/cohort_sizes orotic acid metaanalysis.png",
  plot = combined_plot,
  width = 6.5,      # adjust as needed (inches)
  height = 4,     # adjust as needed (inches)
  dpi = 300
)
