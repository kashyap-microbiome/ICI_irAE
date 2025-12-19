require(GUniFrac)
require(openxlsx)
require(tidyverse)
require(phyloseq)
require(data.table)
require(vegan)
require(MatrixGenerics)
require(ggplot2)


#setwd to where this file is
file_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file_dir) 


source("helper functions/zicoseq_volcano_plot_functional.R")
source("helper functions/zicoseq_heatmap_functional.R")
source("helper functions/MannWhitney_function.R")
source("helper functions/make_boxplot_figure.R")


#-------------------------------------------------------------------------------

rmarkdown::render("1_ICI_cohort_description.Rmd")

filename <- "../data/Supplementary Table 1.xlsx"
getSheetNames(filename)


#-------------------------------------------------------------------------------
#load KEGG mapping file to update names of ECs and KOs

kegg_ann_df <- as.data.frame(read.xlsx("../data/Kegg annotations.xlsx", sheet=1, rowNames = F, colNames = F, startRow = 1))
colnames(kegg_ann_df) <- c("KO", "category", "chemical_action", "extra_info", "name")
kegg_ann_df$name_split <- sapply(kegg_ann_df$name, function(x) strsplit(x, "  ")[[1]][2])
kegg_ann_df$name_split_KO <- paste(kegg_ann_df$KO, kegg_ann_df$name_split)
kegg_ann_df$EC_number <- sapply(kegg_ann_df$name, function(x) strsplit(x, "  ")[[1]][1])


#-------------------------------------------------------------------------------
#functional metagenome

func_paths_df <- as.data.frame(read.xlsx(filename, sheet="pathabundance_df", rowNames = T, startRow = 1))

#filter rows with more than 10% 0's
cutoff <- ncol(func_paths_df) * 0.1
func_paths_df <- func_paths_df[which(rowSums(func_paths_df > 0) > cutoff),]

cbind(rownames(clin_meta), colnames(func_paths_df))
dim(clin_meta); dim(func_paths_df)
#[1] 414 177


#-------------------------------------------------------------------------------
#do PERMANOVA for functional pathways
#pathways, corrected for batch

beta_div <- as.matrix(vegdist(t(func_paths_df), method = "bray"))

set.seed(1)
# Permanova - Distance based multivariate analysis of variance
uni_1a <-  GUniFrac::adonis3(beta_div ~ Batch + IMDC, data = clin_meta, permutations = 999)$aov.tab
uni_1b <-  GUniFrac::adonis3(beta_div ~ Batch + hepatitis, data = clin_meta, permutations = 999)$aov.tab
uni_1c <-  GUniFrac::adonis3(beta_div ~ Batch + pneumonitis, data = clin_meta, permutations = 999)$aov.tab
uni_1d <-  GUniFrac::adonis3(beta_div ~ Batch + other_irAE, data = clin_meta, permutations = 999)$aov.tab
uni_4 <-  GUniFrac::adonis3(beta_div ~ Batch + Sex, data = clin_meta, permutations = 999)$aov.tab
uni_5 <-  GUniFrac::adonis3(beta_div ~ Batch + Elixhauser_elix.sum, data = clin_meta, permutations = 999)$aov.tab
uni_6 <-  GUniFrac::adonis3(beta_div ~ Batch + Age, data = clin_meta, permutations = 999)$aov.tab
uni_7 <-  GUniFrac::adonis3(beta_div ~ Batch + abx_last_month, data = clin_meta, permutations = 999)$aov.tab
uni_8 <-  GUniFrac::adonis3(beta_div ~ Batch + BMI, data = clin_meta, permutations = 999)$aov.tab
uni_9 <-  GUniFrac::adonis3(beta_div ~ Batch + Bristol_score_imputed, data = clin_meta, permutations = 999)$aov.tab


beta_res_list <- list(uni_1a, uni_1b, uni_1c, uni_1d, uni_4, uni_5, uni_6, uni_7, uni_8, uni_9)

R2_vals <- unlist(lapply(beta_res_list, function(x) x[,"R2"][2]))
p_vals <- unlist(lapply(beta_res_list, function(x) x[,"Pr(>F)"][2]))
beta_res_df <- as.data.frame(cbind(R2_vals, p_vals))
rownames(beta_res_df) <- as.character(lapply(beta_res_list, function(x) rownames(x)[2]))

beta_res_df <- beta_res_df[order(beta_res_df$R2_vals, decreasing = F),]


#varying color for p values
#n of group size
par(mar=c(4,14,1,1))
barplot(beta_res_df$R2_vals, horiz = T, las=1, col = beta_res_df$p_vals < 0.05, names.arg = rownames(beta_res_df))
mtext("R2 % variance explained", line=2.5, side=1)


filename <- "../figures/func barplot pathways Bray manuscript.png"
png(filename, width = 8*330, height = 5*330, res = 330)
par(mar=c(4,14,1,1))
barplot(beta_res_df$R2_vals, horiz = T, las=1, col = beta_res_df$p_vals < 0.05, names.arg = rownames(beta_res_df))
mtext("R2 % variance explained", line=2.5, side=1)
dev.off()


write.csv(beta_res_df, file = "../data/beta_res_df bray mgx pathways batch cor.csv")


#-------------------------------------------------------------------------------
#adjust for other significant variables

set.seed(1)
# Permanova - Distance based multivariate analysis of variance
uni_1a <-  GUniFrac::adonis3(beta_div ~ Batch + Age + Elixhauser_elix.sum + IMDC, data = clin_meta, permutations = 999)$aov.tab
uni_1b <-  GUniFrac::adonis3(beta_div ~ Batch + Age + Elixhauser_elix.sum + hepatitis, data = clin_meta, permutations = 999)$aov.tab
uni_1c <-  GUniFrac::adonis3(beta_div ~ Batch + Age + Elixhauser_elix.sum + pneumonitis, data = clin_meta, permutations = 999)$aov.tab
uni_1d <-  GUniFrac::adonis3(beta_div ~ Batch + Age + Elixhauser_elix.sum + other_irAE, data = clin_meta, permutations = 999)$aov.tab

beta_res_list <- list(uni_1a, uni_1b, uni_1c, uni_1d)

R2_vals <- unlist(lapply(beta_res_list, function(x) x[,"R2"][4]))
p_vals <- unlist(lapply(beta_res_list, function(x) x[,"Pr(>F)"][4]))
beta_res_df <- as.data.frame(cbind(R2_vals, p_vals))
rownames(beta_res_df) <- as.character(lapply(beta_res_list, function(x) rownames(x)[4]))

beta_res_df <- beta_res_df[order(beta_res_df$R2_vals, decreasing = F),]

#varying color for p values
#n of group size
par(mar=c(4,14,1,1))
barplot(beta_res_df$R2_vals, horiz = T, las=1, col = beta_res_df$p_vals < 0.05, names.arg = rownames(beta_res_df))
mtext("R2 % variance explained", line=2.5, side=1)

write.csv(beta_res_df, file = "../data/beta_res_df_adj bray mgx pathways batch Age Elix Bristol cor.csv")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#IMDC
#zicoseq pathways

data_df <- as.matrix(func_paths_df)
meta <- clin_meta

set.seed(42)
ZicoSeq.obj <- ZicoSeq(meta.dat = meta, feature.dat = data_df, 
                       grp.name = "IMDC", feature.dat.type = "other",
                       adj.name = c('Age', 'Elixhauser_elix.sum'), 
                       prev.filter = 0.2, mean.abund.filter = 0,  
                       #max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,  #99 for testing, 999 for production
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

ZicoSeq.obj.path.IMDC <- ZicoSeq.obj

zico.obj <- ZicoSeq.obj.path.IMDC

grp.name = 'IMDC' # used for extract coefficient
p.cutoff = 0.05 # design our expected significance level
top.n = 8 # For example we want to show significant p.adj.fdr with top 5 effect size(R2)
p_value_class <- "raw" #or "adj"

#get FC
res_df_IMDC_path <- MannWhitney_function(tax_df = data_df, groups = "IMDC", meta = meta)
head(res_df_IMDC_path)
#inspecting
head(res_df_IMDC_path[order(res_df_IMDC_path$pvals, decreasing = F),])

output_plot <- zicoseq_volcano_plot_func(zico.obj, grp.name, p.cutoff, top.n, p_value_class, 
                                         log2FC = res_df_IMDC_path$log2FC, update_names="pathways", save_plot=TRUE)
output_plot$output_plot
output_plot$res_df[output_plot$res_df$sig == "sig",]

dim(output_plot$res_df[output_plot$res_df$sig == "sig",])
#13 pathways

#Zicoseq results for pyrimidine, 3/13 pathways
#                    otu                                                  p.raw p.adj.fdr           R2       coef          R2.2 -log10(p-value) 
#236  PWY-6545: pyrimidine deoxyribonucleotides de novo biosynthesis III  0.034 0.9981310 0.0128180230  2.4128494  0.0128180230   1.468521   0.000812464
#278  PWY-7184: pyrimidine deoxyribonucleotides de novo biosynthesis I    0.032 0.9981310 0.0134080973  2.5425906  0.0134080973   1.494850   0.000812464
#287  PWY-7210: pyrimidine deoxyribonucleotides biosynthesis from CTP     0.028 0.9981310 0.0169048151  2.4879663  0.0169048151   1.552842   0.000812464



#Mann Whitney results, 2 significant
res_df_IMDC_path[grep("pyrimidine", rownames(res_df_IMDC_path)),]

#Mann Whitney results for FCs                                                               pvals  stats adj.pvals      IMDC0      IMDC1      log2FC
#PWY-6545: pyrimidine deoxyribonucleotides de novo biosynthesis III                       0.04705764  910.0 0.9744686  34.488597  75.498325 -1.13032518
#PWY-7184: pyrimidine deoxyribonucleotides de novo biosynthesis I                         0.04705764  910.0 0.9744686  36.540346  80.745675 -1.14389469


sel_colors_full_names <- brewer.pal(n = 8, "Dark2")

#barplots for these pathways
data_vec <- data_df["PWY-6545: pyrimidine deoxyribonucleotides de novo biosynthesis III",]
plotname <- "Pyrimidine pathway III"
make_boxplot_figure(data_vec, clin_meta$IMDC, plottitle = plotname, ylable = "gene counts", 
                    project = "../figures/func IMDC", left_margin=4, add_significance=T, q_or_p = "p",
                    q.vals=0.034, file_type="png",  symbol=T)

data_vec <- data_df["PWY-7184: pyrimidine deoxyribonucleotides de novo biosynthesis I",]
plotname <- "Pyrimidine pathway I"
make_boxplot_figure(data_vec, clin_meta$IMDC, plottitle = plotname, ylable = "gene counts", 
                      project = "../figures/func IMDC", left_margin=4, add_significance=T, q_or_p = "p",
                      q.vals=0.032, file_type="png",  symbol=T)


data_vec <- data_df["PWY-7210: pyrimidine deoxyribonucleotides biosynthesis from CTP",]
plotname <- "Pyrimidine biosynthesis from CPT"
make_boxplot_figure(data_vec, clin_meta$IMDC, plottitle = plotname, ylable = "gene counts", 
                    project = "../figures/func IMDC", left_margin=4, add_significance=T, q_or_p = "p",
                    q.vals=0.028, file_type="png",  symbol=T)


#sum up to check;
data_vec1 <- data_df["PWY-6545: pyrimidine deoxyribonucleotides de novo biosynthesis III",]
data_vec2 <- data_df["PWY-7184: pyrimidine deoxyribonucleotides de novo biosynthesis I",]
data_vec3 <- data_df["PWY-7210: pyrimidine deoxyribonucleotides biosynthesis from CTP",]
#these are very correlated
cor.test(data_vec1, data_vec3)
cor.test(data_vec1, data_vec2)
cor.test(data_vec2, data_vec3)


#-------------------------------------------------------------------------------
#function for boxplot cutting

p_val <- 0.028  # your p value
p_label <- "*"

data_vec <- data_df["PWY-7210: pyrimidine deoxyribonucleotides biosynthesis from CTP",]

# Set your y-axis label and title here
y_label <- "pathway counts"  # Change this to whatever you want
plot_title <- "PWY-7210:\npyrimidine deoxyribonucleotides\nbiosynthesis from CTP"  # Change this to whatever you want


plot_data <- data.frame(
  value = data_vec,
  group = factor(clin_meta$IMDC,
                 levels = c(0, 1),
                 labels = c("no IMDC", "IMDC"))
)

y_max <- max(plot_data$value, na.rm = TRUE)


set.seed(42)
p <- ggplot(plot_data, aes(x = group, y = value, color = group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, 
               fill = "white", linewidth = 1.2) +  # white fill, thicker lines
  geom_jitter(width = 0.15, size = 6, alpha = 0.7) +
  scale_color_manual(values = c("no IMDC" = "darkgrey",
                                "IMDC"   = "#0072B2")) +
  ## significance bar: horizontal line only
  annotate("segment",
           x = 1, xend = 2,
           y = y_max * 1.05, yend = y_max * 1.05,
           size = 0.8, colour = "black") +
  annotate("text",
           x = 1.5, y = y_max * 1.07,
           label = "*", size = 9, colour = "black") +
  coord_cartesian(ylim = c(min(plot_data$value, na.rm = TRUE),
                           y_max * 1.1)) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(colour = "black", size = 20),
    axis.text = element_text(colour = "black", size = 20),
    axis.title = element_text(colour = "black", size = 20),
    plot.title = element_text(colour = "black", size = 20)  # same size as base
  ) +
  labs(x = "", y = y_label, title = plot_title)

print(p)
           
ggsave("../figures/func IMDC PWY-7210.png", p, width = 5, height = 8)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#zicoseq pathways hepatitis

data_df <- as.matrix(func_paths_df)
meta <- clin_meta

set.seed(42)
ZicoSeq.obj <- ZicoSeq(meta.dat = meta, feature.dat = data_df, 
                       grp.name = "hepatitis", feature.dat.type = "other", 
                       adj.name = c('Age', 'Elixhauser_elix.sum'), 
                       prev.filter = 0.2, mean.abund.filter = 0,  
                       #max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,  #99 for testing, 999 for production
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

ZicoSeq.obj.path.hepatitis <- ZicoSeq.obj
zico.obj <- ZicoSeq.obj.path.hepatitis

grp.name = 'hepatitis' # used for extract coefficient
p.cutoff = 0.01 # design our expected significance level
top.n = 5 # For example we want to show significant p.adj.fdr with top 5 effect size(R2)
p_value_class <- "raw" #or "adj"

#get FC
res_df_hepatitis_path <- MannWhitney_function(tax_df = data_df, groups = "hepatitis", meta = meta)

output_plot <- zicoseq_volcano_plot_func(zico.obj, grp.name, p.cutoff, top.n, p_value_class, 
                                         res_df_hepatitis_path$log2FC, update_names="pathways", save_plot=T)
output_plot$output_plot


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#zicoseq pathways pneumonitis

data_df <- as.matrix(func_paths_df)
meta <- clin_meta

set.seed(42)
ZicoSeq.obj <- ZicoSeq(meta.dat = meta, feature.dat = data_df, 
                       grp.name = "pneumonitis", feature.dat.type = "other", 
                       adj.name = c('Age', 'Elixhauser_elix.sum'), 
                       prev.filter = 0.2, mean.abund.filter = 0,  
                       #max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,  #99 for testing, 999 for production
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

ZicoSeq.obj.path.pneumonitis <- ZicoSeq.obj

zico.obj <- ZicoSeq.obj.path.pneumonitis

grp.name = 'pneumonitis' # used for extract coefficient
p.cutoff = 0.01 # design our expected significance level
top.n = 5 # For example we want to show significant p.adj.fdr with top 5 effect size(R2)
p_value_class <- "raw" #or "adj"

#get FC
res_df_pneumonitis_path <- MannWhitney_function(tax_df = data_df, groups = "pneumonitis", meta = meta)

output_plot <- zicoseq_volcano_plot_func(zico.obj, grp.name, p.cutoff, top.n, p_value_class, 
                                         res_df_pneumonitis_path$log2FC, update_names="pathways", save_plot=T)
output_plot$output_plot



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#zicoseq pathways other irAE

data_df <- as.matrix(func_paths_df)
meta <- clin_meta

set.seed(42)
ZicoSeq.obj <- ZicoSeq(meta.dat = meta, feature.dat = data_df, 
                       grp.name = "other_irAE", feature.dat.type = "other", 
                       adj.name = c('Age', 'Elixhauser_elix.sum'), 
                       prev.filter = 0.2, mean.abund.filter = 0,  
                       #max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 999,  strata = NULL,  #99 for testing, 999 for production
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = T)

ZicoSeq.obj.path.other_irAE <- ZicoSeq.obj

zico.obj <- ZicoSeq.obj.path.other_irAE

grp.name = 'other_irAE' # used for extract coefficient
p.cutoff = 0.01 # design our expected significance level
top.n = 5 # For example we want to show significant p.adj.fdr with top 5 effect size(R2)
p_value_class <- "raw" #or "adj"

#get FC
res_df_other_irAE_path <- MannWhitney_function(tax_df = data_df, groups = "other_irAE", meta = meta)

output_plot <- zicoseq_volcano_plot_func(zico.obj, grp.name, p.cutoff, top.n, p_value_class, 
                                         res_df_other_irAE_path$log2FC, update_names="pathways", save_plot=T)
output_plot$output_plot


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#combined heatmaps Pathways

#add option to filter by (relative effect size as well)

zico_list <- list(ZicoSeq.obj.path.IMDC, ZicoSeq.obj.path.hepatitis, 
                  ZicoSeq.obj.path.pneumonitis, ZicoSeq.obj.path.other_irAE)
grp.names <- unlist(lapply(zico_list, function(x) rownames(x$coef.list[[1]])[nrow(x$coef.list[[1]])]))
grp.labels <- grp.names
sign_cutoffs <- c(0.01, 0.005, 0.001) #for stars in plots
p_or_q <- "p" #"q"

heatmap_list <- zicoseq_heatmap_func(zico_list, grp.names, grp.labels, sign_cutoffs, single_cutoff=F, 
                                     p_or_q, rotate=T, column_name_rot = 40, left_margin=8)
heatmap_list$heatmap_plot

filename <- "../figures/functional pathways ICI heatmap irAEs.png"
png(filename, width = 18*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()



#----------------------------------
#using FDR cutoffs only pneumonitis pathways show up

sign_cutoffs <- c(0.1, 0.05, 0.01) #for stars in plots
p_or_q <- "q" #"q"

heatmap_list <- zicoseq_heatmap_func(zico_list, grp.names, grp.labels, sign_cutoffs, single_cutoff=F, 
                                     p_or_q, rotate=T, column_name_rot = 40, left_margin=8)
heatmap_list$heatmap_plot

filename <- "../figures/functional pathways ICI heatmap irAEs fdr.png"
png(filename, width = 10*330, height = 6*330, res = 330)
heatmap_list$heatmap_plot
dev.off()


