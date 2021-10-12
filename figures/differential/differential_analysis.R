library(Seurat)
library(ribor)
library(reshape2)
library(edgeR)
library(smatr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(dplyr)

source('../qc/Ribo_Summary_Function.R')
source('../qc/rename_experiments.R')

human_ribo_file = '../../../../itp/human-itp_v4.ribo'
mouse_ribo_file = '../../../mouse-itp_v5.ribo'

mouse_rnaseq_count_file = '../../cds_counts.csv.gz'

human = Ribo(human_ribo_file, rename_default)
mm    = Ribo(mouse_ribo_file, rename = rename_default)

################################################################################

BURNT_ORANGE = "#bf5700"
UT_BLUE      = "#005f86"

MOUSE_MIN_LENGTH = 29
MOUSE_MAX_LENGTH = 35

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 9

PDF_resolution = 600
FIGURE_FONT    = "helvetica"


ribo_orange = rgb(228,88,10 , maxColorValue = 255)
rna_blue   = rgb(55,135,192, maxColorValue = 255)

################################################################################

get_output_file_path = function(file_name, output_folder = "pdf"){
  this_path = paste( output_folder, file_name, sep = "/"  )
  return(this_path)
}

save_plot_pdf = function(filename, this_plot, width = NA, height = NA){
  this_file = get_output_file_path(filename)
  print(this_file)
  ggsave(this_file, 
         plot   = this_plot, 
         device = cairo_pdf, 
         width  = width,
         height = height,
         dpi    = PDF_resolution )
  
}

################################################################################

plot_pairwise_relationships = function (counts_w, 
                                        id1, id2, 
                                        xlab    = "", 
                                        ylab    = "", 
                                        num_bin = 52, 
                                        xrange  = 100000, 
                                        yrange  = 100000  ) { 
  
  sp = ggscatter(counts_w, x = id1, y = id2,
                 #                add = "reg.line", conf.int = FALSE,     
                 #                add.params = list(color = "blue", size = 0.5),
                 font.family = "Helvetica", 
                 size        = 0.2,
                 color       = "gray", 
                 alpha       = 0.3, 
                 ggtheme     = theme_bw()) 
  
  formatted =   sp +   
    scale_x_log10(labels = scales::label_number_si(), limits = c(0.3, xrange)) +   
    scale_y_log10(labels = scales::label_number_si(), limits = c(0.3, yrange)) + 
    labs (x=xlab, y = ylab) +
    stat_cor(method        = "spearman", 
             aes(label     = ..r.label..), 
             cor.coef.name = "rho", 
             digits        = 2)  + 
    geom_hex(bins= num_bin, aes(alpha=log10(..count..) ), fill="#bf5700" ) +
    theme( axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           plot.title       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)
    )
  return (formatted)  
}


################################################################################

human_ribo_rc <- get_region_counts(human,
                                   range.lower = MOUSE_MIN_LENGTH,
                                   range.upper = MOUSE_MAX_LENGTH,
                                   length      = TRUE,
                                   transcript  = FALSE,
                                   tidy        = FALSE,
                                   alias       = TRUE,
                                   region      = c("CDS"), 
                                   compact     = FALSE)

human_rcw = dcast(human_ribo_rc, transcript ~ experiment)  

################################################################################


# mouse_exp_name_mapper = make_name_mapper( mouse_ribo@experiments  )

#experiment_name_mapper = mouse_exp_name_mapper

human_exp_name_mapper = make_name_mapper( human@experiments  )

human_rename_experiments = function(experiment_name){
  return(human_exp_name_mapper[[experiment_name]])
}

human_rename_experiments = Vectorize(human_rename_experiments, USE.NAMES = FALSE)

colnames(human_rcw) = c( colnames(human_rcw)[1], human_rename_experiments( colnames(human_rcw)[-1]  )  )

human_experiment_names = colnames(human_rcw)[-1]

################################################################################
#### Replicate Clustering By Spearman
####  H E A T M A P

## Replicate clustering by Spearman
## Input is read counts; first column is gene names
replicate_clustering_spearman = function (rcw, 
                                          cpm_threshold = 1, 
                                          breaks_manual = seq(0,1,.01),  
                                          filter        = FALSE, 
                                          clustering    = T){ 
  if (filter) { 
    # Define set of expressed genes before calculating spearman
    expressed = rowSums( cpm (rcw[,-1]) > cpm_threshold) > (dim(rcw)[2] / 3)
    cor_matrix = cor(rcw[expressed,-1], method = "spearman")
  }
  else { 
    cor_matrix = cor(rcw[,-1], method = "spearman")
  }
  
  heatmap_correlation = pheatmap (
    cor_matrix, 
    color        =  colorRampPalette(brewer.pal( 9 ,"Greens"))(length(breaks_manual)),
    #color       =  colorRampPalette( rev(grDevices::heat.colors(100) )  )(length(breaks_manual)),
    cluster_rows = clustering, 
    cluster_cols = F,
    breaks       = breaks_manual,
    fontsize          = 7,
    fontsize_col      = 7,
    fontsize_row      = 7)
  print(heatmap_correlation)
  return(list(cor_matrix = cor_matrix, heatmap_correlation = heatmap_correlation) )
}

################################################################################

cpm_threshold = 1
detected      = rowSums( cpm (human_rcw[,-1]) > cpm_threshold) > (dim(human_rcw)[2] / 3)






plot(means_sd_cpm$Mean_log2CPM_100, sqrt(means_sd_cpm$Sd_log2CPM_100), pch=19, cex = 0.2)
points(means_sd_cpm$Mean_log2CPM_Monosome, sqrt(means_sd_cpm$Sd_log2CPM_Monosome), pch=19, cex = 0.2, col = "red")
### To be continued with Can
### We'll fir LOESS here
itp_loess = loess( sqrt(means_sd_cpm$Sd_log2CPM_100) ~ means_sd_cpm$Mean_log2CPM_100 )
predict (itp_loess, se = T )
lines ( itp_loess)
#####

average_correlation = cor.test(means_sd_cpm$Mean_log2CPM_100, means_sd_cpm$Mean_log2CPM_Monosome, 
                               method = "spearman")
# plot(means_sd_cpm$Mean_log2CPM_100, means_sd_cpm$Mean_log2CPM_Monosome , pch = 19, cex = .5) 

## Spearman Correction -> True correlation Coefficient
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005206

# We might want to think of how to estimate reliability better. 
# I am going to use the highest observed correlation coefficient (Change: I pciked 0.94 from mnosome data)

r_true = average_correlation$estimate / 0.94

# We can use the SMA residuals to plot against txn length, GC-content, etc
s1 = sma (Mean_log2CPM_100  ~ Mean_log2CPM_Monosome, 
          data = means_sd_cpm) 
plot(s1, pch = 19, cex = 0.3, las = 1, xlab = "Conventional\n Mean log2 CPM", ylab = "PAC-ITPl\n Mean log2 CPM")
sma_res = residuals(s1)


################################################################################
################################################################################
################################################################################
################################################################################
###                       M O U S E   F I G U R E S                          ###

min_len = MOUSE_MIN_LENGTH
max_len = MOUSE_MAX_LENGTH

mouse_exp_name_mapper = make_name_mapper( mm@experiments  )


experiment_name_mapper = mouse_exp_name_mapper

rename_experiments = function(experiment_name){
  return(experiment_name_mapper[[experiment_name]])
}

rename_experiments = Vectorize(rename_experiments, USE.NAMES = FALSE)


experiments_to_use = experiments(mm)


ribo_rc <- get_region_counts(mm,
                             experiment = experiments_to_use,
                             range.lower = min_len,
                             range.upper = max_len,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy = F,
                             alias       = TRUE,
                             region      = c("CDS"), 
                             compact = F)

ribo_rc$experiment = rename_experiments(ribo_rc$experiment)

rcw = dcast(ribo_rc, transcript ~ experiment)  
rcw_2 = dcast(ribo_rc, transcript ~ experiment)
 
# GV:    3, 4, 5
# MII:   2, 3, 4
# 1cell: 1, 2, 5
selected_samples = c(2,3,6, 23, 24, 25, 21, 20, 19)
selected_rcw     = rcw[, selected_samples]





## Mouse correlations
rcw_columns_reordered = rcw[, c(1, 17:26, 2:16)]
mouse_cors = replicate_clustering_spearman( rcw_columns_reordered, breaks_manual = seq(0,1,.01), filter = T, clustering = F)



# my matrix subsampling gives very similar results
set.seed(3)
test_normalize_ribo = NormalizeData(as.matrix(rcw_columns_reordered[,-1]), 
                                    scale.factor = 30000, 
                                    normalization.method= "CLR",
                                    margin = 2,
                                    verbose = TRUE)
#test_normalize_ribo = LogNormalize(as.matrix(rcw_columns_reordered[,-1]), scale.factor = 30000, verbose = TRUE)
variables_ribo           = FindVariableFeatures(test_normalize_ribo, selection.method = "vst")
rcw_columns_reordered[variables_ribo$vst.variance.standardized > 5,] 
## Spin1 is a great one.  -> https://journals.biologists.com/dev/article/124/2/493/39566/Spindlin-a-major-maternal-transcript-expressed-in
#Rfpl4 seems very interesting as well
## AU022751 -> Ovary specific

## Zbed3 protein accumulates by 8-cell but translation is over by 4-cell stage
# https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31795-3


feature_color_palette = colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"), space="Lab")




variance_threshold = 5
variance_threshold = 4.125
dim(test_normalize_ribo[variables_ribo$vst.variance.standardized> variance_threshold, ])

distinct_genes_heatmap_from_normalization = 
  pheatmap(test_normalize_ribo[variables_ribo$vst.variance.standardized> variance_threshold, ], 
           labels_row     = strip_extension(rcw_columns_reordered[variables_ribo$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           cutree_rows    = 8,
           fontsize       = 8,
           fontsize_col   = 7,
           #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
           color        = feature_color_palette(100),
           fontsize_row   = 7)

features_from_normalization = rcw_columns_reordered[variables_ribo$vst.variance.standardized > variance_threshold,1]
print(length(features_from_normalization))

# 100 genes
variance_threshold = 3.68
dim(test_normalize_ribo[variables_ribo$vst.variance.standardized> variance_threshold, ])

ribo_distinct_genes_heatmap_from_normalization_main = 
  pheatmap(test_normalize_ribo[variables_ribo$vst.variance.standardized> variance_threshold, ], 
           labels_row     = strip_extension(rcw_columns_reordered[variables_ribo$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           #cutree_rows    = 7,
           fontsize       = 8,
           fontsize_col   = 6,
           treeheight_row = 0, 
           border_color = NA,
           #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
           color        = feature_color_palette(100),
           show_rownames = FALSE,
           fontsize_row   = 6)

features_from_normalization = rcw_columns_reordered[variables_ribo$vst.variance.standardized > variance_threshold,1]
print(length(features_from_normalization))


averaged_log_normalized_ribo = 
  data.frame( GV        = apply(test_normalize_ribo[, 1:5],   1, mean) ,
              MII       = apply(test_normalize_ribo[, 6:10],  1, mean),
              "1cell"   = apply(test_normalize_ribo[, 11:15], 1, mean),
              "2cell"   = apply(test_normalize_ribo[, 16:18], 1, mean),
              "4cell"   = apply(test_normalize_ribo[, 19:21], 1, mean),
              "8cell" = apply(test_normalize_ribo[, 22:25], 1, mean)
              )

variables_ribo_averaged           = FindVariableFeatures(averaged_log_normalized_ribo, selection.method = "vst")

col_names_of_averaged = c("GV", "MII", "1cell", "2cell", "4cell", "8cell")

# 100 genes
variance_threshold = 3.645
#variance_threshold = 2.905
dim(test_normalize_ribo[variables_ribo_averaged$vst.variance.standardized> variance_threshold, ])

averaged_heatmap_top_200 = 
pheatmap(averaged_log_normalized_ribo[variables_ribo_averaged$vst.variance.standardized> variance_threshold, ], 
         labels_row     = strip_extension(rcw_columns_reordered[variables_ribo_averaged$vst.variance.standardized > variance_threshold,1]), 
         cluster_cols   = FALSE, 
         #cutree_rows    = 7,
         fontsize       = 8,
         fontsize_col   = 6,
         treeheight_row = 0, 
         border_color = NA,
         #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
         color        = feature_color_palette(100),
         show_rownames = FALSE,
         fontsize_row   = 6,
         labels_col     = col_names_of_averaged )


# 100 genes
variance_threshold = 4.22
#variance_threshold = 2.905
dim(test_normalize_ribo[variables_ribo_averaged$vst.variance.standardized> variance_threshold, ])

averaged_heatmap_top_100 = 
  pheatmap(averaged_log_normalized_ribo[variables_ribo_averaged$vst.variance.standardized> variance_threshold, ], 
           labels_row     = strip_extension(rcw_columns_reordered[variables_ribo_averaged$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           #cutree_rows    = 7,
           fontsize       = 8,
           fontsize_col   = 6,
           treeheight_row = 0, 
           border_color = NA,
           #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
           color        = feature_color_palette(100),
           show_rownames = FALSE,
           fontsize_row   = 6,
           labels_col     = col_names_of_averaged )




#variance_threshold = 4.125
variance_threshold = 4.8
dim(test_normalize_ribo[variables_ribo_averaged$vst.variance.standardized> variance_threshold, ])
#dim(test_normalize_ribo[variables_ribo_averaged$mvp.dispersion.scaled> variance_threshold, ])

distinct_genes_heatmap_from_normalization_averaged = 
  pheatmap(averaged_log_normalized_ribo[variables_ribo_averaged$vst.variance.standardized> variance_threshold, ], 
           labels_row     = strip_extension(rcw_columns_reordered[variables_ribo_averaged$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           cutree_rows    = 6,
           #cutree_cols    = 3,
           fontsize       = 8,
           fontsize_col   = 7,
           #color         =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
           color          = feature_color_palette(100),
           fontsize_row   = 7,
           treeheight_row = 0, 
           labels_col     = col_names_of_averaged)


################################################################################
########                    M O U S E     RNA-Seq                      #########
mouse_rna_cds = read.csv(mouse_rnaseq_count_file)
## Remove 20210607.RNAseq.4cell.cross.B
mouse_rna_cds = mouse_rna_cds[,-11]

colnames(mouse_rna_cds)  = gsub(colnames(mouse_rna_cds), pattern = ".", fixed = T, replacement = "-")

mouse_rnaseq_old_names = colnames(mouse_rna_cds)[-1]

mouse_rna_name_mapper = make_name_mapper(mouse_rnaseq_old_names)

rename_rnaseq_experiments = function(experiment_name){
  return( paste( mouse_rna_name_mapper[[experiment_name]], "RNAseq" , sep = "-") )
}

# rename_riboseq_experiments = function(experiment_name){
#   return( paste( mouse_exp_name_mapper[[experiment_name]], "Ribo" , sep = "-") )
# }
# 
# colnames(rcw_2) = c( colnames(rcw)[1], rename_riboseq_experiments( colnames(rcw)[-1]))  

rename_rnaseq_experiments = Vectorize(rename_rnaseq_experiments, USE.NAMES = FALSE)
# rename_ribo_experiments = Vectorize(rename_riboseq_experiments, USE.NAMES = FALSE)

new_rnaseq_exp_names = rename_rnaseq_experiments( mouse_rnaseq_old_names )

colnames(mouse_rna_cds) = c(colnames(mouse_rna_cds)[1], new_rnaseq_exp_names)


for (id in 1:length(mouse_rna_cds$transcript)) { 
  mouse_rna_cds$transcript[id] = rename_default(mouse_rna_cds$transcript[id])
}

## Reorder the experiments
mouse_rna_cds = mouse_rna_cds[, c(1,16:23,2:15)]



## RNA-Seq correlation
# mouse_rna_cors = replicate_clustering_spearman( mouse_rna_cds, 
#                                                 breaks_manual = seq(0,1,.01), 
#                                                 filter        = T,
#                                                 clustering    = F )

## RNA-Seq clustering by variable genes
## Based on the above we might want to change this as well
set.seed(3)
test_normalize_rna = LogNormalize(mouse_rna_cds[,-1], scale.factor = 250000, verbose = TRUE)
variables          = FindVariableFeatures(test_normalize_rna, selection.method = "vst")
mouse_rna_cds[variables$vst.variance.standardized > 7, 1] 


# pheatmap(log10 (test_normalize_rna[variables$vst.variance.standardized > 7,-1] + 1 ) , 
#          #labels_row = mouse_rna_cds[variables$vst.variance.standardized > 7,1],
#          cluster_cols   = FALSE, 
#          cutree_rows    = 5,
#          fontsize       = 8,
#          fontsize_col   = 7,
#          #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
#          color        = feature_color_palette(100),
#          labels_row   = unlist( lapply(strsplit( mouse_rna_cds[variables$vst.variance.standardized > 7,1], split = "-" ) , "[[", 1 ) ),
#          fontsize_row   = 7
# )
# 
# # This should give top 200 genes
# variance_threshold = 5.125
# 
# variance_threshold = 4
# 
# dim(mouse_rna_cds[variables$vst.variance.standardized > variance_threshold ,-1])
# 
# pheatmap(log10 (test_normalize_rna[variables$vst.variance.standardized > variance_threshold,-1] + 1 ) , 
#          #labels_row = mouse_rna_cds[variables$vst.variance.standardized > 7,1],
#          cluster_cols   = FALSE, 
#          cutree_rows    = 6,
#          fontsize       = 8,
#          treeheight_row = 0, 
#          show_rownames  = FALSE, 
#          #color         =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
#          color          = feature_color_palette(100),
#          labels_row     = unlist( lapply(strsplit( mouse_rna_cds[variables$vst.variance.standardized > 7,1], split = "-" ) , "[[", 1 ) ),
#          fontsize_row   = 7,
#          fontsize_col   = 6
# )

################################################################################
##########   P R O T E O M I C S    D A T A               ######################




################################################################################
####### Page Layout

generate_blank_plot = function(bg = "white"){
  
  this_background = element_blank()
  
  if(bg != "white"){
    this_background = element_rect(fill = bg)
  }  
  
  df <- data.frame()
  p = ggplot(df) + geom_point() + 
    xlim(0, 4) + ylim(0, 100) + 
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = this_background,
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.x     = element_blank(),
      plot.background  = element_blank(),
      axis.ticks.y     = element_blank()
    ) 
  return(p)
}

separator = generate_blank_plot()

figure_left = plot_grid(distinct_genes_heatmap_from_normalization_averaged [[4]] , 
                          separator, 
                          rel_heights = c(1, 0.5),
                          nrow = 2  )

figure_layout = plot_grid(figure_left,
                          separator,
                          ncol = 2,
                          rel_widths = c(1, 2.8))

save_plot_pdf("diff_figure_layout.pdf", figure_layout, width = unit(7.20, "in"), height = unit(9, "in"))




