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

source('./Ribo_Summary_Function.R')
source('./rename_experiments.R')

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
rna_blue    = rgb(55,135,192, maxColorValue = 255)

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

#human_cors = replicate_clustering_spearman( human_rcw, breaks_manual = seq(0.6,1,.01), clustering = F )

################################################################################
##### Scatter Plots

## Example Pairwise comparisons. These are not the most beautiful but I don't we should try to make them better
sp_10M_1_vs_10M_2 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[1], 
  human_experiment_names[2], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[1], 
  ylab    = human_experiment_names[2])
sp_10M_1_vs_10M_2

sp_10M_2_vs_10M_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[2], 
  human_experiment_names[3], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[2], 
  ylab    = human_experiment_names[3])
sp_10M_2_vs_10M_3


sp_10M_1_vs_10M_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[1], 
  human_experiment_names[3], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[1], 
  ylab    = human_experiment_names[3])
sp_10M_1_vs_10M_3


sp_100_1_vs_100_2 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[4], 
  human_experiment_names[5], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[4], 
  ylab    = human_experiment_names[5])

sp_100_1_vs_100_2


sp_100_2_vs_100_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[5], 
  human_experiment_names[6], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[5], 
  ylab    = human_experiment_names[6])

sp_100_2_vs_100_3

sp_100_1_vs_100_3 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[4], 
  human_experiment_names[6], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[4], 
  ylab    = human_experiment_names[6])

sp_100_1_vs_100_3

sp_100_1_vs_100_3 + theme( legend.position = "none"  )


sp_10M_1_vs_100_1 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[1], 
  human_experiment_names[4], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[1], 
  ylab    = human_experiment_names[4])

sp_10M_1_vs_100_1

sp_10M_2_vs_100_1 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[2], 
  human_experiment_names[4], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[2], 
  ylab    = human_experiment_names[4])

sp_10M_2_vs_100_1

sp_10M_3_vs_100_1 = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[3], 
  human_experiment_names[4], 
  xrange  = 2000, 
  yrange  = 2000, 
  num_bin = 50,
  xlab    = human_experiment_names[3], 
  ylab    = human_experiment_names[4])

sp_10M_3_vs_100_1

sp_grid_alt = 
  plot_grid( sp_10M_1_vs_10M_2 + theme( legend.position = "none"  ),
             sp_10M_2_vs_10M_3 + theme( legend.position = "none"  ),
             sp_10M_1_vs_10M_3 + theme( legend.position = "none"  ),
             sp_100_1_vs_100_2 + theme( legend.position = "none"  ),
             sp_100_2_vs_100_3 + theme( legend.position = "none"  ),
             sp_100_1_vs_100_3 + theme( legend.position = "none"  ),
             sp_10M_1_vs_100_1 + theme( legend.position = "none"  ),
             sp_10M_2_vs_100_1 + theme( legend.position = "none"  ),
             sp_10M_3_vs_100_1 + theme( legend.position = "none"  ),
             align = "hv",
             ncol  = 3)

sp_grid_alt


cpm_threshold = 1
detected      = rowSums( cpm (human_rcw[,-1]) > cpm_threshold) > (dim(human_rcw)[2] / 3)

## Comparing the variance/mean relationship and ITP ~ Monosome similiarity
human_means_cpm = data.frame("Average_100" = apply( cpm(human_rcw[detected,-1])[,4:6]  , 1, mean),
                             "Average_10M" = apply( cpm(human_rcw[detected,-1])[,1:3] , 1, mean)  )

scatter_axis_range = 10000

human_100_vs_10M_means = plot_pairwise_relationships(
  human_means_cpm, 
  "Average_100", 
  "Average_10M", 
  xrange  = scatter_axis_range, 
  yrange  = scatter_axis_range, 
  num_bin = 50,
  xlab    = "100 Average", 
  ylab    = "10M Average")

human_100_vs_10M_means

sp_10M_2_vs_10M_3_new = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[2], 
  human_experiment_names[3], 
  xrange  = scatter_axis_range, 
  yrange  = scatter_axis_range, 
  num_bin = 50,
  xlab    = human_experiment_names[2], 
  ylab    = human_experiment_names[3])

sp_100_1_vs_100_3_new = plot_pairwise_relationships(
  human_rcw, 
  human_experiment_names[4], 
  human_experiment_names[6], 
  xrange  = scatter_axis_range, 
  yrange  = scatter_axis_range, 
  num_bin = 50,
  xlab    = human_experiment_names[4], 
  ylab    = human_experiment_names[6])

sp_grid = 
  plot_grid( 
    sp_10M_2_vs_10M_3_new      +  ggtitle("Conventional") +  theme( legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()),
    sp_100_1_vs_100_3_new      +  ggtitle("PAC-ITP") +  theme( legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()),
    human_100_vs_10M_means     +  ggtitle("PAC-ITP vs Conventional") +  theme( legend.position = "none", plot.title = element_text(hjust = 0.5),  axis.title.x = element_blank(), axis.title.y = element_blank()),
    align = "hv",
    ncol  = 3)

sp_grid

################################################################################

cpm_threshold = 1
detected      = rowSums( cpm (human_rcw[,-1]) > cpm_threshold) > (dim(human_rcw)[2] / 3)

## Comparing the variance/mean relationship and ITP ~ Monosome similiarity
means_sd_cpm = data.frame(Mean_log2CPM_100 = apply( log2(cpm(human_rcw[detected,-1])[,4:6] + 0.5) , 1, mean),
                          Mean_log2CPM_Monosome = apply( log2(cpm(human_rcw[detected,-1])[,1:3] + 0.5 ), 1, mean),
                          Sd_log2CPM_100 = apply( log2(cpm(human_rcw[detected,-1])[,4:6] + 0.5) , 1, sd),
                          Sd_log2CPM_Monosome = apply( log2(cpm(human_rcw[detected,-1])[,1:3] + 0.5 ), 1, sd)
)

means_sd_cpm_filtered = means_sd_cpm[ means_sd_cpm$Mean_log2CPM_100 > 2 & means_sd_cpm$Mean_log2CPM_Monosome > 2, ]

mean_vs_var_100_color_point = rgb(124,203,162 , maxColorValue = 255 )
mean_vs_var_100_color_line  = rgb(4,82,117 , maxColorValue = 255 )
mean_vs_var_10M_color_point = rgb(240,116,110 , maxColorValue = 255 ) 
mean_vs_var_10M_color_line = rgb(110,0,95 , maxColorValue = 255 )

mean_vs_var_colors = c( "100_point" = mean_vs_var_100_color_point,
                        "100_line" = mean_vs_var_100_color_line,
                        "10M_point" = mean_vs_var_10M_color_point,
                        "10M_line" = mean_vs_var_10M_color_line)

mean_vs_variance_human_plot = 
  ggplot( means_sd_cpm_filtered, 
          aes(x= Mean_log2CPM_100, y = sqrt(Sd_log2CPM_100) )  ) +
  
  geom_smooth( method = "loess", se = FALSE , color = mean_vs_var_100_color_line) + 
  geom_smooth(aes(x = Mean_log2CPM_Monosome, y = Sd_log2CPM_Monosome)  , method = "loess", se = FALSE , color = mean_vs_var_10M_color_line) + 
  geom_point(aes(x = Mean_log2CPM_Monosome, y = Sd_log2CPM_Monosome ),size = 0.4, alpha =0.1, color = mean_vs_var_10M_color_point ) + 
  geom_point(alpha =0.1, color = mean_vs_var_100_color_point, size = 0.4  ) + 
  xlab("Mean") + ylab( "Sqrt(Standard Deviation)" ) + ggtitle("Log2(CPM)") + 
  theme_bw() +
  xlim( c(1,13.5) ) 





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
################################################################################
#########         M o u s e     S c a t t e r     P l o t s           ##########

mouse_scatter_axis_range = 1000

sp_1cell_2_5 = plot_pairwise_relationships(
  rcw, 
  "1cell-2", 
  "1cell-5", 
  xrange  = mouse_scatter_axis_range, 
  yrange  = mouse_scatter_axis_range, 
  num_bin = 50,
  xlab    = "1cell-2", 
  ylab    = "1cell-5")

sp_1cell_2_5 + theme( legend.position = "none"  )

sp_mii_3_4 = plot_pairwise_relationships(
  rcw, 
  "MII-3", 
  "MII-4", 
  xrange  = mouse_scatter_axis_range, 
  yrange  = mouse_scatter_axis_range, 
  num_bin = 50,
  xlab    = "MII-3", 
  ylab    = "MII-4")

sp_mii_3_4 + theme( legend.position = "none"  )


sp_gv_3_4 = plot_pairwise_relationships(
  rcw, 
  "GV-3", 
  "GV-4", 
  xrange  = mouse_scatter_axis_range, 
  yrange  = mouse_scatter_axis_range, 
  num_bin = 50,
  xlab    = "GV-3", 
  ylab    = "GV-4")

sp_gv_3_4 + theme( legend.position = "none"  )

mouse_sp_grid = 
  plot_grid( 
    sp_gv_3_4      +  ggtitle("GV") +  theme( legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()),
    sp_mii_3_4    +   ggtitle("MII") +  theme( legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()),
    sp_1cell_2_5  +   ggtitle("1 cell") +  theme( legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()),
    align = "hv",
    ncol  = 3)

mouse_sp_grid

################################################################################

# GV:    3, 4, 5
# MII:   2, 3, 4
# 1cell: 1, 2, 5
selected_samples = c(2,3,6, 23, 24, 25, 21, 20, 19)
selected_rcw     = rcw[, selected_samples]

set.seed(3)
sub_40k = subsampleMatrix( selected_rcw, 40000)
sub_30k = subsampleMatrix( selected_rcw, 30000)
sub_20k = subsampleMatrix( selected_rcw, 20000)
sub_10k = subsampleMatrix( selected_rcw, 10000)
sub_5k  = subsampleMatrix( selected_rcw, 5000)
sub_0   = subsampleMatrix( selected_rcw, 0)

# Plot number of genes detected as a function of UMI count
# We should plot a line graph that shows cell to cell relationship as well the sum across cells. 
genes_detected_df = data.frame( 
  sample_depth = c( colSums(selected_rcw ), colSums(sub_40k), colSums(sub_30k), 
                    colSums(sub_20k), colSums(sub_10k), colSums(sub_5k), colSums(sub_0)), 
  genes_detected = c(   apply (selected_rcw, 2, function(x){sum(x > 0)}), 
                        apply (sub_40k, 2, function(x){sum(x > 0)}), 
                        apply (sub_30k, 2, function(x){sum(x > 0)}), 
                        apply (sub_20k, 2, function(x){sum(x > 0)}), 
                        apply (sub_10k, 2, function(x){sum(x > 0)}),
                        apply (sub_5k, 2, function(x){sum(x > 0)}), 
                        apply (sub_0, 2, function(x){sum(x > 0)})
  ), 
  sample_id =  rep(colnames(selected_rcw), 7)
  
)

## Ggplot will be finalized
umis_vs_genes_detected_plot = 
  ggplot(
    data = genes_detected_df, 
    aes(x = sample_depth, y = genes_detected, group=sample_id) ) +
  geom_line(linetype = "dashed", aes(col = sample_id) )  +  ylim(c( 0, 7000) ) +
  geom_point(size = 0.7, aes(color = sample_id) ) + 
  xlab("Number of CDS Mapping UMIs") + 
  ylab("Genes Detected") + 
  theme_bw() + 
  theme(legend.title      = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  #scale_x_continuous( expand = c(0.05, 0.05) )
  scale_x_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si()) +
  scale_y_continuous( expand = c(0.05, 0.05), minor_breaks = NULL, labels = scales::label_number_si())



## Mouse correlations
rcw_columns_reordered = rcw[, c(1, 17:26, 2:16)]
mouse_cors = replicate_clustering_spearman( rcw_columns_reordered, breaks_manual = seq(0,1,.01), filter = T, clustering = F)

## We can add a few example pairwise comparisons: 
plot_pairwise_relationships(rcw, 
                            "GV-3", 
                            "GV-5", 
                            xrange = 3000, yrange = 3000, num_bin = 50,
                            xlab = "GV-3", ylab = "GV-5")


plot_pairwise_relationships(rcw, 
                            "1cell-3", 
                            "1cell-5", 
                            xrange = 1000, yrange = 1000, num_bin = 50,
                            xlab = "1cell-3", ylab = "1cell-5")

# my matrix subsampling gives very similar results
set.seed(3)
test_normalize = LogNormalize(as.matrix(rcw_columns_reordered[,-1]), scale.factor = 30000, verbose = TRUE)
variables      = FindVariableFeatures(test_normalize, selection.method = "vst")
rcw_columns_reordered[variables$vst.variance.standardized > 5,] 
## Spin1 is a great one.  -> https://journals.biologists.com/dev/article/124/2/493/39566/Spindlin-a-major-maternal-transcript-expressed-in
#Rfpl4 seems very interesting as well
## AU022751 -> Ovary specific

## Zbed3 protein accumulates by 8-cell but translation is over by 4-cell stage
# https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31795-3


feature_color_palette = colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"), space="Lab")




variance_threshold = 4.3
distinct_genes_heatmap_from_normalization = 
  pheatmap(test_normalize[variables$vst.variance.standardized> variance_threshold, ], 
           labels_row     = strip_extension(rcw_columns_reordered[variables$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           cutree_rows    = 8,
           fontsize       = 8,
           fontsize_col   = 7,
           #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
           color        = feature_color_palette(100),
           fontsize_row   = 7)

features_from_normalization = rcw_columns_reordered[variables$vst.variance.standardized > variance_threshold,1]
print(length(features_from_normalization))



#  We decided not to use the following sampling. Instead we use "LogNormalize" above
# which does the sampling (via scale.factor)
#
# set.seed(3)
# sample_all_30k           = subsampleMatrix( rcw_columns_reordered[, -1], 30000)
# sample_all_30k           = LogNormalize(as.matrix(sample_all_30k), scale.factor = 30000, verbose = TRUE)
# subsample_variables      = FindVariableFeatures(sample_all_30k, selection.method = "vst")
# variance_threshold_subsample = 3.86
# 
# 
# feature_color_palette = colorRampPalette(brewer.pal(9,"YlGn")) 
# feature_color_palette = colorRampPalette(c("#3362A5", "dodgerblue3", "dodgerblue2", "dodgerblue1", "deepskyblue", "lightblue", "lightgray", "gold", "orange","firebrick2"), space="Lab")
# 
# distinct_genes_heatmap_from_subsampling = 
# pheatmap(sample_all_30k[subsample_variables$vst.variance.standardized> variance_threshold_subsample, ], 
#          labels_row     = strip_extension( rcw_columns_reordered[subsample_variables$vst.variance.standardized > variance_threshold_subsample,1]), 
#          cluster_cols   = FALSE, 
#          cutree_rows    = 7,
#          fontsize       = 8,
#          fontsize_col   = 7,
#          color        = feature_color_palette(100),
#          fontsize_row   = 7)
# 
# features_from_subsampling = rcw_columns_reordered[subsample_variables$vst.variance.standardized > variance_threshold_subsample,1]
# print(length(features_from_subsampling))
# 
# length(intersect(features_from_normalization, features_from_subsampling) )



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
mouse_rna_cors = replicate_clustering_spearman( mouse_rna_cds, 
                                                breaks_manual = seq(0,1,.01), 
                                                filter        = T,
                                                clustering    = F )

## RNA-Seq clustering by variable genes
## Based on the above we might want to change this as well
set.seed(3)
test_normalize_rna = LogNormalize(mouse_rna_cds[,-1], scale.factor = 250000, verbose = TRUE)
variables          = FindVariableFeatures(test_normalize_rna, selection.method = "vst")
mouse_rna_cds[variables$vst.variance.standardized > 7, 1] 


pheatmap(log10 (mouse_rna_cds[variables$vst.variance.standardized > 7,-1] + 1 ) , 
         #labels_row = mouse_rna_cds[variables$vst.variance.standardized > 7,1],
         cluster_cols   = FALSE, 
         cutree_rows    = 5,
         fontsize       = 8,
         fontsize_col   = 7,
         #color        =  colorRampPalette(brewer.pal( 9 ,"YlGnBu"))(100),
         color        = feature_color_palette(100),
         labels_row   = unlist( lapply(strsplit( mouse_rna_cds[variables$vst.variance.standardized > 7,1], split = "-" ) , "[[", 1 ) ),
         fontsize_row   = 7
)

## Combined data correlation
all_counts = merge(mouse_rna_cds, rcw, by = "transcript")

new_combined_ordered_indeces = c(39:43, 2:5, 
                                 44:48, 6:9,
                                 24:28, 10:13,
                                 29:31, 14:17,
                                 32:34, 18:19,
                                 35:38, 20:23)

all_counts_arranged = all_counts[, new_combined_ordered_indeces]


## Might want to think about which genes are included in the subsampling. 
## It might be good to filter based on expression before subsampling. 
set.seed(3)
allcounts_cds_subsample = subsampleMatrix (all_counts_arranged[,-1], desiredDepth = 28000)

combined_cors = 
  replicate_clustering_spearman( 
    allcounts_cds_subsample, 
    breaks_manual = seq(0,1,.01), 
    filter = T,
    clustering = FALSE)


all_counts_ribo_rna_separated = all_counts[ colnames( all_counts )[c(1:23, 39:48, 24:38)] ]

combined_cors_rna_ribo_separate = 
  replicate_clustering_spearman( 
    all_counts_ribo_rna_separated, 
    breaks_manual = seq(0,1,.01), 
    filter = T,
    clustering = FALSE)

se = function(x) { sd(x) / sqrt (length(x) )}

plot_cpm_across_conditions = function(counts, gene, stages = c("all"), ymax = 4, plot_type = "point") { 
  # counts is a numeric matrix of raw read counts + gene_ids
  # gene to be plotted. Formatted as "Obox2"
  # Stages are 1cell, 2cell, MII, GV, etc
  set.seed(3)
  colnames(counts) = c( "transcript" ,  paste ( colnames(counts)[-1], 'Ribo', sep = "-" ) )
  
  # A monkey patch for the y-axis aesthetics
  ybreaks = waiver()
  
  if(ymax < 4){
    ybreaks = 1:ymax
  }
  
  if (stages[1] == "all" )  {  
    normalizedcounts = NormalizeData(counts[,-1], 
                                     scale = 10000, 
                                     normalization.method = "CLR", 
                                     margin = 2)
  } 
  else { 
    selected_samples = grep(paste(stages, collapse  = "|"), colnames(counts))
    normalizedcounts = NormalizeData(counts[,selected_samples], scale = 10000,
                                     normalization.method = "CLR", 
                                     margin = 2)
  }
  
  tidy_norm        = melt(normalizedcounts[ strip_extension(as.character(counts$transcript)) %in% gene,  ] )
  tidy_norm$stage  = sapply(strsplit(as.character(colnames(normalizedcounts)), split = "-"), "[[", 1) 
  tidy_norm$stage  = factor( tidy_norm$stage, levels = c("GV", "MII", "1cell",  "2cell",  "4cell", "8cell")  )
  tidy_norm$method = sapply(strsplit(as.character(colnames(normalizedcounts)), split = "-"), "[[", 3) 
  
  p =
    ggplot(tidy_norm, 
           aes( x =  strip_extension(as.character(gene)), y = value) ) + 
    # stat_summary(aes( shape = stage, color = method), fun.data = "mean_se", size = 1, 
    #              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), 
    #              show.legend = FALSE, inherit.aes = FALSE) +
    geom_jitter(aes( shape = stage, color = method),
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1),
                size = 2 ) +
    # scale_color_manual(values =  c('#089099', '#d12959')) + 
    theme_bw()+
    theme(
      axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
    ) + 
    scale_y_continuous( limits = c(0, ymax) , breaks = ybreaks ) + 
    scale_color_manual(values = c("Ribo" = ribo_orange, "RNAseq" = rna_blue)) + 
    scale_shape_manual(values = c(4, 16, 3, 6)) + 
    ylab("Normalized Count") + xlab("") 
  
  if(plot_type == "point"){
    return(p)
  }
  else if (plot_type == "barplot"){
    mean_se_tidy = tidy_norm %>% group_by(stage, method) %>%
      summarise(mean = mean(value), se = se(value))
    
    p =
      ggplot(mean_se_tidy, aes(x=stage , y=mean,  fill=method)) + 
      geom_bar(position=position_dodge(), stat="identity") +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                    width=.25, # Width of the error bars
                    position=position_dodge(.9)) + 
      theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
            panel.border     = element_blank(),
            panel.grid       = element_blank(),
            panel.background = element_blank(),
            axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            legend.title     = element_blank(),
            legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            legend.key.size  = unit(0.15, 'in'),
      ) +
      scale_y_continuous( limits = c(0, ymax) , expand = c(0, 0), breaks = ybreaks ) + 
      scale_fill_manual(values = c("Ribo" = ribo_orange, "RNAseq" = rna_blue)) +  
      labs(title = gene) + 
      ylab("Normalized Count")
    
    return(p)
  }
  else { 
    print("Incorrect plot type")}
}


### Higher in MII
Lbr_barplot = 
  plot_cpm_across_conditions(all_counts, "Lbr", stages = c("1cell", "MII"), plot_type = "barplot", ymax = 4 ) 
Lbr_barplot

Lbr_pointplot = 
  plot_cpm_across_conditions(all_counts, "Lbr", 
                             stages    = c("1cell", "MII"), 
                             plot_type = "point", ymax = 4 ) 
Lbr_pointplot


Cept1_barplot = 
  plot_cpm_across_conditions(all_counts, "Cept1", 
                             stages = c("1cell", "MII"), 
                             plot_type = "barplot", ymax = 5 ) 
Cept1_barplot

Cept1_pointplot = 
  plot_cpm_across_conditions(all_counts, "Cept1", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Cept1_pointplot


H1f8_barplot = 
  plot_cpm_across_conditions(all_counts, "H1f8", 
                             stages = c("1cell", "MII"), 
                             plot_type = "barplot", ymax = 8 ) 
H1f8_barplot

H1f8_pointplot = 
  plot_cpm_across_conditions(all_counts, "H1f8", stages = c("1cell", "MII"), plot_type = "point", ymax = 8 ) 
H1f8_pointplot


Taldo1_barplot = 
  plot_cpm_across_conditions(all_counts, "Taldo1", 
                             stages = c("1cell", "MII"), 
                             plot_type = "barplot", ymax = 5 ) 
Taldo1_barplot

Taldo1_pointplot = 
  plot_cpm_across_conditions(all_counts, "Taldo1", stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Taldo1_pointplot

Gdf9_barplot = 
  plot_cpm_across_conditions(all_counts, "Gdf9", 
                             stages = c("1cell", "MII"), 
                             plot_type = "barplot", ymax = 8 ) 
Gdf9_barplot

Gdf9_pointplot = 
  plot_cpm_across_conditions(all_counts, "Gdf9", 
                             stages    = c("1cell", "MII"), 
                             plot_type = "point", 
                             ymax      = 8 ) 
Gdf9_pointplot

Amfr_barplot = 
  plot_cpm_across_conditions(all_counts, "Amfr", 
                             stages = c("1cell", "MII"), 
                             plot_type = "barplot", ymax = 4 ) 
Amfr_barplot

Amfr_pointplot = 
  plot_cpm_across_conditions(all_counts, "Amfr", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Amfr_pointplot


Pdcd6_barplot = 
  plot_cpm_across_conditions(all_counts, "Pdcd6", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 5 ) 
Pdcd6_barplot

Pdcd6_pointplot = 
  plot_cpm_across_conditions(all_counts, "Pdcd6", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Pdcd6_pointplot


Tada2a_barplot = 
  plot_cpm_across_conditions(all_counts, "Tada2a", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 4 ) 
Tada2a_barplot

Tada2a_pointplot = 
  plot_cpm_across_conditions(all_counts, "Tada2a", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Tada2a_pointplot


## Higher in 1cell

Anapc5_barplot = 
  plot_cpm_across_conditions(all_counts, "Anapc5", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 4 ) 
Anapc5_barplot

Anapc5_pointplot = 
  plot_cpm_across_conditions(all_counts, "Anapc5", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Anapc5_pointplot


Spry4_barplot = 
  plot_cpm_across_conditions(all_counts, "Spry4", 
                             stages = c("1cell", "MII"), plot_type = "barplot" , ymax = 4 ) 
Spry4_barplot

Spry4_pointplot = 
  plot_cpm_across_conditions(all_counts, "Spry4", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5) 
Spry4_pointplot


Mapk1_barplot = 
  plot_cpm_across_conditions(all_counts, "Mapk1", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 3 ) 
Mapk1_barplot

Mapk1_pointplot = 
  plot_cpm_across_conditions(all_counts, "Mapk1", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 3 ) 
Mapk1_pointplot


Abcf3_barplot = 
  plot_cpm_across_conditions(all_counts, "Abcf3", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 4 ) 
Abcf3_barplot

Abcf3_pointplot = 
  plot_cpm_across_conditions(all_counts, "Abcf3", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Abcf3_pointplot


Smarca4_barplot = 
  plot_cpm_across_conditions(all_counts, "Smarca4", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 3 ) 
Smarca4_barplot

Smarca4_pointplot = 
  plot_cpm_across_conditions(all_counts, "Smarca4", stages = c("1cell", "MII"), plot_type = "point", ymax = 3 ) 
Smarca4_pointplot


Mapk1_barplot = 
  plot_cpm_across_conditions(all_counts, "Mapk1", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 3 ) 
Mapk1_barplot

Mapk1_pointplot = 
  plot_cpm_across_conditions(all_counts, "Mapk1", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 3 ) 
Mapk1_pointplot


Rybp_barplot = 
  plot_cpm_across_conditions(all_counts, "Rybp", 
                             stages = c("1cell", "MII"), plot_type = "barplot", ymax = 3 ) 
Rybp_barplot

Rybp_pointplot = 
  plot_cpm_across_conditions(all_counts, "Rybp", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 3 ) 
Rybp_pointplot

################################################################################

pointplot_legend = get_legend(Lbr_pointplot)

supp_control_genes_mii_1cell = 
  plot_grid(Lbr_pointplot     + theme(legend.position = "none"),
            Cept1_pointplot   + theme(legend.position = "none"),
            H1f8_pointplot    + theme(legend.position = "none"),
            Taldo1_pointplot  + theme(legend.position = "none"),
            Gdf9_pointplot    + theme(legend.position = "none"),
            Amfr_pointplot    + theme(legend.position = "none"),
            Pdcd6_pointplot   + theme(legend.position = "none"),
            Tada2a_pointplot  + theme(legend.position = "none"),
            Anapc5_pointplot  + theme(legend.position = "none"),
            Spry4_pointplot   + theme(legend.position = "none"),
            Mapk1_pointplot   + theme(legend.position = "none"),
            Abcf3_pointplot   + theme(legend.position = "none"),
            Smarca4_pointplot + theme(legend.position = "none"),
            Rybp_pointplot    + theme(legend.position = "none"), 
            pointplot_legend,
            ncol = 5)


supp_control_genes_mii_1cell

save_plot_pdf("supp_control_genes_mii_1cell_poin_plots.pdf", 
              supp_control_genes_mii_1cell,
              width = 6, height = 5)


supp_barplot_legend = get_legend(Lbr_barplot)

supp_control_genes_mii_1cell_barplots = 
  plot_grid(Lbr_barplot     + theme(legend.position = "none"),
            Cept1_barplot   + theme(legend.position = "none"),
            H1f8_barplot    + theme(legend.position = "none"),
            Taldo1_barplot  + theme(legend.position = "none"),
            Gdf9_barplot    + theme(legend.position = "none"),
            Amfr_barplot    + theme(legend.position = "none"),
            Pdcd6_barplot   + theme(legend.position = "none"),
            Tada2a_barplot  + theme(legend.position = "none"),
            Anapc5_barplot  + theme(legend.position = "none"),
            Spry4_barplot   + theme(legend.position = "none"),
            Mapk1_barplot   + theme(legend.position = "none"),
            Abcf3_barplot   + theme(legend.position = "none"),
            Smarca4_barplot + theme(legend.position = "none"),
            Rybp_barplot    + theme(legend.position = "none"), 
            supp_barplot_legend,
            ncol = 5)


supp_control_genes_mii_1cell_barplots

save_plot_pdf("supp_control_genes_mii_1cell_bar_plots.pdf", 
              supp_control_genes_mii_1cell_barplots,
              width = 6, height = 5)
## Supp Grid

legend_of_genes_accross_stages_supp_figure_bar = get_legend(Gdf9_barplot)


supp_selected_genes_in_mii_1cell_barplot_grid = 
  plot_grid(
    Gdf9_barplot + theme(legend.position = "none") + ylab(""),
    Amfr_barplot + theme(legend.position = "none") + ylab(""),
    Smarca4_barplot + theme(legend.position = "none") + ylab(""),
    Mapk1_barplot + theme(legend.position = "none") + ylab(""),
    Pdcd6_barplot + theme(legend.position = "none") + ylab(""),
    Tada2a_barplot + theme(legend.position = "none") + ylab(""),
    Abcf3_barplot + theme(legend.position = "none") + ylab(""), 
    Rybp_barplot + theme(legend.position = "none") + ylab(""),
    nrow = 2
    
  )


save_plot_pdf("supp_selected_genes_in_mii_1cell_barplot_grid.pdf", supp_selected_genes_in_mii_1cell_barplot_grid,
              width = 3.4, height = 2.9)

################################################################################
#### Combining Plots

legend_of_genes_accross_stages_main_figure_bar = get_legend(Lbr_barplot)

genes_accross_stages_main_figure_bar_pre = 
  plot_grid(  Lbr_barplot + theme(legend.position = "none"),
              Cept1_barplot + theme(legend.position = "none") + ylab(""),
              H1f8_barplot + theme(legend.position = "none") + ylab(""),
              Taldo1_barplot + theme(legend.position = "none") + ylab(""),
              Anapc5_barplot + theme(legend.position = "none") +  ylab(""),
              Spry4_barplot + theme(legend.position = "none") +  ylab(""),
              #Mapk1_barplot + theme(legend.position = "none") +  ylab(""),
              nrow = 1)

genes_accross_stages_main_figure_bar = plot_grid(
  genes_accross_stages_main_figure_bar_pre,
  legend_of_genes_accross_stages_main_figure_bar,
  nrow = 1,
  rel_widths = c(5,0.8)
)

genes_accross_stages_main_figure_bar

save_plot_pdf("genes_accross_bar_main.pdf", genes_accross_stages_main_figure_bar,
              width = 7.2, height = 2)



##### MAPK1 goes to the supp.
## Mapk1_barplot + theme(legend.position = "none") +  ylab(""),

### Selected plots for  

Ppig_point_plot = 
  plot_cpm_across_conditions(all_counts, "Ppig", stages = c("1cell", "MII"), plot_type = "point", ymax = 3 ) 
Ppig_point_plot

Nlrp14_point_plot = 
  plot_cpm_across_conditions(all_counts, "Obox7", stages = c("1cell", "MII"), plot_type = "point", ymax = 3 ) 
Nlrp14_point_plot

################################################################################
# To be deleted
# plot_cpm_across_conditions(all_counts, "Nfrkb", stages = c("1cell", "MII") ) 
# 
# nfkrb_stages_plot = 
# plot_cpm_across_conditions(all_counts, "Nfrkb", stages = c("1cell", "MII") ) 
# 
# taldo1stages_plot = 
#   plot_cpm_across_conditions(all_counts, "Spry4", stages = c("1cell", "MII"), plot_type = "barplot" ) 
# taldo1stages_plot
# 
# 
# save_plot_pdf("nfkrb_stages_plot.pdf", nfkrb_stages_plot, 
#               width = unit(3, "in"), height = unit(2.5, "in") )
# 
# ppig_stages_plot = 
# plot_cpm_across_conditions(all_counts, "Ppig", stages = c("1cell", "MII", "2cell", "4cell") ) 
# 
# nop58_stages_plot = 
# plot_cpm_across_conditions(all_counts, "Nop58", stages = c("1cell", "MII", "2cell", "4cell") )  ## Proteomics
# 
# bottom_legend = get_legend(nop58_stages_plot)
# 
# selected_genes_accross_stages =
#                           plot_grid(  nfkrb_stages_plot + theme(legend.position = "none"), 
#                           ppig_stages_plot  + theme(legend.position = "none"), 
#                           nop58_stages_plot + theme(legend.position = "none"), 
#                           nfkrb_stages_plot + theme(legend.position = "none"), 
#                           bottom_legend,
#                           rel_widths = c(1, 1, 1, 1, 0.5),
#                           nrow = 1 )

################################################################################
#######                P D F    P R O D U C T I O N                      #######



################################################################################
#######                   P A G E     L A Y O U T                         ######

## Moved elsewhere
################################################################################

save_plot_pdf("supp_pairwise_correlation.pdf", sp_grid, width = unit(7.05, "in"), height = unit(2.35, "in"))

save_plot_pdf("mouse_correlation_map.pdf", mouse_cors$heatmap_correlation[[4]], width = unit(4.7, "in"), height = unit(4.7, "in"))

save_plot_pdf("supp_mouse_pairwise_correlation.pdf", mouse_sp_grid, width = unit(7.05, "in"), height = unit(2.35, "in"))

save_plot_pdf("distinct_genes_heatmap_from_normalization.pdf", distinct_genes_heatmap_from_normalization, 
              width = unit(5, "in"), height = unit(8.5, "in") )

save_plot_pdf("mean_vs_variance_human_plot_add_legend.pdf", mean_vs_variance_human_plot, 
              width = unit(3, "in"), height = unit(3, "in") )

save_plot_pdf("combined_cors.pdf", combined_cors$heatmap_correlation[[4]], 
              width = unit(7.2, "in"), height = unit(7.2, "in") )

save_plot_pdf("umis_vs_genes_detected.pdf", umis_vs_genes_detected_plot, 
              width = unit(3.45, "in"), height = unit(2.7, "in") )

combined_cors_rna_ribo_separate

save_plot_pdf("combined_cors_rna_ribo_separate.pdf", combined_cors_rna_ribo_separate$heatmap_correlation[[4]], 
              width = unit(7.2, "in"), height = unit(7.2, "in") )



#############################################################################

Apc = 
  plot_cpm_across_conditions(all_counts, "Apc", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Apc

Kit = 
  plot_cpm_across_conditions(all_counts, "Kit", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Kit

#Myh9 = 
#  plot_cpm_across_conditions(all_counts, "Myh9", 
#                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
#Myh9

Pcm1 = 
  plot_cpm_across_conditions(all_counts, "Pcm1", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Pcm1


Ptk2 = 
  plot_cpm_across_conditions(all_counts, "Ptk2", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Ptk2

# Rock1 = 
#   plot_cpm_across_conditions(all_counts, "Rock1", 
#                              stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
# Rock1


Hnrnpa2b1 = 
  plot_cpm_across_conditions(all_counts, "Hnrnpa2b1", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Hnrnpa2b1

Qk = 
  plot_cpm_across_conditions(all_counts, "Qk", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Qk

Celf1 = 
  plot_cpm_across_conditions(all_counts, "Celf1", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Celf1

Cnot4 = 
  plot_cpm_across_conditions(all_counts, "Cnot4", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Cnot4



Dazl = 
  plot_cpm_across_conditions(all_counts, "Dazl", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Dazl


#Pa2g4 = 
#  plot_cpm_across_conditions(all_counts, "Pa2g4", 
#                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
#Pa2g4

Stau1 = 
  plot_cpm_across_conditions(all_counts, "Stau1", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Stau1


#Ythdc1 = 
#  plot_cpm_across_conditions(all_counts, "Ythdc1", 
#                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
#Ythdc1

Calr = 
  plot_cpm_across_conditions(all_counts, "Calr", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
Calr


# Eprs = 
#   plot_cpm_across_conditions(all_counts, "Eprs", 
#                              stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
# Eprs


# Fxr1 = 
#   plot_cpm_across_conditions(all_counts, "Fxr1", 
#                              stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
# Fxr1


# Larp7 = 
#   plot_cpm_across_conditions(all_counts, "Larp7", 
#                              stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
# Larp7


# Ncl = 
#   plot_cpm_across_conditions(all_counts, "Ncl", 
#                              stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
# Ncl

# 
# Snw1 = 
#   plot_cpm_across_conditions(all_counts, "Snw1", 
#                              stages = c("1cell", "2cell"), plot_type = "point", ymax = 5 ) 
# Snw1


Ythdf2 = 
  plot_cpm_across_conditions(all_counts, "Ythdf2", 
                             stages = c("1cell", "2cell"), plot_type = "point", ymax = 3 ) 
Ythdf2

Camsap1 = 
  plot_cpm_across_conditions(all_counts, "Camsap1", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Camsap1

Cep120 = 
  plot_cpm_across_conditions(all_counts, "Cep120", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Cep120

Numa1 = 
  plot_cpm_across_conditions(all_counts, "Numa1", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 4 ) 
Numa1

Ezh2 = 
  plot_cpm_across_conditions(all_counts, "Ezh2", 
                           stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Ezh2




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

point_legend_1 = get_legend(Apc)
point_legend_2 = get_legend(Ythdf2)

separator = generate_blank_plot()

selected_genes_point_plots = 
  plot_grid( point_legend_1,
             point_legend_2,
             separator,
             Apc       + theme(legend.position = "none"),
             Kit       + theme(legend.position = "none"),
             Ptk2      + theme(legend.position = "none"),
             Ythdf2    + theme(legend.position = "none"),
             Hnrnpa2b1 + theme(legend.position = "none"),
             rel_widths = c(1,1,0.3,1,1,1,1,1),
             nrow = 1)



save_plot_pdf("selected_genes_point_plots.pdf", selected_genes_point_plots, 
              width = unit(7.2, "in"), height = unit(2, "in") )

APC_related_genes_point_plots = plot_grid(
  Camsap1 + theme(legend.position = "none"),
  Cep120 + theme(legend.position = "none"),
  Numa1 + theme(legend.position = "none"),
  Ezh2 + theme(legend.position = "none"),
  nrow = 1
)

save_plot_pdf("APC_related_genes_point_plots.pdf", APC_related_genes_point_plots, 
              width = unit(5.2, "in"), height = unit(2, "in") )




Ncoa6 = 
  plot_cpm_across_conditions(all_counts, "Ncoa6", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 5 ) 
Ncoa6

Brd2 = 
  plot_cpm_across_conditions(all_counts, "Brd2", 
                             stages = c("1cell", "MII"), plot_type = "point", ymax = 3.2 ) 
Brd2

save_plot_pdf("Brd2.pdf", Brd2, 
              width = unit(3, "in"), height = unit(2.4, "in") )


Alkhb5 =   plot_cpm_across_conditions(all_counts, "Alkbh5", 
                                        stages = c("1cell", "2cell"), plot_type = "point", ymax = 3.2 ) 

Alkhb5

Alkhb5_bar =   plot_cpm_across_conditions(all_counts, "Alkbh5", 
                                      stages = c("1cell", "2cell"), plot_type = "barplot", ymax = 3.2 ) 

Alkhb5_bar

save_plot_pdf("Alkhb5.pdf", Alkhb5, 
              width = unit(3, "in"), height = unit(2.4, "in") )

################################################################################



plot_cpm_across_conditions_ribo = function(counts, gene, stages = c("all"), ymax = 4, plot_type = "point") { 
  # counts is a numeric matrix of raw read counts + gene_ids
  # gene to be plotted. Formatted as "Obox2"
  # Stages are 1cell, 2cell, MII, GV, etc
  set.seed(3)
  colnames(counts) = c( "transcript" ,  paste ( colnames(counts)[-1], 'Ribo', sep = "-" ) )
  
  # A monkey patch for the y-axis aesthetics
  ybreaks = waiver()
  
  if(ymax < 4){
    ybreaks = 1:ymax
  }
  
  if (stages[1] == "all" )  {  
    normalizedcounts = NormalizeData(counts[,-1], 
                                     scale = 10000, 
                                     normalization.method = "CLR", 
                                     margin = 2)
  } 
  else { 
    selected_samples = grep(paste(stages, collapse  = "|"), colnames(counts))
    normalizedcounts = NormalizeData(counts[,selected_samples], scale = 10000,
                                     normalization.method = "CLR", 
                                     margin = 2)
  }
  
  tidy_norm        = melt(normalizedcounts[ strip_extension(as.character(counts$transcript)) %in% gene,  ] )
  tidy_norm$stage  = sapply(strsplit(as.character(colnames(normalizedcounts)), split = "-"), "[[", 1) 
  tidy_norm$stage  = factor( tidy_norm$stage, levels = c("GV", "MII", "1cell",  "2cell",  "4cell", "8cell")  )
  tidy_norm$method = sapply(strsplit(as.character(colnames(normalizedcounts)), split = "-"), "[[", 3) 
  
  p =
    ggplot(tidy_norm, 
           aes( x =  strip_extension(as.character(gene)), y = value) ) + 
    # stat_summary(aes( shape = stage, color = method), fun.data = "mean_se", size = 1, 
    #              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1), 
    #              show.legend = FALSE, inherit.aes = FALSE) +
    geom_jitter(aes( shape = stage, color = method),
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1),
                size = 2 ) +
    # scale_color_manual(values =  c('#089099', '#d12959')) + 
    theme_bw()+
    theme(
      axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      legend.text       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
    ) + 
    scale_y_continuous( limits = c(0, ymax) , breaks = ybreaks ) + 
    scale_color_manual(values = c("Ribo" = ribo_orange, "RNAseq" = rna_blue)) + 
    scale_shape_manual(values = c(4, 16, 3, 6)) + 
    ylab("Normalized Count") + xlab("") 
  
  if(plot_type == "point"){
    return(p)
  }
  else if (plot_type == "barplot"){
    mean_se_tidy = tidy_norm %>% group_by(stage, method) %>%
      summarise(mean = mean(value), se = se(value))
    
    p =
      ggplot(mean_se_tidy, aes(x=stage , y=mean,  fill=method)) + 
      geom_bar(position=position_dodge(), stat="identity", width = 0.6) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                    width=.25, # Width of the error bars
                    position=position_dodge(.9)) + 
      theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
            panel.border     = element_blank(),
            panel.grid       = element_blank(),
            panel.background = element_blank(),
            axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            legend.title     = element_blank(),
            legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
            legend.key.size  = unit(0.15, 'in'),
            legend.position  = "none", 
            axis.line.y      = element_line(colour = "black", size = 0.35),
      ) +
      scale_y_continuous( limits = c(0, ymax) , expand = c(0, 0), breaks = ybreaks ) + 
      scale_fill_manual(values = c("Ribo" = ribo_orange, "RNAseq" = rna_blue)) +  
      #labs(title = gene) + 
      ylab("Normalized Count") + 
      xlab(gene)
    
    return(p)
  }
  else { 
    print("Incorrect plot type")}
  
}

plot_gv_versus_mii = function(counts, gene, ymax){
  return( plot_cpm_across_conditions_ribo(counts, gene, 
                                          stages = c("MII", "GV"), 
                                          plot_type = "barplot", ymax = ymax )  )
}



# For the main figure
Bub1b = plot_gv_versus_mii(rcw, "Bub1b", ymax = 6 ) 
Cdc20 = plot_gv_versus_mii(rcw, "Cdc20", ymax = 4 )
Dcp1a = plot_gv_versus_mii(rcw, "Dcp1a", ymax = 5 )

Cdc27 = plot_gv_versus_mii(rcw, "Cdc27", ymax = 4 ) 
Cdc40 = plot_gv_versus_mii(rcw, "Cdc40", ymax = 3 ) 
Lin7c = plot_gv_versus_mii(rcw, "Lin7c", ymax = 3 ) 
Mos   = plot_gv_versus_mii(rcw, "Mos", ymax = 3 ) 
Pcm1  = plot_gv_versus_mii(rcw, "Pcm1", ymax = 4 ) 
Pum2 = plot_gv_versus_mii(rcw, "Pum2", ymax = 3 ) 
Suz12 = plot_gv_versus_mii(rcw, "Suz12", ymax = 3 ) 

# High in GV vs MII
Cpeb1 = plot_gv_versus_mii(rcw, "Cpeb1", ymax = 4 ) 
Eif4a3 = plot_gv_versus_mii(rcw, "Eif4a3", ymax = 3 ) 
Rpl5 = plot_gv_versus_mii(rcw, "Rpl5", ymax = 5 )


Fzr1  = plot_gv_versus_mii(rcw, "Fzr1", ymax = 3 ) 
Rpl26 = plot_gv_versus_mii(rcw, "Rpl26", ymax = 3 ) 
Rpl36 = plot_gv_versus_mii(rcw, "Rpl36", ymax = 2 ) 
Rps17 = plot_gv_versus_mii(rcw, "Rps17", ymax = 3 ) 
Actb = plot_gv_versus_mii(rcw, "Actb", ymax = 3 ) 
Cd151 = plot_gv_versus_mii(rcw, "Cd151", ymax = 2 ) 
Gapdh = plot_gv_versus_mii(rcw, "Gapdh", ymax = 1 ) 
Mapk3 = plot_gv_versus_mii(rcw, "Mapk3", ymax = 3 ) 
Med7 = plot_gv_versus_mii(rcw, "Med7", ymax = 2 ) 
Miox = plot_gv_versus_mii(rcw, "Miox", ymax = 4 ) 
Paip2 = plot_gv_versus_mii(rcw, "Paip2", ymax = 4 ) 
Pfdn5 = plot_gv_versus_mii(rcw, "Pfdn5", ymax = 3 ) 
Qars = plot_gv_versus_mii(rcw, "Qars", ymax = 4 ) 
Rimkla = plot_gv_versus_mii(rcw, "Rimkla", ymax = 2 )
Stx5a  = plot_gv_versus_mii(rcw, "Stx5a", ymax = 3 )
Tinf2 = plot_gv_versus_mii(rcw, "Tinf2", ymax = 3 )

main_figure_gv_mii = plot_grid(Bub1b, Cdc20, Dcp1a , 
                               Cpeb1, Eif4a3, Rpl5,
                               nrow = 1)

save_plot_pdf("gv_vs_mii_main.pdf", main_figure_gv_mii, height = 1.25, width = 5)

supp_figure_gv_mii = plot_grid(Cdc27 ,
                               Cdc40 , 
                               Lin7c ,
                               Mos   ,
                               Pcm1  ,
                               Pum2 ,
                               Suz12 ,
                               Fzr1  ,
                               Rpl26 ,
                               Rpl36 ,
                               Rps17 ,
                               Actb ,
                               Cd151 ,
                               Gapdh ,
                               Mapk3 ,
                               Med7 ,
                               Miox ,
                               Paip2 ,
                               Pfdn5 ,
                               Qars ,
                               Rimkla ,
                               Stx5a  ,
                               Tinf2 ,
                               ncol = 8) 

supp_figure_gv_mii

save_plot_pdf("supp_figure_gv_mii.pdf", supp_figure_gv_mii, height = 3.75, width = 6.7)
