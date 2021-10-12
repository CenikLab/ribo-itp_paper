## Code to merge with Hakan
library(plotly)
library(ribor)
library(Seurat)
library(reshape2)
library(pheatmap)

mm = Ribo('../../mouse-itp_v5.ribo', rename = rename_default)
min_len = 29
max_len = 35

source('../figures/qc/Ribo_Summary_Function.R')

ribo_rc <- get_region_counts(mm,
                    
                             range.lower = min_len,
                             range.upper = max_len,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy = F,
                             alias       = TRUE,
                             region      = c("CDS"), 
                             compact = F)

rcw = dcast(ribo_rc, transcript ~ experiment)  
mouse_rna_cds = read.csv('../../cds_counts.csv.gz')
colnames(mouse_rna_cds)  = gsub(colnames(mouse_rna_cds), pattern = ".", fixed = T, replacement = "-")

for (id in 1:length(mouse_rna_cds$transcript)) { 
  mouse_rna_cds$transcript[id] = rename_default(mouse_rna_cds$transcript[id])
}
## Remove 20210607.RNAseq.4cell.cross.B
mouse_rna_cds = mouse_rna_cds[,-11]

# Mouse_RNA_CDS as before and doesn't include "X20210607-RNAseq-4cell-cross-B"
all_counts = merge(mouse_rna_cds, rcw, by = "transcript")

# Convert ribo and rna-seq into density 
cds_len = get_region_lengths(mm, alias = T)[,c(1,4)]
# Density corresponds to Reads per 1kb 

allcounts_density = all_counts
for (col_id in 2:48) { 
  allcounts_density[,col_id] = 1000* allcounts_density[,col_id] / 
    cds_len$CDS[match( allcounts_density$transcript, cds_len$transcript )]  
}

# CLR normalize density and average across replicates
set.seed(3)
allcounts_cds_clr_density = NormalizeData(allcounts_density[,-1], scale.factor = 10000, normalization.method= "CLR",verbose = TRUE)
rownames(allcounts_cds_clr_density) = all_counts[,1]

stage = unlist(lapply ( strsplit(colnames(all_counts)[-1], split = "-" ) , "[[", 3 )  ) 
method = unlist(lapply ( strsplit(colnames(all_counts)[-1], split = "-" ) , "[[", 2 )  ) 
groups = as.factor(paste(method, stage,  sep= "_") ) 

mean_clr_density = t(apply(allcounts_cds_clr_density, 1, FUN = function (x){ tapply (as.numeric(x), groups, mean)}) ) 
prots = read.csv ('../../Protein_MouseEmbryo_PMID_29281840.csv')

## Add protein expression
prots_rna_ribo = mean_clr_density[strip_extension(rownames(mean_clr_density) ) %in% prots$Gene.Name, ] 
prots_rna_ribo = cbind( prots_rna_ribo, 
                        prots[match(strip_extension(row.names(prots_rna_ribo)) , prots$Gene.Name), ]  ) 

## Remove annotations/blastocyst from the proteomics data
prots_rna_ribo = prots_rna_ribo[, -c(19:22)]

# Note that 451 ids in proteomics data did not match ids in ribo/rna
dim(prots_rna_ribo)

# Subset to non-zero expression; 29 transcripts with no detectable rna or ribo 
prots_rna_ribo = prots_rna_ribo[apply(prots_rna_ribo[,1:12], 1, sum) > 0, ]
dim(prots_rna_ribo)


## We should be using Spearman's correction before reporting these. 
## Estimation of the reliability of the assays based on replicate correlation: 
rna_reliability = 0.79
ribo_reliability = 0.71
prot_reliability = 0.8


## Note that there are some negative correlation that will get replaced with zero in the following figure
## If using this version might either change color scheme or label as less than zero
## Will split this into RNA - Prot and Ribo-Prot 
rna_prot = prots_rna_ribo[ , c(13, 11, 12, 7:10, 14:18)]
ribo_prot = prots_rna_ribo[, c(13, 5, 6, 1:4,14:18 )]

rna_prot_cors = cor(rna_prot[,-1], method = "spearman") / sqrt(prot_reliability * rna_reliability)
ribo_prot_cors = cor(ribo_prot[,-1], method = "spearman") / sqrt(prot_reliability * ribo_reliability)

rna_prot_cors = rna_prot_cors[1:6, 7:11]
ribo_prot_cors = ribo_prot_cors[1:6, 7:11]

rna_prot_cors
ribo_prot_cors



feature_color_palette = colorRampPalette(c("#3361A5", "#2884E7", "#1BA7FF", "#76CEFF", "#FFFFFF", "#FFE060", "#FA8E24", "#DA2828", "#A31D1D"), space="Lab")
number_cols= 50

color_palette_prot = feature_color_palette(number_cols)
breaks_manual = seq(-0.75, 0.75, length.out = number_cols)
values_to_remove = 2:25
breaks_manual = breaks_manual[-values_to_remove]
color_palette_prot = color_palette_prot[-(values_to_remove-1)]

rna_prot_heatmap = 
  pheatmap(rna_prot_cors, 
           #labels_row     = strip_extension(rcw_columns_reordered[variables$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           cluster_rows   = FALSE, 
           fontsize       = 8,
           fontsize_col   = 7,
           color        = color_palette_prot,
          breaks = breaks_manual, 
           fontsize_row   = 7)

ribo_prot_heatmap = 
  pheatmap(ribo_prot_cors, 
           #labels_row     = strip_extension(rcw_columns_reordered[variables$vst.variance.standardized > variance_threshold,1]), 
           cluster_cols   = FALSE, 
           cluster_rows   = FALSE, 
           fontsize       = 8,
           fontsize_col   = 7,
           color        = color_palette_prot,
           breaks = breaks_manual, 
           fontsize_row   = 7)

plot_grid(ribo_prot_heatmap[[4]], rna_prot_heatmap[[4]])

## UMAP clustering of the rank based data
library(M3C, include.only = "umap")
ranking = apply ( prots_rna_ribo_density[, -c(13,19)], 2, rank)
umap(ranking, dotsize = 2, axistextsize = 12, legendtextsize = 12, textlabelsize = 4, 
     legendtitle = "", text = c("", "", "4 cell", "8 cell", rep("",12), "Morula"),
     labels = as.factor(c(rep("Ribo", 6) , rep( "RNA", 6), rep ("Prot", 5) )), 
) 

plotMDS(ranking)
cluster::pam(t(ranking), k = 5)


## Sankey Diagram of the Correlations
exp_stages = c("GV", "MII", "1cell", "2cell", "4cell", "8cell")
prot_stages = c("1cell", "2cell", "4cell", "8cell", "morula")
prots_rna_ribo_density_ordered = prots_rna_ribo[,c(5,6, 1:4, 14:18, 11, 12, 7:10)]


# Remove 8-cell protein add blastocyst
ribo_orange = rgb(228,88,10 , maxColorValue = 255)
rna_blue   = rgb(55,135,192, maxColorValue = 255)

# All columns
selected_ribo = 1:6
selected_rna = 13:17
selected_protein = c(7:9, 11)

# Split into two graphs. One for GV - 1cell; One drawback is the ratio between two is warped. 
selected_ribo = 1:3
selected_rna = 13:15
selected_protein = c(7:9)

## 
# The other for 4,8 cell to morula
selected_ribo = 5:6
selected_rna = 16:17
selected_protein = 11

selected_cols = c(selected_ribo, selected_protein, selected_rna)

## We can create two link_lists and make two diagram to get better spacing
link_list = list( 
      source = c(), 
      target = c(), 
      value = c()
)
  # We will only calculate ITP-Prot; RNA-Prot correlation with Spearman's corrrection
  # We will update link list with >0 correlations
for (column_id in c(selected_ribo, selected_rna)  ) { 
  for (prot_id in selected_protein) { 
    if (column_id < 7) { 
        spearman = cor(prots_rna_ribo_density_ordered[,column_id], prots_rna_ribo_density_ordered[,prot_id], method = "spearman") / sqrt(prot_reliability * ribo_reliability)
      } else { 
        spearman = cor(prots_rna_ribo_density_ordered[,column_id], prots_rna_ribo_density_ordered[,prot_id], method = "spearman") / sqrt(prot_reliability * rna_reliability)
      }
    if (spearman >  0 ) {
      if (column_id < 7) { 
        link_list[["source"]] = c(link_list[["source"]], column_id-1)
        link_list[["target"]] = c(link_list[["target"]], prot_id-1)
        link_list[["value"]] = c(link_list[["value"]], round(spearman, 2) ) 
      } else { 
        link_list[["source"]] = c(link_list[["source"]], prot_id-1)
        link_list[["target"]] = c(link_list[["target"]], column_id-1)
        link_list[["value"]] = c(link_list[["value"]], round(spearman, 2) ) 
      }
    }
      
  }
}

cor_breaks = seq(0,1, by=.1)
cor_breaks = c(0,0.2, 0.4, 0.5, 0.6, 0.7, 0.8)
# Based on quantile of link_list value
cor_breaks = c(0,0.18, 0.32, 0.46, 0.54, 0.731)

# Based on hand selected 
cor_breaks = c(0,0.35, 0.42, 0.47, 0.52, 0.58, 0.731)


link_list$color = cut(link_list$value, breaks = cor_breaks, 
                    labels = colorRampPalette(
                      colors = c("lightgray", "darkblue"))(length(cor_breaks)-1 ) )


node_list_all = list(
  label = colnames(prots_rna_ribo_density_ordered),
  color = c(rep(ribo_orange, 6), rep("green", 5), rep(rna_blue, 6) ),
  y = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 
        0.3, 0.4, 0.6, 0.8,
        0.1, 0.2, 0.3, 0.4, 0.6, 0.7), 
  x = c(rep(0.1, 6), 
        rep(0.3, 4), 
        rep(0.6, 6)), 
  pad = 40,
  thickness = 20,
  line = list(
    color = "black",
    width = 0.5
  )
)


node_list_48 = list(
  label = colnames(prots_rna_ribo_density_ordered),
  color = c(rep(ribo_orange, 6), rep("green", 5), rep(rna_blue, 6) ),
  y = c(0.1, 0.5, 0.55, 0.1, 0.5), 
  x = c(0.1, 0.1, 
        0.3, 
        0.6, 0.6), 
  pad = 40,
  thickness = 16,
  line = list(
    color = "black",
    width = 0.5
  )
)

node_list_GV_2 = list(
  label = colnames(prots_rna_ribo_density_ordered),
  color = c(rep(ribo_orange, 6), rep("green", 5), rep(rna_blue, 6) ),
  y = c(0.1, 0.34, 0.55, 
        0.3, 0.56, 0.77, 
        0.05, 0.24, 0.45), 
  x = c(0.1, 0.1, 0.1, 
        0.3, 0.3, 0.3,
        0.6, 0.6, 0.6), 
  pad = 40,
  thickness = 16,
  line = list(
    color = "black",
    width = 0.5
  )
)


## Node color added manually
fig <- plot_ly(
  type = "sankey",
  orientation = "v",
  textfont = list( family = "Arial"),

#  node = node_list_all , 
  node = node_list_48 , 
#  node =node_list_GV_2,
  link = link_list
)

fig <- fig %>% layout(
  title = "Spearman Correlation",
  font = list(
    size = 10
  )
)

fig

## Orca requires the commandline-tools separately installed
orca(fig, "sankey.pdf")


