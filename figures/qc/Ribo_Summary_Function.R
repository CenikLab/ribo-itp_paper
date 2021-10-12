## This function takes a ribo file as input and creates the usual output graphs
## Determine min_len; max_len based on region counts
# source('~/Documents/jungla_source/cancer_jungla/mapCode-colors.r')

subsampleMatrix =  function (counts, desiredDepth) 
{
  rns <- rownames(counts)
  n = nrow(counts)
  proportion =  desiredDepth / colSums(counts)
  if (sum(proportion >1) ) { stop(paste("DesiredDepth must be less than library size for all!"))}
  ret = counts
  for (experiment in 1:ncol(counts)){ 
    ret[,experiment] = rbinom(n, counts[,experiment], proportion[experiment])  
  }
  rownames(ret) <- rns
  ret
}

summarize_ribo_results = function (ribo_file, min_len = 28, max_len = 35) { 
  library(ribor)
  library(edgeR)
  library(reshape2)
  library(pheatmap)
  p2 = plot_length_distribution(x           = ribo_file,
                           region      = "CDS",
                           range.lower = 15,
                           range.upper = 40,
                           fraction    = TRUE)
  print(p2)
  
  for(seq_len in 15:40) {
    p1 = plot_region_counts(x           = ribo_file,
                            range.lower = seq_len,
                            range.upper = seq_len, 
                            title= seq_len)
    print (p1)
  }
  
  
  p3 = plot_region_counts(x           = ribo_file,
                     range.lower = min_len,
                     range.upper = max_len)
  print(p3)
  
  p4 = plot_metagene(ribo_file,
                site        = "stop",
                normalize   = TRUE,
                title       = "Stop Site Coverage",
                range.lower = min_len,
                range.upper = max_len)
  print(p4)
  
  p5 = plot_metagene(ribo_file,
                site        = "start",
                normalize   = TRUE,
                title       = "Start Site Coverage",
                range.lower = min_len,
                range.upper = max_len)
  print(p5)
  
  ribo_rc <- get_region_counts(ribo_file,
                               range.lower = min_len,
                               range.upper = max_len,
                               length      = TRUE,
                               transcript  = FALSE,
                               tidy = F,
                               alias       = TRUE,
                               region      = c("CDS"), 
                               compact = F)
  
  rcw = dcast(ribo_rc, transcript ~ experiment)  
  expressed_ribo = rowSums ( cpm(rcw[,-1]) > 1) > 2
  c3 = cor(rcw[expressed_ribo,-1], method = "spearman")
  pheatmap(c3, main  = "Ribosome Profiling")
  
  rnaseq <- get_rnaseq(ribo.object = ribo_file,
                       tidy        = F,
                       compact = F, 
                       region = "CDS",
                       alias       = TRUE)
  rnaseq_w = dcast(rnaseq, transcript ~ experiment)
  
}

strip_extension = function(genenames) {
  return(sapply ( strsplit( genenames, split = "-") , "[[", 1 ) )
}


plot_pairwise_relationships <- function (counts_w, id1, id2, 
                                         xlab = "", ylab= "", num_bin = 52, 
                                         xrange = 100000, yrange = 100000  ) { 
    library(ggpubr)
    sp <- ggscatter(counts_w, x = id1, y = id2,
        #                add = "reg.line", conf.int = FALSE,     
        #                add.params = list(color = "blue", size = 0.5),
                  font.family = "Helvetica", size = 0.2,
                  color = "gray", alpha = 0.3, ggtheme = theme_bw()) 
  formatted =   sp +   
    scale_x_log10(labels = scales::label_number_si(), limits = c(0.3, xrange)) +   
    scale_y_log10(labels = scales::label_number_si(), limits = c(0.3, yrange)) + 
    labs (x=xlab, y = ylab) +
    stat_cor(method = "spearman", 
             aes(label = ..r.label..), 
             cor.coef.name = "rho", digits = 2)  + 
    geom_hex(bins= num_bin, aes(alpha=log10(..count..) ), fill="#bf5700" )   
  return (formatted)  
}

