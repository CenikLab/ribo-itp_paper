#SNP Figure

riboseq_overall_percentage_file = "../snp/snp_dataframes/riboseq_snp_percentages.csv"
rnaseq_overall_percentage_file  = "../snp/snp_dataframes/rnaseq_snp_percentages.csv"

overall_ribo_to_rna_ratios_file = "../snp/snp_dataframes/refined_ribo_to_rna_ratios.csv"
overall_ribo_to_rna_binary_file = "../snp/snp_dataframes/refined_ribo_to_rna_binary_significance.csv"



detailed_riboseq_count_file             = "../snp/snp_dataframes/riboseq_detailed_snps.csv.gz"
detailed_rnaseq_count_file              = "../snp/snp_dataframes/rnaseq_detailed_snps.csv.gz"

output_folder           = "./pdfs"

################################################################################
########                    L I B R A R I E S                          #########

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(cowplot)
library(RColorBrewer)

library(Cairo)

#Heatmap related packages
library(pheatmap)

# We don't need the following two libraries.
#library(seriation)
#library(dendextend)

################################################################################
#########                  C O L O R I N G                             #########

BURNT_ORANGE = "#bf5700"
UT_BLUE      = "#005f86"

PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#f0756e"
PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#7ccba2"
PERCENTAGE_BARPLOT_OTHER_COLOR    = "#045175"

#PERCENTAGE_BARPLOT_PATERNAL_COLOR =  "#999999"
#PERCENTAGE_BARPLOT_MATERNAL_COLOR = "#E69F00"
#PERCENTAGE_BARPLOT_OTHER_COLOR    = "#56B4E9"

PERCENTAGE_BARPLOT_COLORS = c(PERCENTAGE_BARPLOT_OTHER_COLOR,
                              PERCENTAGE_BARPLOT_MATERNAL_COLOR,
                              PERCENTAGE_BARPLOT_PATERNAL_COLOR
                              )




#MAIN_PERCENTAGE_RNASEQ_FILL_COLOR     = "skyblue"
#MAIN_PERCENTAGE_RIBOSEQ_FILL_COLOR    = "#abe39d"
#MAIN_PERCENTAGE_ERRORBAR_COLOR = "#413bb8"

FOURCELL_BACKGROUND_COLOR = "#d9d9d9"

COUNT_NORMALIZATION_FACTOR = 10000

#RATIO_RIBO_COLOUR = "orange"
#RATIO_RNA_COLOUR  = "blue"

RATIO_RIBO_COLOUR = BURNT_ORANGE
RATIO_RNA_COLOUR  = "navy"

PERCENTAGE_DASHED_COLOR = "#088f99"

################################################################################
#########                 F O N T   S I Z E S                          #########

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 10

PDF_resolution = 600
FIGURE_FONT = "helvetica"

################################################################################
##### Experiments



# rename_experiments = function(raw_names){
#   plain_names = unlist(lapply ( strsplit(  raw_names, split = "-" ), "[[", 3) )
#   
#   uniq_names = unique( plain_names )
#   
#   enumerator = c()
#   
#   for(n in uniq_names){
#     this_count = length( grep( n, plain_names )  )
#     enumerator = c(enumerator, 1:this_count)
#   }
#   
#   actual_names = paste( plain_names, enumerator, sep = "-" ) 
#   
#   return(actual_names)
# }




# rename_experiments_2 = function(raw_names){
#   plain_names = unlist(lapply ( strsplit(  raw_names, split = "-" ), "[[", 3) )
#   
#   uniq_names = unique( plain_names )
#   
#   enumerator = c()
#   
#   for(n in uniq_names){
#     this_count = length( grep( n, uniq_names )  )
#     enumerator = c(enumerator, 1:this_count)
#   }
#   
#   actual_names = paste( plain_names, enumerator, sep = "-" ) 
#   
#   return(actual_names)
# }


make_name_mapper = function(experiment_names){
  exp_names        = sort( experiment_names )
  unique_raw_names = unique(exp_names )
  group_names      = sort(unique(unlist( lapply(strsplit( unique_raw_names, split = "-" ) , "[[", 3 ) ) ) )
  
  original_exp_names     = c()
  new_experiment_names   = c()
  list_counter           = 1
  experiment_name_mapper = list()
  
  for(g in group_names){
    elements_of_group = unique_raw_names[grep(  g, unique_raw_names )]
    size_of_group     = length( elements_of_group  )
    
    for( i in seq(1:size_of_group) ){
      this_basename = strsplit( elements_of_group[i], split = "-"   )[[1]][3]
      new_name = paste(  this_basename, i , sep = "-" )
      print( c( elements_of_group[i], new_name )  )
      
      original_exp_names   = c(original_exp_names, elements_of_group[i])
      new_experiment_names = c(new_experiment_names, new_name) 
      
      experiment_name_mapper[[list_counter]] = new_name
      list_counter = list_counter + 1
    }
    
  }
  
  names(  experiment_name_mapper ) = original_exp_names
  
  return(experiment_name_mapper)
}

temp_riboseq_df = read.csv(riboseq_overall_percentage_file, row.names = 1)
temp_rnaseq_df  = read.csv(rnaseq_overall_percentage_file, row.names = 1)


ribo_exp_name_mapper = make_name_mapper( row.names(temp_riboseq_df)  )
rna_exp_name_mapper  = make_name_mapper( row.names(temp_rnaseq_df) )

experiment_name_mapper = c(ribo_exp_name_mapper, rna_exp_name_mapper)

rename_experiments = function(experiment_name){
  return(experiment_name_mapper[[experiment_name]])
}

rename_experiments = Vectorize(rename_experiments, USE.NAMES = FALSE)




################################################################################

rename_genes = function(gene_names, split = "."){
  new_names = unlist(lapply ( strsplit(  gene_names, split = split, fixed = TRUE ), "[[", 1) )
  return(new_names)
}


# rename_genes = function(gene_names){
#   new_names = unlist(lapply ( strsplit(  gene_names, split = ".", fixed = TRUE ), "[[", 1) )
#   return(new_names)
# }


################################################################################

##### PERCENTAGES SUPPLEMENTARY

riboseq_overall_percentages_pre = read.csv( riboseq_overall_percentage_file, 
                                        row.names = 1 )

rnaseq_overall_percentages_pre  = read.csv( rnaseq_overall_percentage_file, 
                                        row.names = 1 )

riboseq_raw_experiment_names               = row.names(riboseq_overall_percentages_pre)
riboseq_overall_percentages_pre$Experiment = rename_experiments(riboseq_raw_experiment_names)

rnaseq_raw_experiment_names               = row.names(rnaseq_overall_percentages_pre)
rnaseq_overall_percentages_pre$Experiment = rename_experiments(rnaseq_raw_experiment_names)

# Prep the Data for ggplot
riboseq_overall_percentages = melt(riboseq_overall_percentages_pre, value.name = "Percentage")
rnaseq_overall_percentages  = melt(rnaseq_overall_percentages_pre, value.name = "Percentage")

# Reorder the columns so that 1-cell appears on top
riboseq_overall_percentages$Experiment = 
    factor(riboseq_overall_percentages$Experiment, 
           levels = unique(rev( riboseq_overall_percentages$Experiment )) )

rnaseq_overall_percentages$Experiment = 
  factor(rnaseq_overall_percentages$Experiment, 
         levels = unique(rev( rnaseq_overall_percentages$Experiment )) )

################################################################################
## SUPPLEMENTARY FIGURE FOR DETAILED ALLELE PERCENTAGES ###############################

# Ribosome Profiling

supplementary_allele_percentages_ribo = 
ggplot(data=riboseq_overall_percentages, aes(x=Experiment, y=Percentage, fill = factor(variable, levels = c("Other", "Maternal", "Paternal")  ) ) )  +
  geom_bar(stat="identity" ) +
  coord_flip() + 
  guides(fill=guide_legend(title="Type")) + 
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
        )  + 
  scale_fill_manual(values = PERCENTAGE_BARPLOT_COLORS ) + 
  geom_hline(yintercept = 50, linetype="dotted", 
             color = PERCENTAGE_DASHED_COLOR, size=0.6) + 
  labs(title = "Allele Percentages of Ribosome Footprints")


supplementary_allele_percentages_ribo

# RNA-Seq

supplementary_allele_percentages_rnaseq = 
ggplot(data=rnaseq_overall_percentages, aes(x=Experiment, y=Percentage, fill = factor(variable, levels = c("Other", "Maternal", "Paternal")  ) ) )  +
  geom_bar(stat="identity" ) +
  coord_flip() + 
  guides(fill=guide_legend(title="Type")) + 
  theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(family = "helveticassss", face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title     = element_blank(),
        legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE))  + 
  scale_fill_manual(values = PERCENTAGE_BARPLOT_COLORS ) + 
  geom_hline(yintercept = 50, linetype="dotted", 
             color = PERCENTAGE_DASHED_COLOR, size=0.6) + 
  labs(title = "Allele Percentages")


supplementary_allele_percentages_rnaseq

################################################################################
#####       M A I N    F I G U R E   P E R C E N T A G E S      ################

get_percentage_averages = function( input_df, label ){
  result = 
    input_df %>%
    filter(variable == "Paternal") %>%
    select(Experiment, variable, Percentage) %>%
    group_by( stage = unlist(lapply ( strsplit(  as.vector(Experiment), split = "-" ), "[[", 1) )  ) %>%
    mutate( average_percentage = mean( (Percentage)  )  ) %>%
    filter(stage != "GV") %>%
    select(stage, average_percentage, Percentage) %>%
    summarise(average_percentage = mean( (Percentage)), sd_percentage = sd(Percentage) ) %>%
    mutate(experiment_type = label)
  
  return(result)
}

rnaseq_percentage_averages  = get_percentage_averages(rnaseq_overall_percentages, "rna")
riboseq_percentage_averages = get_percentage_averages(riboseq_overall_percentages, "ribo")    

percentage_averages = bind_rows( riboseq_percentage_averages, rnaseq_percentage_averages  )

percentage_averages$stage = factor( percentage_averages$stage, levels = c("MII", "1cell", "2cell", "4cell", "8cell") )

### We won't use this plot!!!!
### See the corrected version below
ggplot(data=percentage_averages, aes(x=stage, y=average_percentage, fill = experiment_type )  )  + 
  geom_bar(position = "dodge", stat="identity", alpha = 0.7 ) +
  geom_errorbar( aes(  x    = stage, 
                       ymin = average_percentage - sd_percentage, 
                       ymax = average_percentage + sd_percentage), 
                 width    = 0.4, 
                 colour   = MAIN_PERCENTAGE_ERRORBAR_COLOR, 
                 alpha    = 0.9, 
                 size     = 1.3,
                 position = position_dodge(width = 0.8)) + 
  
  theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title     = element_blank(),
        legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  
  scale_fill_manual( values = c(MAIN_PERCENTAGE_RIBOSEQ_FILL_COLOR,
                                MAIN_PERCENTAGE_RNASEQ_FILL_COLOR) ) + 
  
  guides(fill=guide_legend(title="")) + 
  
  labs(title = "Percentage of paternal alleles", 
       x     = "Stage", 
       y     = "Percentage")

################################################################################




get_stage_percentages = function( group, type ){
  if(type == "ribo"){
    this_df = riboseq_overall_percentages
  }
  else{
    this_df = rnaseq_overall_percentages
  }
  
  raw_percentages = this_df %>%
    filter(variable == "Paternal") %>%
    select(Experiment, variable, Percentage) %>%
    group_by( stage = unlist(lapply ( strsplit(  as.vector(Experiment), split = "-" ), "[[", 1) )  ) %>%
    filter(stage == group)
  
  return(raw_percentages$Percentage)
}








p_observed_to_p_real = function(p_observed, error){
  #All values must be percentages
  
  p_real = ( (300 * p_observed) - (error * 100) ) / (300 - 4*error)
  
  return(p_real)
}


## The actual error is 3 times the paternal ratio
## Because there are 2 other nucleotides 
## that the maternal SNP nucleotide can accidentally go to
ribo_percentage_error = mean(get_stage_percentages( group = "MII", type = "ribo" ) ) * 3
rna_percentage_error  = mean(get_stage_percentages( group = "MII", type = "rna" ) ) *3

corrected_ribo_percentages = 
  riboseq_overall_percentages %>% filter(variable == "Paternal") %>%
  mutate( Percentage = p_observed_to_p_real(Percentage, ribo_percentage_error) )

corrected_rna_percentages = 
  rnaseq_overall_percentages %>% filter(variable == "Paternal") %>%
  mutate( Percentage = p_observed_to_p_real(Percentage, rna_percentage_error) )


corrected_rnaseq_percentage_averages  = get_percentage_averages(corrected_rna_percentages, "rna")
corrected_riboseq_percentage_averages = get_percentage_averages(corrected_ribo_percentages, "ribo")    

corrected_percentage_averages = bind_rows( corrected_riboseq_percentage_averages, corrected_rnaseq_percentage_averages  )

corrected_percentage_averages$stage = factor( corrected_percentage_averages$stage, levels = c("MII", "1cell", "2cell", "4cell", "8cell") )


main_paternal_percentage_figure = 
ggplot(data=corrected_percentage_averages, aes(x=stage, y=average_percentage, fill = experiment_type )  )  + 
  geom_bar(position = "dodge", stat="identity", alpha = 0.6 ) +
  geom_errorbar( aes(  x    = stage, 
                       ymin = average_percentage - sd_percentage, 
                       ymax = average_percentage + sd_percentage,
                       color = (experiment_type)), 
                 width    = 0.4, 
                 alpha    = 1, 
                 size     = 1.3,
                 position = position_dodge(width = 0.8)) + 
  
  theme(plot.title       = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border     = element_blank(),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.title     = element_blank(),
        legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)) + 
  
  scale_fill_manual( name = "experiment_type", values = c(RATIO_RIBO_COLOUR,
                                RATIO_RNA_COLOUR) ) + 
  scale_color_manual(name = "experiment_type",  values = c(RATIO_RIBO_COLOUR,
                                RATIO_RNA_COLOUR) ) + 
  labs(title = "Percentages of paternal alleles", 
       x     = "Stage", 
       y     = "Percentage")

main_paternal_percentage_figure


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#####################      H E A T M A P      ##################################

heatmap_df_raw             = read.csv(overall_ribo_to_rna_ratios_file, row.names = 1)
colnames( heatmap_df_raw ) = rename_genes(colnames( heatmap_df_raw ) )
heatmap_df                 = t(heatmap_df_raw[-1,])

heatmap_binary_raw             = read.csv(overall_ribo_to_rna_binary_file, row.names = 1)
colnames( heatmap_binary_raw ) = rename_genes(colnames( heatmap_binary_raw ) )
heatmap_binary                 = t(heatmap_binary_raw[-1,])
  
# Adapted from https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
paletteLength = 100

#myColor       = colorRampPalette(c("navy", "white", "orange"))(paletteLength)
myColor       = colorRampPalette(c(RATIO_RNA_COLOUR, "white", RATIO_RIBO_COLOUR))(paletteLength)
#myColor       = colorRampPalette(c("navy", "white", "#f8971f"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(heatmap_df_raw, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1), 
#              seq(max(heatmap_df_raw, na.rm = TRUE)/paletteLength, max(heatmap_df_raw, na.rm = TRUE), length.out=floor(paletteLength/2)))

myBreaks <- c(seq( -1 ,              0, length.out = ceiling(paletteLength/2) + 1), 
              seq( 1 /paletteLength, 1, length.out = floor(paletteLength/2)))


### Trimmed down
### For the main figure
main_heatmap_figure = 
pheatmap( heatmap_df , 
        # clustering_method = "median",
        #clustering_distance_rows = "correlation",
        #clustering_method = "centroid",
        #clustering_method = "complete",
        clustering_method = "ward.D",
        show_rownames     = FALSE,
        cutree_rows       = 4,

        cellwidth         = 40,
        treeheight_row    = 0,
        treeheight_col    = 0,
        cluster_cols      = FALSE, 
        color             = myColor,
        breaks            = myBreaks,
        na_col            = "white",
        angle_col         = 0,
        labels_col        = c("1cell", "2cell", "4cell", "8cell"),
        main              = "Paternal Ratio Difference",
        fontsize          = 10,
        fontsize_col      = FONT_LABEL_SIZE,
        fontsize_row      = 8,
        fontsize_number   = FONT_LABEL_SIZE
        )

main_heatmap_figure

#### Detailed with row names
#### For the supplementary Figure
supplementary_heatmap_figure = 
pheatmap( heatmap_df , 
          #clustering_method = "centroid",
          #clustering_method = "median",
          #clustering_method = "complete",
          clustering_method = "ward.D",
          legend_breaks     = c(-1,  -0.5, 0, 0.5, 1),
          legend            = T,
          show_rownames     = TRUE,
          cellwidth         = 40,
          cutree_rows       = 4,
          treeheight_row    = 0,
          cluster_cols      = FALSE, 
          color             = myColor,
          breaks            = myBreaks,
          na_col            = "white",
          angle_col         = 0,
          labels_col        = c("1cell", "2cell", "4cell", "8cell"),
          main              = "Paternal Ratio Difference",
          display_numbers   = heatmap_binary,
          fontsize          = 10,
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 8,
          fontsize_number   = FONT_LABEL_SIZE
          )

supplementary_heatmap_figure

################################################################################
################################################################################
################################################################################
################################################################################
########        G E N E  S P E C I F I C    P L O T S                   ########

detailed_riboseq_table_raw = read.csv(detailed_riboseq_count_file)
detailed_rnaseq_table_raw  = read.csv(detailed_rnaseq_count_file)

detailed_riboseq_table_raw$experiment = rename_experiments( detailed_riboseq_table_raw$experiment  )
detailed_riboseq_table_raw$transcript = rename_genes( detailed_riboseq_table_raw$transcript , split = "-" )

detailed_rnaseq_table_raw$experiment = rename_experiments( detailed_rnaseq_table_raw$experiment  )
detailed_rnaseq_table_raw$transcript = rename_genes( detailed_rnaseq_table_raw$transcript , split = "-" )

detailed_riboseq_table_raw$group = unlist(  lapply( strsplit(detailed_riboseq_table_raw$experiment, "-")  , "[[", 1 )  )
detailed_rnaseq_table_raw$group  = unlist(  lapply( strsplit(detailed_rnaseq_table_raw$experiment, "-")  , "[[", 1 )  )

experiment_groups = unique( detailed_riboseq_table_raw$group )
  
## Add counts per thousand to the dataframes
detailed_riboseq_table = 
  detailed_riboseq_table_raw %>% 
  group_by(experiment) %>%
  mutate( experiment_total = sum(paternal + maternal) ) %>%
  mutate( paternal_per_k = (paternal / experiment_total)* COUNT_NORMALIZATION_FACTOR,
          maternal_per_k = (maternal / experiment_total)*COUNT_NORMALIZATION_FACTOR ) %>%
  mutate(type = "ribo")


detailed_rnaseq_table = 
  detailed_rnaseq_table_raw %>% 
  group_by(experiment) %>%
  mutate( experiment_total = sum(paternal + maternal) ) %>%
  mutate( paternal_per_k = (paternal / experiment_total)* COUNT_NORMALIZATION_FACTOR,
          maternal_per_k = (maternal / experiment_total)*COUNT_NORMALIZATION_FACTOR ) %>%
  mutate(type = "rna")


detailed_table = rbind(detailed_riboseq_table, detailed_rnaseq_table)

###
# To handle graphs easier, we add dummy experiments to 2,4,8 cell groups
# so that they have equal (4) number of experiments.
# The new experiments have 0 counts.
detailed_riboseq_table_with_dummies =
  detailed_riboseq_table %>%
  filter(group == "2cell" | group == "4cell" | group == "8cell")

tmp_dummy_entry = detailed_riboseq_table %>% 
                    filter( experiment == "2cell-1" ) %>%
                    mutate( experiment = "2cell-4" ) %>%
                    mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_riboseq_table_with_dummies = rbind(detailed_riboseq_table_with_dummies, tmp_dummy_entry)


tmp_dummy_entry = detailed_riboseq_table %>% 
  filter( experiment == "4cell-1" ) %>%
  mutate( experiment = "4cell-4" ) %>%
  mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_riboseq_table_with_dummies = rbind(detailed_riboseq_table_with_dummies, tmp_dummy_entry)



detailed_rnaseq_table_with_dummies =
  detailed_rnaseq_table %>%
  filter(group == "2cell" | group == "4cell" | group == "8cell")

tmp_dummy_entry = detailed_rnaseq_table %>%
  filter(experiment == "4cell-1") %>%
  mutate(experiment = "4cell-3") %>%
  mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_rnaseq_table_with_dummies = rbind( detailed_rnaseq_table_with_dummies, tmp_dummy_entry )

tmp_dummy_entry = detailed_rnaseq_table %>%
  filter(experiment == "4cell-1") %>%
  mutate(experiment = "4cell-4") %>%
  mutate( paternal = 0, maternal = 0 , paternal_per_k = 0, maternal_per_k = 0 )

detailed_rnaseq_table_with_dummies = rbind( detailed_rnaseq_table_with_dummies, tmp_dummy_entry )

detailed_table_with_dummies = rbind(  detailed_riboseq_table_with_dummies, detailed_rnaseq_table_with_dummies )
##############################################################

get_gene_percentages = function(this_df, this_gene){
  ## Add percentage column to the table
  
  gene_data = this_df %>% 
                          filter(transcript == this_gene) %>%
                          group_by(  experiment, position ) %>%
                          mutate( maternal_total = sum(maternal), paternal_total = sum(paternal) ) %>%
                          mutate( paternal_percentage = (paternal / (paternal_total + maternal_total )) * 100  )

  return(gene_data)
}

##############################################################

##############################################################
##############################################################
##############################################################



################################################################################
################################################################################
################################################################################




##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

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

format_y_axis = function(x){
  return( sprintf("%.0f", x) )
}

make_y_axis = function(ymax){
  this_multiplier    = floor(ymax / 5)
  
  if( (ymax %% 5) >= 1 ) {
    this_multiplier = this_multiplier + 1
  } 
  
  yticks = seq(0, 5*this_multiplier, 5)
  
  return(yticks)
  
}

adjust_ymax = function( y ){
  this_multiplier    = floor(y / 5)
  
  if(y == 0){
    return(5)
  }
  
  if( (y %% 5) > 0 ) {
    this_multiplier = this_multiplier + 1
  } 
  
  return( 5*this_multiplier)
}

################################################################################


generate_legend_unit = function(gene, experiment_type){
  
  gene_data = 
    detailed_table_with_dummies %>%
    filter( transcript == gene & type == experiment_type ) %>%
    filter( group == "2cell" | group == "4cell" | group == "8cell")
  
  number_of_snps = length( unique(gene_data$position  ) )
  
  bar_palettes = list( ribo = "Oranges", rna = "Blues"  )
  title_colors = list( ribo = "orange", rna = "blue"  )
  title_texts  = list( ribo = "ribo", rna = "rna"  )
  
  barplot_colors = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps + 5 )[-seq(1,5)] 
  
  df = data.frame( x = rep(1,number_of_snps), y = rep(1,number_of_snps), position = seq(1,number_of_snps))
  
  p = ggplot(data=df, 
             aes(x    = x, 
                 y    = y, 
                 fill = factor(position, levels = sort(unique( df$position  ) )  )  )  )  + 
    geom_bar(stat="identity", width= 0.8, color = "black") + 
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.title.x     = element_blank(),
      #axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = 10),
      legend.position  = "none",
      axis.title.y     = element_blank(),
      axis.text.x = element_blank(),
      #axis.text.x = element_text(family = FIGURE_FONT, face = "plain", size = 10),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      plot.background = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(color = title_colors[[experiment_type]], size=10, face="bold", hjust = 0.5),
    ) + 
    scale_fill_manual(values = barplot_colors)  + 
    ggtitle(title_texts[[experiment_type]])
  
  #if(experiment_type == "rna"){
  #  p = p + scale_y_reverse()
  #}
  
  return(p)
}

generate_legend = function(gene){
  p_ribo = generate_legend_unit(gene, experiment_type = "ribo")
  p_rna  = generate_legend_unit(gene, experiment_type = "rna")
  
  this_blank_pad  = generate_blank_plot()
  this_legend_pre = plot_grid( this_blank_pad, p_ribo, this_blank_pad, p_rna, 
                           nrow = 1,
                           rel_widths = c(0.2, 1 , 0.2, 1 ) )
  this_legend     = plot_grid(this_blank_pad, this_legend_pre, this_blank_pad,
                              ncol = 1, rel_heights = c(1, 0.6, 1)) 
  
  return(this_legend)
}

snp_detailed_legend = generate_legend("Nin")

snp_detailed_legend
generate_legend("Nin")

################################################################################

this_gene = "Ncoa3"

a = 
detailed_table_with_dummies %>%
  filter(group == "2cell" | group == "4cell" | group == "8cell") %>%
  filter(transcript == this_gene) %>% 
  group_by(position) %>%
  filter( sum(paternal_per_k + maternal_per_k) > 0  )



################################################################################
### U N I T     P L O T 

unit_plot = function(gene, experiment_type, experiment_group, allele_type, ymax){
  
  bar_palettes = list( ribo = "Oranges", rna = "Blues"  )
  
  y_scale_expand = c(0,0)
  
  ymax = adjust_ymax(ymax)
  
  gene_data = 
    detailed_table_with_dummies %>%
    filter(group == "2cell" | group == "4cell" | group == "8cell") %>%
    filter(transcript == gene) %>% 
    group_by(position) %>%
    filter( sum(paternal_per_k + maternal_per_k) > 0  ) %>%
    filter( type == experiment_type ) %>%
    filter( group == experiment_group)
  
  number_of_snps = length( unique(gene_data$position  ) )
  
  paternal_color_palette = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps + 5 )[-seq( 1,  5) ]
  #paternal_color_palette = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps  ) 
  barplot_colors         = paternal_color_palette
  
  yticks = make_y_axis(ymax)
  
  if(experiment_group == "4cell"){
    background_canvas = element_rect(fill = FOURCELL_BACKGROUND_COLOR)
  }
  else{
    background_canvas = element_rect(fill = "white")
  }
  
  
  if(allele_type == "paternal"){
    plot_base = ggplot(data=gene_data, 
                       aes(x    = experiment, 
                           y    = paternal_per_k, 
                           fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  ) 
  }else{
    plot_base = ggplot(data=gene_data, 
                       aes(x    = experiment, 
                           y    = maternal_per_k, 
                           fill = factor(position, levels = sort(unique( gene_data$position  ) )  )  )  )
  }
  
  this_plot = 
    plot_base + 
    geom_bar(stat="identity", width= 1) + 
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.title.x     = element_blank(),
      legend.position  = "none",
      # plot.margin      = margin(5.5, 5.5, 0, 5.5),
      axis.title.y     = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      plot.background = background_canvas,
      axis.text.y = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      # axis.ticks.y=element_blank()
    ) +
    scale_fill_manual(values = barplot_colors) 
  
  if(experiment_type == "ribo"){
    if(allele_type == "maternal"){
      this_plot = this_plot + scale_y_reverse(position = "right", 
                                              labels   = format_y_axis, 
                                              breaks   = yticks,
                                              limits   = c(ymax, 0) , 
                                              expand   = y_scale_expand ) 
    }
    else{
      this_plot = this_plot + scale_y_continuous(position = "right", 
                                                 labels   = format_y_axis, 
                                                 breaks   = yticks,
                                                 limits   = c(0, ymax) , 
                                                 expand   = y_scale_expand )
    }
  } 
  else{
    if(allele_type == "maternal"){
      this_plot = this_plot + scale_y_reverse(labels = format_y_axis, 
                                              limits = c(ymax, 0), 
                                              breaks   = yticks,
                                              expand = y_scale_expand )
    }
    else{
      this_plot = this_plot + scale_y_continuous( labels = format_y_axis, 
                                                  limits = c(0, ymax),
                                                  breaks   = yticks,
                                                  expand = y_scale_expand  )
    }
  }
  
  return(this_plot)
}


paternal_maternal_unit_plot = function( gene, experiment_group, experiment_type ){
  
  gene_sums = 
    detailed_table_with_dummies %>%
    filter( transcript == gene & type == experiment_type ) %>%
    filter( group == experiment_group) %>%
    group_by(experiment) %>%
    mutate(paternal_sum = sum(paternal_per_k), maternal_sum = sum(maternal_per_k))

  ymax = max( max(gene_sums$paternal_sum), max( gene_sums$maternal_sum ) )
  
  paternal_plot= unit_plot(gene = gene, experiment_type = experiment_type, experiment_group = experiment_group, allele_type = "paternal", ymax = ymax)
  maternal_plot= unit_plot(gene = gene, experiment_type = experiment_type, experiment_group = experiment_group, allele_type = "maternal", ymax = ymax)
  
  this_plot = plot_grid(paternal_plot, maternal_plot, ncol = 1)
  
  return(this_plot)
}

paternal_maternal_unit_plot(gene = "Nin", experiment_group = "8cell", experiment_type = "rna")


plot_experiment_group = function(gene, experiment_group){
  ribo_plot = paternal_maternal_unit_plot(gene = gene, experiment_group = experiment_group, experiment_type = "ribo")
  rna_plot  = paternal_maternal_unit_plot(gene = gene, experiment_group = experiment_group, experiment_type = "rna")
  
  this_plot = plot_grid(  rna_plot, ribo_plot, ncol = 2 , align = "h", axis = "b" )
  
  return(this_plot)
}

plot_experiment_group(gene = "Nin", experiment_group = "8cell")


################################################################################

plot_snp_gene_detailed = function(gene){
# This is the main function that plots the snps in detail  

  
  plot_2cell_raw = plot_experiment_group(gene = gene, experiment_group = "2cell")
  
  title_2cell = ggdraw() + 
    draw_label(
      "2cell",
      fontfamily = 'helvetica',
      size = FONT_LABEL_SIZE
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  plot_2cell = plot_grid(
    title_2cell, plot_2cell_raw,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 2)
  )
  
  
  plot_4cell_raw =  plot_experiment_group(gene = gene, experiment_group = "4cell")  

  title_4cell = ggdraw() + 
    draw_label(
      "4cell",
      fontfamily = 'helvetica',
      size = FONT_LABEL_SIZE
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  plot_4cell = plot_grid(
    title_4cell, plot_4cell_raw,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 2)
  )  
  
  plot_8cell_raw = plot_experiment_group(gene = gene, experiment_group = "8cell")

  title_8cell = ggdraw() + 
    draw_label(
      "8cell",
      fontfamily = 'helvetica',
      size = FONT_LABEL_SIZE
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  plot_8cell = plot_grid(
    title_8cell, plot_8cell_raw,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 2)
  )  
  
  
  paternal_label = ggdraw() + 
    draw_label(
      "paternal count",
      fontface   = 'plain',
      fontfamily = 'helvetica',
      size       = FONT_LABEL_SIZE,
      angle      = 90
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  
  maternal_label = ggdraw() + 
    draw_label(
      "maternal count",
      fontface   = 'plain',
      fontfamily = 'helvetica',
      size       = FONT_LABEL_SIZE,
      angle      = 90
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  
  x_label = plot_grid(paternal_label, maternal_label, ncol = 1)
    
  separator_1         = generate_blank_plot()
  
  snp_detailed_legend = generate_legend(gene)
  
  raw_plot = plot_grid( x_label,    plot_2cell,  separator_1, 
                        plot_4cell, separator_1, 
                        plot_8cell, separator_1, snp_detailed_legend,
                        rel_widths  = c(0.2, 1, 0.2, 1, 0.2, 1, 0.1, 0.4)  , 
                        rel_heights =  c(0.8, 1, 1, 1, 1, 1, 1, 0.01),
                        ncol        = 8  )
  
  title_main = ggdraw() + 
    draw_label(
      gene,
      fontface   = 'plain',
      size       = FONT_TITLE_SIZE ,
      fontfamily = 'helvetica',
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  this_plot = plot_grid(
    title_main, raw_plot,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
  )
  
  return(this_plot)
}

plot_snp_gene_detailed(gene = "Ncoa3")  


plot_snp_gene_detailed(gene = "Pla2g4c")  

plot_snp_gene_detailed(gene = "Fbxo15")

plot_snp_gene_detailed(gene = "Top1")

plot_snp_gene_detailed(gene = "Top2a")


plot_snp_gene_detailed(gene = "Cacna1b")

################################################################################
################################################################################
################################################################################
################################################################################

compute_sd_of_ratios = function(paternal_counts, total_counts){
  paternal_mean = mean(paternal_counts)
  total_mean    = mean(total_counts)
  
  f = (paternal_mean / total_mean)
  
  this_cov = cov( paternal_counts, total_counts )
  
  paternal_sd = sd(paternal_counts)
  total_sd    = sd(total_counts)
  
  this_sd =  f * sqrt(  (paternal_sd / paternal_mean)**2 + 
                         (total_sd/ total_mean)**2 - 
                         2*(this_cov/ (paternal_mean* total_mean) )  )
  
  return(this_sd)
}


determine_propogated_sd = function(df){
  pre_data = df %>%  group_by(group, experiment) %>%
    summarise(experiment_total = sum(paternal_per_k + maternal_per_k) , experiment_paternal_total = sum(paternal_per_k))
  
  prop_sd = pre_data %>% group_by(group) %>%
    summarise( group_sd =  compute_sd_of_ratios( experiment_paternal_total, experiment_total )  ) %>%
    arrange(group)
  
  return(prop_sd$group_sd)
}



################################################################################

find_paternal_ratios = function(df){
 
  ratios = 
    df %>%
    group_by( group ) %>%
    summarise( paternal_total = sum(paternal), 
               maternal_total = sum(maternal), 
               overall_total  = paternal_total + maternal_total, 
               paternal_ratio = paternal_total / overall_total ) 
  
  return( ratios$paternal_ratio )
}




find_paternal_per_k_ratios = function(df){
  
  ratios = 
    df %>%
    group_by( group ) %>%
    summarise( paternal_total = sum(paternal_per_k), 
               maternal_total = sum(maternal_per_k), 
               overall_total  = paternal_total + maternal_total, 
               paternal_ratio = paternal_total / overall_total ) %>%
    arrange(group)
  
  return( ratios$paternal_ratio )
}

################################################################################

plot_snp_ratios = function(gene, ymax = 0){
  gene_data_main = 
    detailed_table %>% 
    filter( transcript == gene ) %>%
    filter ( group == "2cell" | group == "4cell" | group == "8cell"  )
  
  gene_data_ribo = 
    gene_data_main %>%
    filter( type == "ribo")
  
  gene_data_rna = 
    gene_data_main %>%
    filter( type == "rna")
  
  ribo_ratios = find_paternal_per_k_ratios( gene_data_ribo )
  rna_ratios  = find_paternal_per_k_ratios( gene_data_rna )
  
  ribo_sd = determine_propogated_sd( gene_data_ribo  )
  rna_sd  = determine_propogated_sd(  gene_data_rna )
  
  plot_df = data.frame(  paternal_ratios = c(ribo_ratios, rna_ratios), 
                         stage           = rep( c("2cell", "4cell", "8cell") , 2), 
                         sd              = c(ribo_sd, rna_sd), 
                         type =  c( rep("ribo", 3), rep("rna", 3)  )  )
  
  p = ggplot(data=plot_df, 
             aes(x    = stage, y = paternal_ratios, group = type )  )  +
    geom_point( aes(colour = type, shape = type), position = position_dodge(width = 0.1), size = 4 ) +
    geom_line(aes(linetype=type, colour = type), position = position_dodge(width = 0.1) ) + 
    scale_linetype_manual(values=c("solid", "dashed")) + 
    geom_errorbar(aes(ymin=paternal_ratios-sd, ymax=paternal_ratios+sd, color = type), position = position_dodge(width = 0.1), width = 0.1 ) + 
    scale_colour_manual(values = c(RATIO_RIBO_COLOUR, RATIO_RNA_COLOUR)   ) + 
    theme(
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
      panel.background  = element_blank(),
      axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      legend.title      = element_blank(),
      legend.text        = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)
    ) +
    labs(title = gene, y = "paternal ratio") 
  
  if(ymax > 0){
    p = p + ylim(c(0, ymax))
  }
  
  return(p)
}

plot_snp_ratios("Ncoa3")

plot_snp_ratios("Top1")

plot_snp_ratios("Top2a")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
########            A N A L Y S I S    o f    T O P   S N P S         ##########

top_snp_pick_count = 50

snp_totals_2cell = 
detailed_riboseq_table %>%
  filter(group == '2cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )


snp_totals_4cell = 
  detailed_riboseq_table %>%
  filter(group == '4cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )

snp_totals_8cell = 
  detailed_riboseq_table %>%
  filter(group == '8cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )

transcripts_with_at_least_10_reads = intersect(snp_totals_8cell$transcript, snp_totals_4cell$transcript)
transcripts_with_at_least_10_reads = intersect(transcripts_with_at_least_10_reads, snp_totals_2cell$transcript)


# We decided not to exclude those transcripts
# mii_0_paternal_transcripts = 
#   unique(
#       (
#         detailed_table %>% filter(group == "MII") %>%
#         filter(type == "ribo") %>%  
#       group_by(transcript, group) %>%
#       mutate(paternal_ratio = sum(paternal) /  sum(paternal + maternal) ) %>%
#       filter( paternal_ratio < 0.05 )
#       ) $transcript )


#top_transcripts_combined = intersect(transcripts_with_at_least_10_reads, mii_0_paternal_transcripts)
top_transcripts_combined = transcripts_with_at_least_10_reads


top_snps_df_ribo = detailed_riboseq_table %>%
  filter( group %in% c('1cell', '2cell', '4cell', '8cell') ) %>%
  filter(transcript %in% top_transcripts_combined) %>%
  group_by( group, transcript ) %>%
  mutate( paternal_ratio = sum(paternal) / sum( paternal + maternal) )

top_snps_df_rna = detailed_rnaseq_table %>%
  filter( group %in% c('1cell', '2cell', '4cell', '8cell') ) %>%
  filter(transcript %in% top_transcripts_combined) %>%
  group_by( group, transcript ) %>%
  mutate( paternal_ratio = sum(paternal) / sum( paternal + maternal) )


top_snps_paternal_ratios_ribo = 
  top_snps_df_ribo %>%
    group_by(group, transcript, paternal_ratio) %>%
    summarise() %>% 
    group_by(group) %>% arrange(transcript)

top_snps_paternal_ratios_rna = 
  top_snps_df_rna %>%
  group_by(group, transcript, paternal_ratio) %>%
  summarise() %>% 
  group_by(group) %>% arrange(transcript)

top_snps_paternal_ratio_difference = data.frame( 
  group                     = top_snps_paternal_ratios_ribo$group,
  transcript                = top_snps_paternal_ratios_ribo$transcript,
  paternal_ratio_difference = top_snps_paternal_ratios_ribo$paternal_ratio - top_snps_paternal_ratios_rna$paternal_ratio )




#### Difference Dataframe for Heatmap
top_snps_paternal_ratio_difference_wide            = dcast( top_snps_paternal_ratio_difference,  
                                                            transcript ~ group ,
                                                            value.var = "paternal_ratio_difference" )

row.names(top_snps_paternal_ratio_difference_wide) = top_snps_paternal_ratio_difference_wide$transcript

top_snps_paternal_ratio_difference_wide            = subset(top_snps_paternal_ratio_difference_wide, 
                                                            select = c("2cell","4cell","8cell"))



top_snps_paternal_ratios_ribo_wide = dcast( top_snps_paternal_ratios_ribo,  transcript ~ group ,value.var = "paternal_ratio" )

row.names(top_snps_paternal_ratios_ribo_wide) = top_snps_paternal_ratios_ribo_wide$transcript

top_snps_paternal_ratios_ribo_wide = subset( top_snps_paternal_ratios_ribo_wide,
                                            select = c("2cell","4cell","8cell") )


top_snps_paternal_ratios_rna_wide = dcast( top_snps_paternal_ratios_rna,  transcript ~ group ,value.var = "paternal_ratio" )

row.names(top_snps_paternal_ratios_rna_wide) = top_snps_paternal_ratios_rna_wide$transcript

top_snps_paternal_ratios_rna_wide = subset( top_snps_paternal_ratios_rna_wide,
                                            select = c("2cell","4cell","8cell") )


top_snps_paternal_ratios_rna_wide[is.na(top_snps_paternal_ratios_rna_wide) ] = 0



difference_heatmap_color = colorRampPalette(c(RATIO_RNA_COLOUR, "white", RATIO_RIBO_COLOUR))(paletteLength)
ribo_heatmap_color       = colorRampPalette(c("white", RATIO_RIBO_COLOUR))(paletteLength)
rna_heatmap_color        = colorRampPalette(c("white", RATIO_RNA_COLOUR))(paletteLength)

difference_heatmap_breaks = c(seq( -1 ,              0, length.out = ceiling(paletteLength/2) + 1), 
              seq( 1 /paletteLength, 1, length.out = floor(paletteLength/2)))

ribo_breaks = c(seq(0, 1, length.out = paletteLength ))


prop_test_significant_genes  = row.names(heatmap_binary)
prop_test_intersect_top_snps = intersect(top_transcripts_combined, prop_test_significant_genes)



top_snp_markers            = matrix("", nrow = dim(top_snps_paternal_ratios_ribo_wide)[1], ncol = dim(top_snps_paternal_ratios_ribo_wide)[2] )
row.names(top_snp_markers) = row.names(top_snps_paternal_ratios_ribo_wide)
col.names(top_snp_markers) = col.names(top_snps_paternal_ratios_ribo_wide)

for(g in prop_test_intersect_top_snps){
  top_snp_markers[g,] = heatmap_binary[g,c(2,3,4)]
}


### ribo_paternal_percentage
top_snps_paternal_ratio_ribo_plot = 
pheatmap( top_snps_paternal_ratios_ribo_wide , 
          #clustering_method = "median",
          #clustering_distance_rows = "correlation",
          #clustering_method = "centroid",
          #clustering_method = "complete",
          clustering_method = "ward.D",
          show_rownames     = TRUE,
          display_numbers   = top_snp_markers,
          cellwidth         = 40,
          #treeheight_row    = 100,
          treeheight_col    = 0,
          cluster_cols      = FALSE, 
          cluster_rows      = TRUE, 
          color             = ribo_heatmap_color,
          breaks            = ribo_breaks,
          na_col            = "white",
          angle_col         = 0,
          labels_col        = c( "2cell", "4cell", "8cell"),
          main              = "Paternal Ratio Ribo",
          fontsize          = 10,
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 8,
          fontsize_number   = FONT_LABEL_SIZE
)

ribo_heatmap_row_order = top_snps_paternal_ratio_ribo_plot$tree_row$order

# rna_paternal_percentage
top_snps_paternal_ratio_rna_plot = 
pheatmap( top_snps_paternal_ratios_rna_wide[ribo_heatmap_row_order, ] , 
          #clustering_method = "median",
          #clustering_distance_rows = "correlation",
          #clustering_method = "centroid",
          #clustering_method = "complete",
          #clustering_method = "ward.D",
          show_rownames     = TRUE,
          display_numbers   = top_snp_markers[ribo_heatmap_row_order, ],
          #cutree_rows       = 4,
          cellwidth         = 40,
          treeheight_row    = 50,
          treeheight_col    = 0,
          cluster_cols      = FALSE, 
          cluster_rows      = FALSE, 
          color             = rna_heatmap_color,
          breaks            = ribo_breaks,
          na_col            = "white",
          angle_col         = 0,
          labels_col        = c( "2cell", "4cell", "8cell"),
          main              = "Paternal Ratio RNA",
          fontsize          = 10,
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 8,
          fontsize_number   = FONT_LABEL_SIZE,
          number_color      = RATIO_RIBO_COLOUR 
          
)


### Difference between ribo and rna ratios
### ribo_paternal_percentage - rna_paternal_ercentage
top_snps_paternal_ratio_difference_plot = 
  pheatmap( top_snps_paternal_ratio_difference_wide[ribo_heatmap_row_order,] , 
            #clustering_method = "median",
            #clustering_distance_rows = "correlation",
            #clustering_method = "centroid",
            #clustering_method = "complete",
            #clustering_method = "ward.D",
            show_rownames     = TRUE,
            #cutree_rows       = 4,
            display_numbers   = top_snp_markers[ribo_heatmap_row_order, ],
            cellwidth         = 40,
            treeheight_row    = 0,
            treeheight_col    = 0,
            cluster_cols      = FALSE, 
            cluster_rows      = FALSE, 
            color             = difference_heatmap_color,
            breaks            = difference_heatmap_breaks,
            na_col            = "white",
            angle_col         = 0,
            labels_col        = c( "2cell", "4cell", "8cell"),
            main              = "Paternal Ratio Difference",
            fontsize          = 10,
            fontsize_col      = FONT_LABEL_SIZE,
            fontsize_row      = 8,
            fontsize_number   = FONT_LABEL_SIZE
  )




supplementary_heatmap_grid_plot = 
  plot_grid(top_snps_paternal_ratio_ribo_plot[[4]], 
            top_snps_paternal_ratio_rna_plot[[4]] ,
            top_snps_paternal_ratio_difference_plot[[4]],
            ncol  = 3, 
            align = 'hv')

supplementary_heatmap_grid_plot



################################################################################
################################################################################
################################################################################
################################################################################
#########            P D F    P R O D U C T I O N                      #########

get_output_file_path = function(file_name, folder = output_folder){
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

combine_gene_detail_and_ratio_plots = function(gene){
  line_plot     = plot_snp_ratios(gene)
  detailed_plot = plot_snp_gene_detailed(gene)
  
  row_spacer = generate_blank_plot()
  
  this_plot = plot_grid(line_plot, row_spacer, detailed_plot, 
                        rel_heights = c(1,0.05,1),
                        ncol = 1,
                        align = "h")
  
  return(this_plot)
}

combine_gene_detail_and_ratio_plots("Ncoa3")

################################################################################

save_plot_pdf("supplementary_allele_percentages_ribo.pdf", supplementary_allele_percentages_ribo)

save_plot_pdf("supplementary_allele_percentages_rnaseq.pdf", supplementary_allele_percentages_rnaseq)

save_plot_pdf("main_paternal_percentage_figure.pdf", main_paternal_percentage_figure)

save_plot_pdf("main_heatmap_figure.pdf", main_heatmap_figure)

save_plot_pdf("supplementary_heatmap_figure.pdf", supplementary_heatmap_figure, width = 5, height = 12)

save_plot_pdf("supplementary_heatmap_grid_plot.pdf", supplementary_heatmap_grid_plot, width = 15, height = 12)

################################################################################


combine_gene_detail_and_ratio_plots("Cdt1")

#################################################

interesting_genes = c(
  "Ncoa3", # !!! High paternal ratio at 4cell stage, good reproducibility
  "Rpl38", #-> Small effect size; large counts
  "Tsr1",  #  -> rRNA accumulation associated; strange behavior that doesn’t follow the clear patterns
  "Rps6",  # -> Lagging translation with large counts; paternal activation by 2-cell stage; paternal translation slowly starts at 4-cell and matches by 8-cell stage
  "Eef1b2",  #-> Same story of translation lagging the RNA expression. 
  "Polr1e",
  "Eif3d", # -> Lagged translation
  "Rpl21", #-> Not super interesting but maybe lower translation of the paternal copy
  "Eif3g", #-> RNA degraded in the GV to MII transition; both copies activated by 2-cell stage; translation of maternal copy looks higher in 4-cell stage but likely not a great one to highlight. 
  "Ltv1", # -> Lagged; ribosome biogenesis factor
  "Rpl28", # -> Much stronger translation of the paternal copy
  "Eif4b", # -> Lagging
  "Ddx21", # -> Lagging translation
  "Rpl28", # ->Lagging translation
  "Znhit3", # -> Lagging translation; implicated in pre-ribosomal RNA processing
  "Eprs", # !!!-> Paternal copy seems to be lowly tranlsated consistently at both 4- and 8- cell stage
  "Lyar", # -> Nucleolar gene

  "Tpx2", # - Maternal gene with lagging translation of paternal copy; involved in cell cycle
  "Pa2g4", #  - Not maternal; involved in proliferation; interesting asymmetric transcription activation at 2cell stage which reaches 50:50 in RNA by 4-cell but translation lags
  "Nin", # - Maternal gene; decent read counts; member of pericentriolar matrix
  "Cdc42", # -> lagging translation
  "Mcm7", # !!! Paternal is consistently high at 4 and 8 cell stages
  "Cdk1", # Higher paternal reads at 4 cell stage
  "Ccnb1",
  "Ccnh", # Very high paternal ratio at 8 cell stage
  "Brca2", # -> Might be more interesting in a stage-specific manner
  "Cct6a", #-> Consistently high translation of the paternal copy starting at 2-cell (Note: paternal ratios are almost equal at cell)
  "Dnmt1", #!!! -> Maternal transcript but small expression from paternal copy with exclusive translation of the paternal copy
  "Pgd", # Not much known about in the context of development but has robust counts: Maybe not super interesting, high var.
  "Ppat", # Seems like there are a number of metabolic enzymes in the list but nothing super obvious in terms of known biology
  "Pemt",
  "Folr1", # -> Good counts but there is a fetal version Folr2. Shouldn’t that be more abundant?
  "Zfp296",
  "Lclat1", # Paternal ratio decreases at 8cell
  "Stip1",
  "Gpd1l", ## Almost same ratios on all cases
  "Top1", # Almost same ratios
  "Rsl1d1", #Almost same ratios, decent counts
  "Pttg1", # Almost same ratios, large counts
  "Kdm5b", # goes flat in both cases. High variation though
  "Slc16a6",  # goes flat in both cases. High variation though
  "Cdt1" # Same behavior in both cases
)

### Save gene plots into pdf files

for(g in interesting_genes){
  this_file = paste( "gene_", g, ".pdf", sep = ""  )
  save_plot_pdf(this_file, combine_gene_detail_and_ratio_plots(g), width = 12, height = 12)
}



