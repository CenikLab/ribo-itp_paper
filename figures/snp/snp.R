#SNP Figure

riboseq_overall_percentage_file = "../../snp/notebooks/snp_dataframes/riboseq_snp_percentages.csv"
rnaseq_overall_percentage_file  = "../../snp/notebooks/snp_dataframes/rnaseq_snp_percentages.csv"

overall_ribo_to_rna_ratios_file = "../../snp/notebooks/snp_dataframes/refined_ribo_to_rna_ratios.csv"
overall_ribo_to_rna_binary_file = "../../snp/notebooks/snp_dataframes/refined_ribo_to_rna_binary_significance.csv"



detailed_riboseq_count_file             = "../../snp/notebooks/snp_dataframes/riboseq_detailed_snps.csv.gz"
detailed_rnaseq_count_file              = "../../snp/notebooks/snp_dataframes/rnaseq_detailed_snps.csv.gz"

output_folder           = "./pdfs"

################################################################################
########                    L I B R A R I E S                          #########

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

library(Cairo)



#Heatmap related packages
library(pheatmap)

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

#FOURCELL_BACKGROUND_COLOR = "#d9d9d9"

TWOCELL_BACKGROUND_COLOR = "#d9d9d9"
TWOCELL_BACKGROUND_COLOR  = "gray95"
FOURCELL_BACKGROUND_COLOR = "gray90"
EIGHTCELL_BACKGROUND_COLOR = "gray80"

COUNT_NORMALIZATION_FACTOR = 10000

#RATIO_RIBO_COLOUR = "orange"
#RATIO_RNA_COLOUR  = "blue"

RATIO_RIBO_COLOUR = BURNT_ORANGE
RATIO_RNA_COLOUR  = "navy"

PERCENTAGE_DASHED_COLOR = "#088f99"

heatmap_ribo_orange = rgb(228,88,10 , maxColorValue = 255)
heatmap_ribo_blue   = rgb(55,135,192, maxColorValue = 255)

################################################################################
#########                 F O N T   S I Z E S                          #########

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 9

PDF_resolution = 600
FIGURE_FONT = "helvetica"

DETAILED_BARPLOT_FONT_SIZE = 7
GENE_NAME_FONT_SIZE        = 8

################################################################################
##### Experiments



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
        legend.key.size  = unit(0.15, 'in'),
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
        legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        legend.key.size  = unit(0.15, 'in'),)  + 
  scale_fill_manual(values = PERCENTAGE_BARPLOT_COLORS ) + 
  geom_hline(yintercept = 50, linetype="dotted", 
             color = PERCENTAGE_DASHED_COLOR, size=0.6) + 
  labs(title = "Allele Percentages of RNA-Seq")


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
  
  labs(title = "", 
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
  geom_bar(position = "dodge", stat="identity", alpha = 0.8 ) +
  geom_errorbar( aes(  x    = stage, 
                       ymin = average_percentage - sd_percentage, 
                       ymax = average_percentage + sd_percentage,
                       color = (experiment_type)), 
                 width    = 0.4, 
                 alpha    = 1, 
                 size     = 0.4,
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
        legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.line.y      = element_line(colour = "black", size = 0.35), 
        legend.key.size  = unit(0.18, "in"),
        legend.position  = c(0.3,0.8)) + 
  
  scale_y_continuous( expand = c(0, 0
                                 ), limits = c(-0.5, 50) ) +
  scale_fill_manual( name = "experiment_type", values = c(heatmap_ribo_orange,
                                                          heatmap_ribo_blue) ) + 
  scale_color_manual(name = "experiment_type",  values = c(heatmap_ribo_orange,
                                                           heatmap_ribo_blue) ) + 
  labs(title = "", 
       x     = "Stage", 
       y     = "Percentage")

main_paternal_percentage_figure



standalone_paternal_percentage_figure = 
  ggplot(data=corrected_percentage_averages, aes(x=stage, y=average_percentage, fill = experiment_type )  )  + 
  geom_bar(position = "dodge", stat="identity", alpha = 0.8 ) +
  geom_errorbar( aes(  x    = stage, 
                       ymin = average_percentage - sd_percentage, 
                       ymax = average_percentage + sd_percentage,
                       color = (experiment_type)), 
                 width    = 0.4, 
                 alpha    = 1, 
                 size     = 0.4,
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
        legend.text      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.line.y      = element_line(colour = "black", size = 0.35), 
        legend.key.size  = unit(0.18, "in"),
        ) + 
  
  scale_y_continuous( expand = c(0, 0
  ), limits = c(-0.5, 50) ) +
  scale_fill_manual( name = "experiment_type", values = c(RATIO_RIBO_COLOUR,
                                                          RATIO_RNA_COLOUR) ) + 
  scale_color_manual(name = "experiment_type",  values = c(RATIO_RIBO_COLOUR,
                                                           RATIO_RNA_COLOUR) ) + 
  labs(title = "Paternal Allele Percentages", 
       x     = "Stage", 
       y     = "Percentage")

standalone_paternal_percentage_figure

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
#myColor       = colorRampPalette(c(RATIO_RNA_COLOUR, "white", RATIO_RIBO_COLOUR))(paletteLength)
myColor       = colorRampPalette(c(heatmap_ribo_blue, "white", heatmap_ribo_orange  ))(paletteLength)
#myColor       = colorRampPalette(c("navy", "white", "#f8971f"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(heatmap_df_raw, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1), 
#              seq(max(heatmap_df_raw, na.rm = TRUE)/paletteLength, max(heatmap_df_raw, na.rm = TRUE), length.out=floor(paletteLength/2)))

myBreaks <- c(seq( -1 ,              0, length.out = ceiling(paletteLength/2) + 1), 
              seq( 1 /paletteLength, 1, length.out = floor(paletteLength/2)))


clustering_order = c(
  'Nop14',
  'Tmppe', 
  'Slc13a2',
  'Ppp2ca',
  'Srpk1',
  'Cbx3',
  'Ncoa3',
  'Cdk1',
  'Baz1a',
  'Dyrk3',
  'Lclat1',
  'Lyar',
  'Umps',
  'Tsen2',
  'Ccnh',
  
  'Folr1',
  'Pa2g4',
  'Zfp296',
  'Mrps9',
  'Eif3d',
  'Nin',
  'Ddx21',
  'Bcat1',
  'Mysm1'
)

ordered_heatmap_df =  heatmap_df[match(clustering_order,  row.names(heatmap_df) ), ]



NUMBER_OF_HEATMAP_CLUSTERS = 4
### Trimmed down
### For the main figure
old_main_heatmap_figure = 
pheatmap( t( ordered_heatmap_df) , 
        # clustering_method = "median",
        #clustering_distance_rows = "correlation",
        #clustering_method = "centroid",
        #clustering_method = "complete",
        clustering_method = "ward.D",
        show_rownames     = TRUE,
        cutree_cols       = NUMBER_OF_HEATMAP_CLUSTERS,
        cellwidth         = 15,
        treeheight_row    = 0,
        treeheight_col    = 0,
        cluster_cols      = FALSE, 
        cluster_rows      = FALSE,
        color             = myColor,
        breaks            = myBreaks,
        na_col            = "white",
        angle_col         = 90,
        labels_row        = c("1cell", "2cell", "4cell", "8cell"),
        main              = "",
        fontsize          = 9,
        fontsize_col      = FONT_LABEL_SIZE,
        fontsize_row      = 6,
        fontsize_number   = FONT_LABEL_SIZE,
        )

old_main_heatmap_figure

#### Detailed with row names
#### For the supplementary Figure
standalone_heatmap_figure = 
pheatmap( heatmap_df , 
          #clustering_method = "centroid",
          #clustering_method = "median",
          #clustering_method = "complete",
          clustering_method = "ward.D",
          legend_breaks     = c(-1,  -0.5, 0, 0.5, 1),
          legend            = T,
          show_rownames     = TRUE,
          cellwidth         = 25,
          cutree_rows       = NUMBER_OF_HEATMAP_CLUSTERS,
          treeheight_row    = 0,
          cluster_cols      = FALSE, 
          color             = myColor,
          breaks            = myBreaks,
          na_col            = "white",
          angle_col         = 0,
          labels_col        = c("1cell", "2cell", "4cell", "8cell"),
          main              = "Paternal Ratio Difference",
          display_numbers   = heatmap_binary,
          fontsize          = 9,
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 7,
          fontsize_number   = FONT_LABEL_SIZE
          )

standalone_heatmap_figure


standalone_heatmap_figure_no_stars = 
  pheatmap( heatmap_df , 
            #clustering_method = "centroid",
            #clustering_method = "median",
            #clustering_method = "complete",
            clustering_method = "ward.D",
            legend_breaks     = c(-1,  -0.5, 0, 0.5, 1),
            legend            = T,
            show_rownames     = TRUE,
            cellwidth         = 25,
            cutree_rows       = NUMBER_OF_HEATMAP_CLUSTERS,
            treeheight_row    = 0,
            cluster_cols      = FALSE, 
            color             = myColor,
            breaks            = myBreaks,
            na_col            = "white",
            angle_col         = 0,
            labels_col        = c("1cell", "2cell", "4cell", "8cell"),
            main              = "Paternal Ratio Difference",
            fontsize          = 9,
            fontsize_col      = FONT_LABEL_SIZE,
            fontsize_row      = 7,
            fontsize_number   = FONT_LABEL_SIZE
  )

standalone_heatmap_figure_no_stars


################################################################################
####   N E W   V E R S I O N   O F   T H E   M A I N   H E A T M A P   #########

selected_snps_pseudo_count = 1

riboseq_paternal_ratios_selected = 
  detailed_riboseq_table %>%
    filter( group %in% c( '2cell', '4cell', '8cell') ) %>%
    filter(transcript %in% clustering_order) %>%
    group_by( group, transcript ) %>%
    summarise( paternal_ratio =  ( sum(paternal) ) / (sum( paternal + maternal) + 1*selected_snps_pseudo_count ) )

rnaseq_paternal_ratios_selected = 
  detailed_rnaseq_table %>%
  filter( group %in% c( '2cell', '4cell', '8cell') ) %>%
  filter(transcript %in% clustering_order) %>%
  group_by( group, transcript ) %>%
  summarise( paternal_ratio =  ( sum(paternal)) / (sum( paternal + maternal) + 1*selected_snps_pseudo_count ) )

riboseq_paternal_ratios_selected_wide_unordered = 
  dcast( riboseq_paternal_ratios_selected,  
         transcript ~ group ,
         value.var = "paternal_ratio" )

rnaseq_paternal_ratios_selected_wide_unordered = 
  dcast( rnaseq_paternal_ratios_selected,  
         transcript ~ group ,
         value.var = "paternal_ratio" )

riboseq_paternal_ratios_selected_wide = 
  riboseq_paternal_ratios_selected_wide_unordered[match( clustering_order, riboseq_paternal_ratios_selected_wide_unordered$transcript), ]
rnaseq_paternal_ratios_selected_wide = 
rnaseq_paternal_ratios_selected_wide_unordered[match( clustering_order, rnaseq_paternal_ratios_selected_wide_unordered$transcript), ]

grey_color_map = colorRampPalette(c("white", "green" ,"black"  ))(paletteLength)

RColorBrewer::brewer.pal(9, "BuGn")

grey_color_map = colorRampPalette( RColorBrewer::brewer.pal(9, "BuGn") )(paletteLength)

grey_color_map = colorRampPalette( RColorBrewer::brewer.pal(9, "Greys") )(paletteLength)

riboseq_paternal_ratios_heatmap = 
  pheatmap( t( riboseq_paternal_ratios_selected_wide[, -1]) , 
            # clustering_method = "median",
            #clustering_distance_rows = "correlation",
            #clustering_method = "centroid",
            #clustering_method = "complete",
            clustering_method = "ward.D",
            show_rownames     = TRUE,
            cutree_cols       = NUMBER_OF_HEATMAP_CLUSTERS,
            cellwidth         = 15,
            treeheight_row    = 0,
            treeheight_col    = 0,
            cluster_cols      = FALSE, 
            cluster_rows      = FALSE,
            color             = grey_color_map,
            legend_breaks            = c(0, 0.5, 1),
            breaks = seq(0,1, length.out = 100),
            na_col            = "white",
            angle_col         = 90,
            labels_row        = c( "2cell", "4cell", "8cell"),
            labels_col        = riboseq_paternal_ratios_selected_wide[,1], 
            main              = "",
            fontsize          = 9,
            fontsize_col      = FONT_LABEL_SIZE,
            fontsize_row      = 6,
            fontsize_number   = FONT_LABEL_SIZE,
  )


rnaseq_paternal_ratios_heatmap = 
  pheatmap( t( rnaseq_paternal_ratios_selected_wide[, -1]) , 
            # clustering_method = "median",
            #clustering_distance_rows = "correlation",
            #clustering_method = "centroid",
            #clustering_method = "complete",
            clustering_method = "ward.D",
            show_rownames     = TRUE,
            cutree_cols       = NUMBER_OF_HEATMAP_CLUSTERS,
            cellwidth         = 15,
            treeheight_row    = 0,
            treeheight_col    = 0,
            cluster_cols      = FALSE, 
            cluster_rows      = FALSE,
            color             = grey_color_map,
            legend_breaks            = c(0, 0.5, 1),
            breaks = seq(0,1, length.out = 100), 
            na_col            = "white",
            angle_col         = 90,
            labels_row        = c( "2cell", "4cell", "8cell"),
            labels_col        = rnaseq_paternal_ratios_selected_wide[, 1], 
            main              = "",
            fontsize          = 9,
            fontsize_col      = FONT_LABEL_SIZE,
            fontsize_row      = 6,
            fontsize_number   = FONT_LABEL_SIZE,
  )



ribo_rna_paternal_ratios_combined = cbind( riboseq_paternal_ratios_selected_wide, rnaseq_paternal_ratios_selected_wide[, -1] )
View(ribo_rna_paternal_ratios_combined)


## We can use the clustering here for hte paper.
pheatmap( t( ribo_rna_paternal_ratios_combined[, -1]) , 
          # clustering_method = "median",
          clustering_distance_cols = "correlation",
          #clustering_method = "centroid",
          #clustering_method = "complete",
          #clustering_method = "ward.D",
          show_rownames     = TRUE,
          cutree_cols       = NUMBER_OF_HEATMAP_CLUSTERS + 2,
          cellwidth         = 15,
          treeheight_row    = 0,
          treeheight_col    = 0,
          cluster_cols      = TRUE, 
          cluster_rows      = FALSE,
          color             = grey_color_map,
          legend_breaks            = c(0, 0.5, 1),
          breaks = seq(0,1, length.out = 100), 
          na_col            = "white",
          angle_col         = 90,
          labels_row        = c( "2cell", "4cell", "8cell", "2cell", "4cell", "8cell"),
          labels_col        = ribo_rna_paternal_ratios_combined[, 1], 
          main              = "",
          fontsize          = 9,
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 6,
          fontsize_number   = FONT_LABEL_SIZE,
)

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


make_y_axis = function(ymax){
  if(ymax == 5){
    result = c(0, 5)
  }
  else if(ymax %% 2 == 0){
    result = c(0, ymax / 2, ymax)
  } 
  else if(ymax %% 3 == 0){
    result = c(0, ymax / 3, 2*( ymax / 3 ) , ymax)
  } 
  else{
    result = c(0, ymax)
  }
  
  return(result)
}

################################################################################


generate_legend_unit = function(gene, experiment_type){
  
  gene_data = 
    detailed_table_with_dummies %>%
    filter( transcript == gene & type == experiment_type ) %>%
    filter( group == "2cell" | group == "4cell" | group == "8cell") %>%
    filter(paternal_per_k > 0 | maternal_per_k > 0)
  
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
    geom_bar(stat="identity", width= 0.8, color = "black", size = 0.3) + 
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
      plot.title = element_text(color = title_colors[[experiment_type]], size=7, face="bold", hjust = 0.5),
    ) + 
    scale_fill_manual(values = barplot_colors)  + 
    #ggtitle(title_texts[[experiment_type]])
    ggtitle( paste( title_texts[[experiment_type]], as.character(number_of_snps), sep = " " ) )
  
  #if(experiment_type == "rna"){
  #  p = p + scale_y_reverse()
  #}
  
  return(p + rotate())
}

generate_legend_unit("Ncoa3", "ribo")

generate_legend = function(gene){
  p_ribo = generate_legend_unit(gene, experiment_type = "ribo")
  p_rna  = generate_legend_unit(gene, experiment_type = "rna")
  
  this_blank_pad  = generate_blank_plot()
  this_legend_pre = plot_grid( this_blank_pad, p_ribo, this_blank_pad, p_rna, 
                           nrow = 1,
                           rel_widths = c(0.2, 1 , 0.2, 1 ) )
  this_legend     = plot_grid(this_blank_pad, this_legend_pre, this_blank_pad,
                              ncol = 1, rel_heights = c(0.01, 1, 0.01)) 
  
  return(this_legend)
}

snp_detailed_legend = generate_legend("Nin")

snp_detailed_legend 
generate_legend("Ncoa3")

################################################################################

this_gene = "Ncoa3"

a = 
detailed_table_with_dummies %>%
  filter(group == "2cell" | group == "4cell" | group == "8cell") %>%
  filter(transcript == this_gene) %>% 
  group_by(position) %>%
  filter( sum(paternal_per_k + maternal_per_k) > 0  )


################################################################################
################################################################################
###### Adding color values to the SNP table

get_snp_positions_of_gene = function(gene){
  this_df = detailed_table_with_dummies %>% 
    filter(experiment == "4cell-1" & type == "ribo" ) %>%
    filter(transcript == gene)
  
  return( sort(this_df$position) )
}

a = get_snp_positions_of_gene("Ncoa3")
a

sort(a$position)

produce_color_map_for_gene_positions = function( gene_positions, type){
  
  if(type == "ribo"){
    this_color_name = "Oranges"
  }
  else if(type == "rna"){
    this_color_name = "Blues"
  } 
  else{
    print("Error! Wrong Type!")
    return(NULL)
  }
  
  number_of_snps                 = length(gene_positions)
  this_ribo_color_palette        = (colorRampPalette(brewer.pal(9,this_color_name))( number_of_snps + 5 )[-seq( 1,  5) ])
  names(this_ribo_color_palette) = gene_positions
  
  return( rev(this_ribo_color_palette) )
}

get_gene_vector = function(){
  this_df = detailed_table_with_dummies %>% 
    filter(experiment == "4cell-1" & type == "ribo" ) %>%
    summarise(transcript)
    
  return( sort(unique(this_df$transcript)) )
  
}

produce_color_hash_for_gene = function(gene, type){
  
  snp_positions  = get_snp_positions_of_gene(gene)
  this_color_map = produce_color_map_for_gene_positions(snp_positions, type)
    
  return(this_color_map)
}

a = produce_color_hash_for_gene("Nin", "ribo")
a


k = produce_color_map_for_gene_positions(a, "rna")

typeof(k)

b = 
  detailed_table_with_dummies %>% 
  filter(experiment == "2cell-1" & type == "rna" ) %>%
  filter(transcript == "Ncoa3")

sort(b$position)

c = 
  detailed_table_with_dummies %>% 
  filter(experiment == "8cell-1" & type == "rna" ) %>%
  filter(transcript == "Ncoa3")

sorted_snp_positions = sort(c$position)
length(sorted_snp_positions)

number_of_snps           = length( sorted_snp_positions  )

this_blues_color_palette = colorRampPalette(brewer.pal(9,"Blues"))( number_of_snps + 5 )[-seq( 1,  5) ]

this_ribo_color_mapper = vector("list", number_of_snps)

names(this_ribo_color_mapper) = sorted_snp_positions

for(i in 1:number_of_snps){
  this_ribo_color_mapper[i] = this_ribo_color_palette[i]
}
this_ribo_color_mapper
################################################################################
### U N I T     P L O T 

unit_plot = function(gene, experiment_type, experiment_group, allele_type, ymax){
  
  bar_palettes = list( ribo = "Oranges", rna = "Blues"  )
  
  y_scale_expand = c(0,0)
  ymax           = adjust_ymax(ymax)
  
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
  else if(experiment_group == "2cell"){
    background_canvas = element_rect(fill = TWOCELL_BACKGROUND_COLOR)  
  }
  else if(experiment_group == "8cell"){
    background_canvas = element_rect(fill = EIGHTCELL_BACKGROUND_COLOR)  
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
      axis.text.y = element_text(family = FIGURE_FONT, face = "plain", size = DETAILED_BARPLOT_FONT_SIZE),
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

################################################################################
################################################################################
unit_plot_alternative = function(gene, experiment_type, experiment_group, allele_type, ymax){
  ## We won't use this
  ## We generate the fill colors for each transcript anually
  ## for sanity check reaosns
  bar_palettes = list( ribo = "Oranges", rna = "Blues"  )
  
  y_scale_expand = c(0,0)
  ymax           = adjust_ymax(ymax)
  
  gene_data = 
    detailed_table_with_dummies %>%
    filter(group == "2cell" | group == "4cell" | group == "8cell") %>%
    filter(transcript == gene) %>% 
    group_by(position) %>%
    #filter( sum(paternal_per_k + maternal_per_k) > 0  ) %>%
    filter( type == experiment_type ) %>%
    filter( group == experiment_group)
  
  number_of_snps = length( unique(gene_data$position  ) )
  
  #paternal_color_palette = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps + 5 )[-seq( 1,  5) ]
  #paternal_color_palette = colorRampPalette(brewer.pal(9,bar_palettes[[experiment_type]]))( number_of_snps  ) 
  this_fill_colors = produce_color_hash_for_gene(gene, experiment_type)
  
  yticks = make_y_axis(ymax)
  
  if(experiment_group == "4cell"){
    background_canvas = element_rect(fill = FOURCELL_BACKGROUND_COLOR)
  }
  else if(experiment_group == "2cell"){
    background_canvas = element_rect(fill = TWOCELL_BACKGROUND_COLOR)  
  }
  else if(experiment_group == "8cell"){
    background_canvas = element_rect(fill = EIGHTCELL_BACKGROUND_COLOR)  
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
      axis.text.y = element_text(family = FIGURE_FONT, face = "plain", size = DETAILED_BARPLOT_FONT_SIZE),
      # axis.ticks.y=element_blank()
    ) +
    scale_fill_manual(values = this_fill_colors) 
  
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
  
  plot_4cell_raw =  plot_experiment_group(gene = gene, experiment_group = "4cell")  
  
  plot_8cell_raw = plot_experiment_group(gene = gene, experiment_group = "8cell")

  paternal_label = ggdraw() + 
    draw_label(
      "paternal c.",
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
      "maternal c.",
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
  
  # OLD WORKING
  #snp_detailed_legend = generate_legend(gene)
  
  # raw_plot = plot_grid( x_label,    plot_2cell,  separator_1, 
  #                       plot_4cell, separator_1, 
  #                       plot_8cell, separator_1, snp_detailed_legend,
  #                       rel_widths  = c(0.2, 1, 0.2, 1, 0.2, 1, 0.1, 0.4)  , 
  #                       rel_heights =  c(0.8, 1, 1, 1, 1, 1, 1, 0.01),
  #                       ncol        = 8  )
  
  
  raw_plot = plot_grid( x_label,    plot_2cell_raw,  separator_1, 
                        plot_4cell_raw, separator_1, 
                        plot_8cell_raw,
                        rel_widths  = c(0.2, 1, 0.02, 1, 0.02, 1 )  , 
                        rel_heights =  c(0.8, 1, 1, 1, 1, 1),
                        ncol        = 6  )
  
  title_main = ggdraw() + 
    draw_label(
      gene,
      fontface   = 'plain',
      size       = GENE_NAME_FONT_SIZE ,
      fontfamily = 'helvetica',
    ) 
    # theme(
    #   # add margin on the left of the drawing canvas,
    #   # so title is aligned with left edge of first plot
    #   plot.margin = margin(0, 0, 0, 0)
    # )
  
  this_plot = plot_grid(
    title_main, raw_plot,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
  )
  
  return(raw_plot)
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
  #
  paternal_mean = mean(paternal_counts)
  total_mean    = mean(total_counts)
  
  f = (paternal_mean / total_mean)
  
  this_cov = cov( paternal_counts, total_counts )
  
  paternal_sd = sd(paternal_counts)
  total_sd    = sd(total_counts)
  
  this_sd =  f * sqrt(  (paternal_sd / paternal_mean)**2 + 
                         (total_sd/ total_mean)**2 - 
                         2*(this_cov/ (paternal_mean* total_mean) )  )
  
  return(this_sd / sqrt(length(paternal_counts)))
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

scale_decimal_digit <- function(x) sprintf("%.1f", x)

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
  
  max_ratio = max( c(ribo_ratios, rna_ratios)  )
  
  y_breaks = waiver()
  
  # if(max_ratio > 0.8){
  #   y_breaks = c(0, 0.5 , 1)
  # }
  # else{
  #   y_breaks = waiver()
  # }
  
  ribo_sd = determine_propogated_sd( gene_data_ribo  )
  rna_sd  = determine_propogated_sd(  gene_data_rna )
  
  plot_df = data.frame(  paternal_ratios = c(ribo_ratios, rna_ratios), 
                         stage           = rep( c("2cell", "4cell", "8cell") , 2), 
                         sd              = c(ribo_sd, rna_sd), 
                         type =  c( rep("ribo", 3), rep("rna", 3)  )  )
  
  p = ggplot(data=plot_df, 
             aes(x    = stage, y = paternal_ratios, group = type )  )  +
    geom_point( aes(colour = type, shape = type), position = position_dodge(width = 0.1), size = 1.8 ) +
    geom_line(aes(linetype=type, colour = type), position = position_dodge(width = 0.1) ) + 
    scale_linetype_manual(values=c("solid", "dashed")) + 
    geom_errorbar(aes(ymin=paternal_ratios-sd, ymax=paternal_ratios+sd, color = type), position = position_dodge(width = 0.1), width = 0.1, size=0.4 ) + 
    scale_colour_manual(values = c(RATIO_RIBO_COLOUR, RATIO_RNA_COLOUR)   ) + 
    theme(
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = GENE_NAME_FONT_SIZE),
      panel.background  = element_blank(),
      axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      #axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
      axis.title.x      = element_blank(),
      legend.position   = "none",
      axis.line         = element_line(colour = "black", size = 0.35), 
      #legend.title      = element_blank(),
      #legend.text        = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)
    ) +
    labs(title = gene, y = "paternal ratio") +
    scale_x_discrete( expand = c(0.05, 0.05)) +
    scale_y_continuous( breaks = y_breaks, labels = scale_decimal_digit)

  
  if(ymax > 0){
    p = p + ylim(c(0, ymax))
  }
  
  return(p)
}

plot_snp_ratios("Lyar")

plot_snp_ratios("Ncoa3")

plot_snp_ratios("Top1")

plot_snp_ratios("Top2a")

plot_snp_ratios("Rpa1")

plot_snp_ratios("Tpx2")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
########            A N A L Y S I S    o f    T O P   S N P S         ##########

top_snp_pick_count = 50

snp_totals_2cell_ribo = 
detailed_riboseq_table %>%
  filter(group == '2cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )


snp_totals_4cell_ribo = 
  detailed_riboseq_table %>%
  filter(group == '4cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )

snp_totals_8cell_ribo = 
  detailed_riboseq_table %>%
  filter(group == '8cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )

transcripts_with_at_least_10_reads_ribo = intersect(snp_totals_8cell_ribo$transcript, snp_totals_4cell_ribo$transcript)
transcripts_with_at_least_10_reads_ribo = intersect(transcripts_with_at_least_10_reads_ribo, snp_totals_2cell_ribo$transcript)

snp_totals_2cell_rna = 
  detailed_rnaseq_table %>%
  filter(group == '2cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )


snp_totals_4cell_rna = 
  detailed_rnaseq_table %>%
  filter(group == '4cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )

snp_totals_8cell_rna = 
  detailed_rnaseq_table %>%
  filter(group == '8cell'  ) %>%
  group_by(transcript) %>%
  summarise(transcript_sum = sum(paternal + maternal)   ) %>%
  filter(transcript_sum >= 10) %>%
  arrange(desc( transcript_sum) )

transcripts_with_at_least_10_reads_rna = intersect(snp_totals_8cell_rna$transcript, snp_totals_4cell_rna$transcript)
transcripts_with_at_least_10_reads_rna = intersect(transcripts_with_at_least_10_reads_rna, snp_totals_2cell_rna$transcript)

length(top_transcripts_combined)

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
top_transcripts_combined = intersect(transcripts_with_at_least_10_reads_ribo, transcripts_with_at_least_10_reads_rna)

top_snps_pseudo_count = 0

top_snps_df_ribo = detailed_riboseq_table %>%
  filter( group %in% c('1cell', '2cell', '4cell', '8cell') ) %>%
  filter(transcript %in% top_transcripts_combined) %>%
  group_by( group, transcript ) %>%
  mutate( paternal_ratio =  ( sum(paternal) + top_snps_pseudo_count) / (sum( paternal + maternal) + 2*top_snps_pseudo_count ) )

top_snps_df_rna = detailed_rnaseq_table %>%
  filter( group %in% c('1cell', '2cell', '4cell', '8cell') ) %>%
  filter(transcript %in% top_transcripts_combined) %>%
  group_by( group, transcript ) %>%
  mutate( paternal_ratio =  ( sum(paternal) + top_snps_pseudo_count) / (sum( paternal + maternal) + 2*top_snps_pseudo_count ) )


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



difference_heatmap_color = colorRampPalette(c(heatmap_ribo_blue, "white", heatmap_ribo_orange))(paletteLength)
ribo_heatmap_color       = colorRampPalette(c("white", heatmap_ribo_orange))(paletteLength)
rna_heatmap_color        = colorRampPalette(c("white", heatmap_ribo_blue))(paletteLength)

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


supp_heatmap_cell_width = 20
heatmap_title_size = 7

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
          cellwidth         = supp_heatmap_cell_width,
          treeheight_row    = 35,
          treeheight_col    = 0,
          cluster_cols      = FALSE, 
          cluster_rows      = TRUE, 
          color             = ribo_heatmap_color,
          breaks            = ribo_breaks,
          na_col            = "white",
          angle_col         = 0,
          labels_col        = c( "2cell", "4cell", "8cell"),
          main              = "Paternal Ratio Ribo",
          fontsize          = heatmap_title_size,
          legend_breaks     = c(0, 0.5, 1),
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 7,
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
          cellwidth         = supp_heatmap_cell_width,
          treeheight_row    = 20,
          treeheight_col    = 0,
          cluster_cols      = FALSE, 
          cluster_rows      = FALSE, 
          color             = rna_heatmap_color,
          breaks            = ribo_breaks,
          na_col            = "white",
          angle_col         = 0,
          labels_col        = c( "2cell", "4cell", "8cell"),
          main              = "Paternal Ratio RNA",
          fontsize          = heatmap_title_size,
          legend_breaks     = c(0, 0.5, 1), 
          fontsize_col      = FONT_LABEL_SIZE,
          fontsize_row      = 7,
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
            cellwidth         = supp_heatmap_cell_width,
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
            fontsize          = heatmap_title_size,
            fontsize_col      = FONT_LABEL_SIZE,
            fontsize_row      = 7,
            fontsize_number   = FONT_LABEL_SIZE
  )


supplementary_heatmap_grid_plot = 
  plot_grid(top_snps_paternal_ratio_ribo_plot[[4]], 
            top_snps_paternal_ratio_rna_plot[[4]] ,
            top_snps_paternal_ratio_difference_plot[[4]],
            ncol  = 3, 
            align = 'hv')

supplementary_heatmap_grid_plot


#save_plot_pdf( "heatmap_highly_translated_genes.pdf", supplementary_heatmap_grid_plot, width = 7.25, height = 7.25  )


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
                        rel_heights = c(0.9 ,0.05,1),
                        ncol = 1)
  
  return(this_plot)
}

################################################################################
### G E N E S    F O R   T H E   M A I N F I G U R E

ncoa3_detailed_plot = combine_gene_detail_and_ratio_plots("Ncoa3")
ncoa3_detailed_plot
#save_plot_pdf("ncoa3.pdf", ncoa3_detailed_plot, width = 3.5, height = 3.5)

eif3d_detailed_plot = combine_gene_detail_and_ratio_plots("Eif3d")
eif3d_detailed_plot
#save_plot_pdf("eif3d.pdf", eif3d_detailed_plot, width = 3.5, height = 3.5)

hsp90ab1_detailed_plot = combine_gene_detail_and_ratio_plots("Hsp90ab1")
hsp90ab1_detailed_plot
#save_plot_pdf("hsp90ab1.pdf", hsp90ab1_detailed_plot, width = 3.5, height = 3.5)

folr1_detailed_plot = combine_gene_detail_and_ratio_plots("Folr1")
folr1_detailed_plot
#save_plot_pdf("folr1.pdf", folr1_detailed_plot, width = 3.5, height = 3.5)
################################################################################

################################################################################
###             SNP Legends for the main figure                             ####

snp_legend_width  = 2
snp_legend_height = 0.5

ncoa3_snp_legend = generate_legend("Ncoa3")
ncoa3_snp_legend
save_plot_pdf("ncoa3_legend.pdf", ncoa3_snp_legend, width = snp_legend_width, height = snp_legend_height)

eif3d_snp_legend = generate_legend("Eif3d")
eif3d_snp_legend
save_plot_pdf("eif3d_legend.pdf", eif3d_snp_legend, width = snp_legend_width, height = snp_legend_height)


hsp90ab1_snp_legend = generate_legend("Hsp90ab1")
hsp90ab1_snp_legend
save_plot_pdf("hsp90ab1_legend.pdf", hsp90ab1_snp_legend, width = snp_legend_width, height = snp_legend_height)

folr1_snp_legend = generate_legend( "Folr1" )
folr1_snp_legend
save_plot_pdf("folr1_legend.pdf", folr1_snp_legend, width = snp_legend_width, height = snp_legend_height)


supp_legend_genes = c("Cdk1", "Baz1a", "Lclat1", "Umps", "Mrps9", "Nin",
                      "Aff1", "Pttg1")

for(g in supp_legend_genes){
  this_legend_p = generate_legend( g )
  this_file     = paste( g, "legend.pdf", sep = "_"  )
  save_plot_pdf(this_file, this_legend_p, width = snp_legend_width, height = snp_legend_height)
}


################################################################################
###             SNP Legends for the supp figure                             ####


################################################################################
#########      M A I N    S N P   F I G U R E   L A Y O U T         ############

schematic_label = ggdraw() + 
  draw_label(
    "Schematic",
    fontface   = 'plain',
    fontfamily = 'helvetica',
    size       = FONT_LABEL_SIZE,
    angle      = 45
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

main_top = plot_grid(schematic_label,  main_paternal_percentage_figure , ncol = 2, rel_widths = c(2.5, 1) )

main_middle = plot_grid( 
          hsp90ab1_detailed_plot, ncoa3_detailed_plot,
           eif3d_detailed_plot, folr1_detailed_plot,
          ncol = 2)

#main_bottom = main_heatmap_figure[[4]]
main_bottom = riboseq_paternal_ratios_heatmap[[4]]

main_snp_plot = plot_grid(main_top,  main_bottom, main_middle,  ncol = 1, rel_heights = c(1, 1, 4))




main_snp_plot_ribo = plot_grid(main_top,  riboseq_paternal_ratios_heatmap[[4]], main_middle,  ncol = 1, rel_heights = c(1, 1, 4))

save_plot_pdf( "main_snp_plot_ribo.pdf", main_snp_plot_ribo, width = 7.2, height = 9  )

main_snp_plot_rna = plot_grid(main_top,  rnaseq_paternal_ratios_heatmap[[4]], main_middle,  ncol = 1, rel_heights = c(1, 1, 4))

save_plot_pdf( "main_snp_plot_rna.pdf", main_snp_plot_rna, width = 7.2, height = 9  )

################################################################################
#####  P a t e r n a l    A l l e l e     P e r c e n t a g e s

save_plot_pdf("paternal_percentage.pdf", 
              standalone_paternal_percentage_figure, 
              width = 3.5, height = 3)

################################################################################
#####     H e a t m a p s 

standalone_heatmap_figure_no_stars
standalone_heatmap_figure

save_plot_pdf("standalone_heatmap_no_stars.pdf", 
              standalone_heatmap_figure_no_stars, 
              width = 3, height = 3.5 * 2)

save_plot_pdf("standalone_heatmap.pdf", 
              standalone_heatmap_figure, 
              width = 3, height = 3.5 * 2)



save_plot_pdf( "heatmap_highly_translated_genes.pdf", supplementary_heatmap_grid_plot, width = 7.2, height = 9  )

################################################################################


################################################################################
## Supplementary Allele Percentages Figures

save_plot_pdf("supplementary_allele_percentages_ribo.pdf", supplementary_allele_percentages_ribo,
              width = 3.5, height = 4)

save_plot_pdf("supplementary_allele_percentages_rnaseq.pdf", supplementary_allele_percentages_rnaseq,
              width = 3.5, height = 4)

################################################################################

# Higher paternal Ratio in 4-8 cells
nop14_plot   = combine_gene_detail_and_ratio_plots("Nop14") 
tmppe_plot   = combine_gene_detail_and_ratio_plots("Tmppe") 
slc13a2_plot = combine_gene_detail_and_ratio_plots("Slc13a2") 

# Higher paternal Ratio in 8 cells
mrps9_plot = combine_gene_detail_and_ratio_plots("Mrps9") 
tsen2_plot = combine_gene_detail_and_ratio_plots("Tsen2") 
ccnh_plot  = combine_gene_detail_and_ratio_plots("Ccnh") 

detailed_plot_page_1 = 
  plot_grid(nop14_plot, tmppe_plot,
            slc13a2_plot, mrps9_plot,
            tsen2_plot, ccnh_plot,
            ncol = 2, align = "hv")



# Higher paternal Ratio in 4 cells
cdk1_plot   = combine_gene_detail_and_ratio_plots("Cdk1") 
baz1a_plot  = combine_gene_detail_and_ratio_plots("Baz1a") 

# Lower paternal Ratio in 8 cells
dyrk3_plot  = combine_gene_detail_and_ratio_plots("Dyrk3") 
lclat1_plot  = combine_gene_detail_and_ratio_plots("Lclat1") 
lyar_plot  = combine_gene_detail_and_ratio_plots("Lyar") 
umps_plot  = combine_gene_detail_and_ratio_plots("Umps") 

detailed_plot_page_2 = 
  plot_grid( cdk1_plot, baz1a_plot,
             dyrk3_plot, lclat1_plot,
             lyar_plot, umps_plot,
             ncol = 2, align = "hv")

detailed_plot_page_2

# Lower paternal Ratio in 4-8 cells

ppp2ca_plot  = combine_gene_detail_and_ratio_plots("Ppp2ca") 
srpk1_plot    = combine_gene_detail_and_ratio_plots("Srpk1") 
cbx3_plot  = combine_gene_detail_and_ratio_plots("Cbx3") 

# Lower paternal Ratio in 4 cells
pa2g4_plot   = combine_gene_detail_and_ratio_plots("Pa2g4") 
zfp296_plot  = combine_gene_detail_and_ratio_plots("Zfp296") 
nin_plot     = combine_gene_detail_and_ratio_plots("Nin") 

detailed_plot_page_3 = 
  plot_grid( ppp2ca_plot, srpk1_plot ,
             cbx3_plot , pa2g4_plot,
             zfp296_plot, nin_plot,
             ncol = 2, align = "hv")



ddx21_plot     = combine_gene_detail_and_ratio_plots("Ddx21") 
bcat1_plot     = combine_gene_detail_and_ratio_plots("Bcat1") 
mysm1_plot     = combine_gene_detail_and_ratio_plots("Mysm1") 


pttg1_plot     = combine_gene_detail_and_ratio_plots("Pttg1") 
npm1_plot     = combine_gene_detail_and_ratio_plots("Npm1")
aff1_plot     = combine_gene_detail_and_ratio_plots("Aff1")

detailed_plot_page_4 = 
  plot_grid( ddx21_plot, bcat1_plot,
             mysm1_plot, pttg1_plot,
             npm1_plot, aff1_plot,
             ncol = 2, align = "hv")

#detailed_plot_page_4

save_plot_pdf( "detailed_plot_page_1.pdf", detailed_plot_page_1, width = 7.2, height = 8.9  )
save_plot_pdf( "detailed_plot_page_2.pdf", detailed_plot_page_2, width = 7.2, height = 8.9  )
save_plot_pdf( "detailed_plot_page_3.pdf", detailed_plot_page_3, width = 7.2, height = 8.9  )
save_plot_pdf( "detailed_plot_page_4.pdf", detailed_plot_page_4, width = 7.2, height = 8.9  )

################################################################################


combine_gene_detail_and_ratio_plots("Cdt1")

#################################################

# 
# 
# combine_gene_detail_and_ratio_plots("Rpa1")
# 
# combine_gene_detail_and_ratio_plots("Slc6a8")
# 
# 
# ncoa_plot = combine_gene_detail_and_ratio_plots("Ncoa3") 
# nin_plot = combine_gene_detail_and_ratio_plots("Nin") 
# lyar_plot = combine_gene_detail_and_ratio_plots("Lyar")
# top1_plot = combine_gene_detail_and_ratio_plots("Top1")
# tpx2_plot = combine_gene_detail_and_ratio_plots("Tpx2")
# mcm7_plot = combine_gene_detail_and_ratio_plots("Mcm7")
# rpl38_plot = combine_gene_detail_and_ratio_plots("Rpl38")
# kdm5b_plot = combine_gene_detail_and_ratio_plots("Kdm5b")
# ppat_plot = combine_gene_detail_and_ratio_plots("Ppat")
# 

################################################################################
### t-test for the aggregated paternal ratios

aggregated_riboseq_counts =
detailed_riboseq_table %>% 
  filter(group %in% c("4cell", "8cell")) %>%
  group_by(experiment) %>%
  summarise(aggregated_maternal_raw = sum(maternal),
            aggregated_paternal_raw = sum(paternal))

View(aggregated_riboseq_counts)

ribo_2cell_paternal_percentages = 
  corrected_ribo_percentages %>%
  filter( Experiment %in% c("2cell-1", "2cell-2", "2cell-3") )

rna_2cell_paternal_percentages = 
  corrected_rna_percentages %>%
  filter( Experiment %in% c("2cell-1", "2cell-2", "2cell-3", "2cell-4") )

ribo_4cell_paternal_percentages = 
  corrected_ribo_percentages %>%
  filter( Experiment %in% c("4cell-1", "4cell-2", "4cell-3") )

ribo_8cell_paternal_percentages = 
  corrected_ribo_percentages %>%
  filter( Experiment %in% c("8cell-1", "8cell-2", "8cell-3", "8cell-4") )

rna_4cell_paternal_percentages = 
  corrected_rna_percentages %>%
  filter( Experiment %in% c("4cell-1", "4cell-2") )

rna_8cell_paternal_percentages = 
  corrected_rna_percentages %>%
  filter( Experiment %in% c("8cell-1", "8cell-2", "8cell-3", "8cell-4") )


result_t_test_2cell = 
  t.test(ribo_2cell_paternal_percentages$Percentage, rna_2cell_paternal_percentages$Percentage)
print(paste("2cell comparison p-value is", result_t_test_2cell$p.value, sep = " " ))

result_t_test_4cell = 
    t.test(ribo_4cell_paternal_percentages$Percentage, rna_4cell_paternal_percentages$Percentage)
print(paste("4cell comparison p-value is", result_t_test_4cell$p.value, sep = " " ))

result_t_test_8cell = 
  t.test(ribo_8cell_paternal_percentages$Percentage, rna_8cell_paternal_percentages$Percentage)
print(paste("8cell comparison p-value is", result_t_test_8cell$p.value, sep = " " ))


################################################################################
### Sanity Checks on our filtering

# The intersections should produce empty sets in each case.

ribo_t = 
  detailed_riboseq_table %>%
  filter(group %in% c("2cell", "4cell", "8cell")) %>%
  group_by(group, transcript) %>%
  summarise( maternal_total = sum(maternal), paternal_total = sum(paternal) )

rna_t = 
  detailed_rnaseq_table %>%
  filter(group %in% c("2cell", "4cell", "8cell")) %>%
  group_by(group, transcript) %>%
  summarise( maternal_total = sum(maternal), paternal_total = sum(paternal) )

ribo_intr_genes =  ribo_t %>% filter(maternal_total ==2 & paternal_total >= 8)
rna_intr_genes =  rna_t %>% filter(maternal_total ==2 & paternal_total >= 8)

ribo_intr_genes_2cells = (ribo_intr_genes %>% filter(group == "2cell"))$transcript
rna_intr_genes_2cells  = (rna_intr_genes %>% filter(group == "2cell"))$transcript

intersect(ribo_intr_genes_2cells, rna_intr_genes_2cells)

ribo_intr_genes_4cells = (ribo_intr_genes %>% filter(group == "4cell"))$transcript
rna_intr_genes_4cells  = (rna_intr_genes %>% filter(group == "4cell"))$transcript

intersect(ribo_intr_genes_4cells, rna_intr_genes_4cells)

ribo_intr_genes_8cells = (ribo_intr_genes %>% filter(group == "8cell"))$transcript
rna_intr_genes_8cells  = (rna_intr_genes %>% filter(group == "8cell"))$transcript

intersect(ribo_intr_genes_8cells, rna_intr_genes_8cells)