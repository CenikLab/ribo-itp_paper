#REGION COUNTS

mouse_ribo_file = "../../../mouse-itp_v5.ribo"
human_ribo_file = "../../../../itp/human-itp_v4.ribo"


################################################################################
########                    L I B R A R I E S                          #########

library(ribor)

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

mouse_ribo = Ribo(mouse_ribo_file, rename = rename_default)
human_ribo = Ribo(human_ribo_file, rename = rename_default)

################################################################################

MOUSE_MIN_LENGTH = 29
MOUSE_MAX_LENGTH = 35

################################################################################

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

#CDS_GREEN  = "#00BA38"
#UTR5_BLUE  = "#619CFF"
#UTR3_GREEN = "619CFF"

CDS_GREEN  = rgb(124, 203, 162, maxColorValue = 255)
UTR5_BLUE  = rgb(4,   82,  117, maxColorValue = 255)
UTR3_GREEN = rgb(240, 116, 110,maxColorValue = 255)

################################################################################
#########                 F O N T   S I Z E S                          #########

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 9

PDF_resolution = 600
FIGURE_FONT = "helvetica"


################################################################################
#######                   E X P E R I M E N T S                      ###########


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




mouse_exp_name_mapper = make_name_mapper( mouse_ribo@experiments  )


experiment_name_mapper = mouse_exp_name_mapper

rename_experiments = function(experiment_name){
  return(experiment_name_mapper[[experiment_name]])
}

rename_experiments = Vectorize(rename_experiments, USE.NAMES = FALSE)


human_exp_name_mapper = make_name_mapper( human_ribo@experiments  )

human_rename_experiments = function(experiment_name){
  return(human_exp_name_mapper[[experiment_name]])
}

human_rename_experiments = Vectorize(human_rename_experiments, USE.NAMES = FALSE)
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
###########               R E G I O N     C O U N T S           ################

mouse_region_counts = get_region_counts(ribo.object = mouse_ribo,
                                        region      = c("UTR5", "CDS", "UTR3"),
                                        range.lower = MOUSE_MIN_LENGTH,
                                        range.upper = MOUSE_MAX_LENGTH,
                                        compact     = FALSE)

human_region_counts = get_region_counts(ribo.object = human_ribo,
                                        region      = c("UTR5", "CDS", "UTR3"),
                                        range.lower = MOUSE_MIN_LENGTH,
                                        range.upper = MOUSE_MAX_LENGTH,
                                        compact     = FALSE)

mouse_region_counts$experiment = rename_experiments(mouse_region_counts$experiment)

human_region_counts$experiment = human_rename_experiments(human_region_counts$experiment)

plot_region_counts(x           = human_ribo,
                   range.lower = MOUSE_MIN_LENGTH,
                   range.upper = MOUSE_MAX_LENGTH)


plot_region_counts(x           = mouse_ribo,
                   range.lower = MOUSE_MIN_LENGTH,
                   range.upper = MOUSE_MAX_LENGTH)

###############################################################################

compute_region_percentages = function(df){
  this_df =
    df %>% 
    group_by(experiment) %>%
    mutate(experiment_total = sum(count))
  
  this_df =
    this_df %>%
    mutate( percentage = round(100 * (count / experiment_total), 1 )  ) %>%
    mutate(region = factor(region, levels = c("UTR3", "CDS", "UTR5") ))
    
  
  
  return(this_df)
}

################################################################################

customized_plot_region_counts = function(df, plot_title = "Distribution of RPFs by transcript regions"){
  
  percentages <- replace(df$percentage, 
                         df$region != "CDS", "")
  
  this_df = df
  
  this_plot = 
    ggplot(this_df, 
           aes(x=experiment, y=percentage, fill=region)) +
    geom_col() +
    coord_flip() + 
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
          panel.border = element_blank(),
          panel.grid   = element_blank(),
          axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          legend.title      = element_blank(),
          legend.key.size   = unit(0.15, 'inches')) +
    scale_fill_manual("legend", 
                      values = c("UTR5" = UTR5_BLUE, "UTR3" = UTR3_GREEN, "CDS" = CDS_GREEN), 
                      breaks = c("UTR5", "CDS", "UTR3"), 
                      labels = c("5' UTR", "CDS", "3' UTR")) +
    geom_text(aes(x=.data$experiment, y=50, label=percentages), size=3) +
    scale_x_discrete(limits = rev) + 
    labs(title = plot_title, 
         fill  = "Region", 
         x     = "Experiment", 
         y     = "Percentage")
  
  return(this_plot)
}

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
#####           S U P P L E M E N T A R Y     F I G U R E S             ########    

mouse_region_percentages       = compute_region_percentages(mouse_region_counts)

mouse_region_percentages = 
  mouse_region_percentages %>%
  mutate( group = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) %>%
  group_by(group) %>%
  mutate( replicate_count = length( unique( experiment )  ) ) 

mouse_region_percentages$group = 
  factor(mouse_region_percentages$group , 
         levels = c("GV", "MII", "1cell", "2cell", "4cell", "8cell"))

mouse_region_percentages$experiment = 
  factor( mouse_region_percentages$experiment, 
          levels =  
            c( unique( (mouse_region_percentages %>% filter( group == "GV" ) )$experiment ),
               unique( (mouse_region_percentages %>% filter( group == "MII" ) )$experiment ),
               unique( (mouse_region_percentages %>% filter( group == "1cell" ) )$experiment ),
               unique( (mouse_region_percentages %>% filter( group == "2cell" ) )$experiment ),
               unique( (mouse_region_percentages %>% filter( group == "4cell" ) )$experiment ),
               unique( (mouse_region_percentages %>% filter( group == "8cell" ) )$experiment )) )

human_region_percentages = compute_region_percentages(human_region_counts)

mouse_supplementary_plot = 
  customized_plot_region_counts(mouse_region_percentages)

human_supplementary_plot = 
  customized_plot_region_counts(human_region_percentages)

save_plot_pdf("mouse_region_counts_supp.pdf", mouse_supplementary_plot, width = 3.54, height = 3.54*1.5)

save_plot_pdf("human_region_counts_supp.pdf", human_supplementary_plot, width = 3.54, height = 1.9)

################################################################################


################################################################################
########           A G G R E G A T E D     F I G U R E S             ########### 



 

mouse_region_percentages =
  mouse_region_percentages %>%
  group_by(group, region) %>%
  mutate(average_percentage = mean(percentage) ) %>%
  mutate(standard_error = sd(percentage) / sqrt(replicate_count) )
  

human_region_percentages = 
  human_region_percentages %>%
  mutate( group = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) %>%
  group_by(group) %>%
  mutate( replicate_count = length( unique( experiment )  ) ) 

human_region_percentages =
  human_region_percentages %>% 
  group_by(group, region) %>%
  mutate(average_percentage = mean(percentage) ) %>%
  mutate(standard_error = sd(percentage) / sqrt(replicate_count) )

plot_bar_with_error_bars = function(df, 
                                    plot_title = "Distribution of RPFs by transcript regions"){
  this_df = df
  this_df$region = factor(this_df$region, levels = c("UTR5", "CDS", "UTR3"))
  
  this_plot = 
    ggplot(data=this_df, aes(x= group, y=average_percentage, fill = region )  )  + 
      geom_bar(position = "dodge", stat="identity", alpha = 1 ) + 
      geom_errorbar( aes(  x        = group, 
                           ymin     = average_percentage - standard_error, 
                           ymax     = average_percentage + standard_error), 
                           width    = 0.4, 
                           alpha    = 0.4, 
                           size     = 0.4,
                           position = position_dodge(width = 0.9)) +
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
            legend.key.size  = unit(0.12, 'inches')) + 
      labs(title= plot_title, fill="Region", x="Stage", y="Average percentage") + 
      scale_y_continuous(limits = c(0, 100), expand = c(0,0)) + 
      scale_fill_manual("legend", 
                        values = c("UTR5" = UTR5_BLUE, "UTR3" = UTR3_GREEN, "CDS" = CDS_GREEN), 
                        breaks = c("UTR5", "CDS", "UTR3"), 
                        labels = c("5' UTR", "CDS", "3' UTR")) 
  
  return(this_plot)
}

mouse_region_counts_with_error_bars = 
    plot_bar_with_error_bars(mouse_region_percentages)

save_plot_pdf("mouse_region_counts_with_error_bars.pdf", 
              mouse_region_counts_with_error_bars, 
              width = 3.54, height = 2.6)

human_region_counts_with_error_bars = 
    plot_bar_with_error_bars(human_region_percentages, plot_title = "")

save_plot_pdf("human_region_counts_with_error_bars.pdf", 
              human_region_counts_with_error_bars, 
              width = 2.1, height = 2.6)

################################################################################

# Distribution of translated gene region lengths
################################################################################
### M O U S E 

mouse_region_counts_per_gene = 
  get_region_counts(ribo.object = mouse_ribo,
                    region      = c("UTR5", "CDS", "UTR3"),
                    range.lower = MOUSE_MIN_LENGTH,
                    range.upper = MOUSE_MAX_LENGTH,
                    transcript  = FALSE,
                    compact     = FALSE)

mouse_region_total_counts = 
  mouse_region_counts_per_gene %>%
  group_by(experiment, transcript) %>%
  summarise(total_count = sum(count))

mouse_region_lengths = get_region_lengths(mouse_ribo)

mouse_region_lengths = mouse_region_lengths %>% select(transcript, UTR5, CDS, UTR3)



mouse_region_total_counts_and_lengths = merge(mouse_region_total_counts, mouse_region_lengths)



mouse_region_total_counts_and_lengths = 
  mouse_region_total_counts_and_lengths %>% 
    mutate(total_length = sum(UTR5 + CDS + UTR3)) %>%
    mutate( UTR5_w_ratio = (UTR5 / total_length) * total_count,
            CDS_w_ratio  = (CDS / total_length) * total_count,
            UTR3_w_ratio = (UTR3 / total_length) * total_count )


mouse_region_total_counts_and_lengths$experiment = rename_experiments(mouse_region_total_counts_and_lengths$experiment)

mouse_region_total_counts_and_lengths = 
  mouse_region_total_counts_and_lengths %>%
  mutate( group = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) %>%
  group_by(group) %>%
  mutate( replicate_count = length( unique( experiment )  ) ) 




  
mouse_region_len_percentages_pre = 
  mouse_region_total_counts_and_lengths %>%
  group_by(experiment, group) %>%
  mutate(UTR5_w_ratio_total = sum(UTR5_w_ratio),
         CDS_w_ratio_total  = sum(CDS_w_ratio),
         UTR3_w_ratio_total = sum(UTR3_w_ratio),
         all_total = UTR5_w_ratio_total + CDS_w_ratio_total + UTR3_w_ratio_total) %>%
  summarise( UTR5_percentage =  round( median((UTR5_w_ratio_total / all_total) * 100) , 1) ,
             CDS_percentage  =  round( median((CDS_w_ratio_total / all_total) * 100) , 1),
             UTR3_percentage =  round( median((UTR3_w_ratio_total / all_total) * 100 ), 1)  )%>%
  rename_at( "UTR5_percentage", ~"UTR5" ) %>%
  rename_at( "CDS_percentage", ~"CDS" ) %>%
  rename_at( "UTR3_percentage", ~"UTR3" )


mouse_region_len_percentages  = 
  melt(mouse_region_len_percentages_pre, value.name = "percentage") %>% 
  rename_at("variable", ~"region") %>%
  mutate(region = factor(region, levels = c("UTR3", "CDS", "UTR5") ))


mouse_region_len_percentages$group = 
  factor(mouse_region_len_percentages$group , 
         levels = c("GV", "MII", "1cell", "2cell", "4cell", "8cell"))

mouse_region_len_percentages$experiment = 
  factor( mouse_region_len_percentages$experiment, 
          levels =  
            c( unique( (mouse_region_len_percentages %>% filter( group == "GV" ) )$experiment ),
               unique( (mouse_region_len_percentages %>% filter( group == "MII" ) )$experiment ),
               unique( (mouse_region_len_percentages %>% filter( group == "1cell" ) )$experiment ),
               unique( (mouse_region_len_percentages %>% filter( group == "2cell" ) )$experiment ),
               unique( (mouse_region_len_percentages %>% filter( group == "4cell" ) )$experiment ),
               unique( (mouse_region_len_percentages %>% filter( group == "8cell" ) )$experiment )) )



mouse_region_len_percentage_barplot = 
  customized_plot_region_counts(mouse_region_len_percentages) + 
  ggtitle("Distribution of Region Lengths")

###

mouse_region_len_percentages = 
  mouse_region_len_percentages %>%
  group_by(group) %>%
  mutate( replicate_count = length( unique( experiment )  ) ) 

mouse_region_len_percentages_with_error = 
  mouse_region_len_percentages %>%
    group_by(group, region) %>%
    mutate(average_percentage = mean(percentage) ) %>%
    mutate(standard_error = sd(percentage) / sqrt(replicate_count) )

mouse_region_len_percentages_with_error_bars = 
  plot_bar_with_error_bars(mouse_region_len_percentages_with_error,
                           "Region lengths weighted by ribosome occupancy")

# EXPORT THIS FOR THE MAIN FIGURE

legend_of_mouse_region_counts_comperative_plot = get_legend(mouse_region_counts_with_error_bars)

mouse_region_counts_comperative_plot = 
  plot_grid( mouse_region_counts_with_error_bars + theme(legend.position = "none"),
             mouse_region_len_percentages_with_error_bars + theme(legend.position = "none"),
             legend_of_mouse_region_counts_comperative_plot,
             ncol = 3,
             rel_widths = c(1,1,0.2))



################################################################################
### H U M A N





human_region_counts_per_gene = 
  get_region_counts(ribo.object = human_ribo,
                    region      = c("UTR5", "CDS", "UTR3"),
                    range.lower = MOUSE_MIN_LENGTH,
                    range.upper = MOUSE_MAX_LENGTH,
                    transcript  = FALSE,
                    compact     = FALSE)


################################################################################
####             Chi-squared Test for Region Counts (HUMAN)                 ####

human_counts_wide = dcast( human_region_counts_per_gene,  experiment + transcript ~ region)
colnames(human_counts_wide) = c(colnames(human_counts_wide)[1:2], "CDS_count", "UTR3_count", "UTR5_count")

human_region_total_counts = 
  human_region_counts_per_gene %>%
  group_by(experiment, transcript) %>%
  summarise(total_count = sum(count))

human_region_lengths = get_region_lengths(human_ribo)

human_region_lengths = human_region_lengths %>% select(transcript, UTR5, CDS, UTR3)

human_counts_and_lengths_wide = merge( human_counts_wide, human_region_lengths  )

human_region_total_counts_and_lengths = merge(human_region_total_counts, human_region_lengths)

human_expected_counts = human_counts_and_lengths_wide %>%
  group_by(experiment, transcript) %>%
  mutate(total_count = sum(CDS_count, UTR3_count, UTR5_count), total_length = sum(CDS, UTR3, UTR5)) %>%
  mutate(CDS_expected = total_count * (CDS / total_length),
         UTR5_expected = total_count * (UTR5 / total_length),
         UTR3_expected = total_count * (UTR3/ total_length))

human_expected_totals = 
  human_expected_counts %>% 
  group_by(experiment) %>%
  summarise( CDS_total           = sum(CDS_count),
             UTR3_total          = sum(UTR3_count),
             UTR5_total          = sum(UTR5_count),
             CDS_expected_total  = sum(CDS_expected),
             UTR5_expected_total = sum(UTR5_expected),
             UTR3_expected_total = sum(UTR3_expected),
             ) 

human_chisquared_res = human_expected_totals %>%
  group_by(experiment) %>%
  summarise( pval = chisq.test( x = c(CDS_total, UTR3_total, UTR5_total),
                                p = c(CDS_expected_total, UTR3_expected_total, UTR5_expected_total),
                                rescale.p = T)$p.value  )


rep1_chi = 
chisq.test( x = unlist(human_expected_totals[4,2:4]),
            p = unlist(human_expected_totals[4,5:7]),
            rescale.p = T)

print(human_expected_totals[4,1])
rep1_chi

rep2_chi = 
  chisq.test( x = unlist(human_expected_totals[5,2:4]),
              p = unlist(human_expected_totals[5,5:7]),
              rescale.p = T)

print(human_expected_totals[5,1])
rep2_chi

rep3_chi = 
  chisq.test( x = unlist(human_expected_totals[6,2:4]),
              p = unlist(human_expected_totals[6,5:7]),
              rescale.p = T)

print(human_expected_totals[6,1])
rep3_chi


################################################################################
####             Chi-squared Test for Region Counts (MOUSE)                 ####

mouse_counts_wide           = dcast( mouse_region_counts_per_gene,  experiment + transcript ~ region)
colnames(mouse_counts_wide) = c(colnames(mouse_counts_wide)[1:2], "CDS_count", "UTR3_count", "UTR5_count")

mouse_region_total_counts = 
  mouse_region_counts_per_gene %>%
  group_by(experiment, transcript) %>%
  summarise(total_count = sum(count))

mouse_region_lengths = get_region_lengths(mouse_ribo)

mouse_region_lengths = mouse_region_lengths %>% select(transcript, UTR5, CDS, UTR3)

mouse_counts_and_lengths_wide = merge( mouse_counts_wide, mouse_region_lengths  )

mouse_region_total_counts_and_lengths = merge(mouse_region_total_counts, mouse_region_lengths)

mouse_expected_counts = mouse_counts_and_lengths_wide %>%
  group_by(experiment, transcript) %>%
  mutate(total_count = sum(CDS_count, UTR3_count, UTR5_count), total_length = sum(CDS, UTR3, UTR5)) %>%
  mutate(CDS_expected = total_count * (CDS / total_length),
         UTR5_expected = total_count * (UTR5 / total_length),
         UTR3_expected = total_count * (UTR3/ total_length))


mouse_expected_totals = 
  mouse_expected_counts %>% 
  group_by(experiment) %>%
  summarise( CDS_total           = sum(CDS_count),
             UTR3_total          = sum(UTR3_count),
             UTR5_total          = sum(UTR5_count),
             CDS_expected_total  = sum(CDS_expected),
             UTR5_expected_total = sum(UTR5_expected),
             UTR3_expected_total = sum(UTR3_expected),
  ) 



################################################################################

human_region_total_counts_and_lengths$experiment = human_rename_experiments(human_region_total_counts_and_lengths$experiment)

human_region_total_counts_and_lengths = 
  human_region_total_counts_and_lengths %>%
  mutate( group = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) %>%
  group_by(group) %>%
  mutate( replicate_count = length( unique( experiment )  ) ) 





human_region_len_percentages_pre = 
  human_region_total_counts_and_lengths %>%
  group_by(experiment, group) %>%
  mutate(UTR5_w_ratio_total = sum(UTR5_w_ratio),
         CDS_w_ratio_total  = sum(CDS_w_ratio),
         UTR3_w_ratio_total = sum(UTR3_w_ratio),
         all_total = UTR5_w_ratio_total + CDS_w_ratio_total + UTR3_w_ratio_total) %>%
  summarise( UTR5_percentage =  round( median((UTR5_w_ratio_total / all_total) * 100) , 1) ,
             CDS_percentage  =  round( median((CDS_w_ratio_total / all_total) * 100) , 1),
             UTR3_percentage =  round( median((UTR3_w_ratio_total / all_total) * 100 ), 1)  )%>%
  rename_at( "UTR5_percentage", ~"UTR5" ) %>%
  rename_at( "CDS_percentage", ~"CDS" ) %>%
  rename_at( "UTR3_percentage", ~"UTR3" )


human_region_len_percentages  = 
  melt(human_region_len_percentages_pre, value.name = "percentage") %>% 
  rename_at("variable", ~"region") %>%
  mutate(region = factor(region, levels = c("UTR3", "CDS", "UTR5") ))

# Our ultimate plot for mouse region lengths
human_region_len_percentage_barplot = 
  customized_plot_region_counts(human_region_len_percentages) + 
  ggtitle("Distribution of Region Lengths")


save_plot_pdf("mouse_region_len_percentage_barplot_supp.pdf", mouse_region_len_percentage_barplot, width = 3.54, height = 3.54*1.5)

save_plot_pdf("human_region_len_percentage_barplot_supp.pdf", human_region_len_percentage_barplot, width = 3.54, height = 1.9)
