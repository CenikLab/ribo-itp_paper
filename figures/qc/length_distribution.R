#SNP Figure

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
###########      L E N G T H      D I S T R I B U T I O N       ################      

mouse_length_distribution = 
  get_length_distribution(ribo.object = mouse_ribo,
                          region      = "CDS",
                          compact     = FALSE
                          )

# Rename the experiments
mouse_length_distribution$experiment = rename_experiments( mouse_length_distribution$experiment  )

# Add experiment group column "stage"
mouse_length_distribution = 
  mouse_length_distribution %>%
  mutate( stage = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) 
  
# Sum lengths of experiments and find the percentages
mouse_stage_length_dist = 
  mouse_length_distribution %>%
  group_by(experiment) %>%
  mutate( experiment_sum = sum(count) ) %>%
  mutate( experiment_percentage = (count / experiment_sum) * 100 ) 

# Add number of replicates
mouse_stage_length_dist = 
  mouse_stage_length_dist %>%
    group_by(stage) %>%
    mutate(replicate_count = length(unique(experiment)) )


# Add sd after grouping by stage and length
mouse_stage_length_dist =
  mouse_stage_length_dist%>%
  group_by(stage, length) %>%
  summarise( average_percentage = mean(experiment_percentage) , 
             experiment_se_high = average_percentage + 
                                  ( sd(experiment_percentage) / sqrt(replicate_count) ),
             experiment_se_low = average_percentage - 
                                  (sd(experiment_percentage) / sqrt(replicate_count) ) )



plot_stage_length_distribution = function(df, this_stage, color, plot_title, ymax = 10){
  this_plot = 
    ggplot(data=
             df %>% filter(stage == this_stage) , 
           aes(x=length, y=average_percentage )  )  +
      geom_line(aes(color = stage)) + 
      geom_ribbon(aes(ymin = experiment_se_low, ymax = experiment_se_high, fill = stage), alpha = 0.15) + 
      theme(
        panel.border      = element_blank(),
        panel.grid        = element_blank(),
        plot.title        = element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.background  = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.line         = element_line(colour = "black", size = 0.35), 
        legend.position   = "none"
      ) + 
      scale_colour_manual(values = c(color), name = "stage"   ) +
      scale_y_continuous(limits = c(0,ymax), breaks = seq(0,ymax,5), expand = c(0, 0.2)) +
      scale_x_continuous(expand = c(0,0), limits = c(21,40), breaks = c(21,25,30,35,40)) +
      scale_fill_manual(values=c(color) ) + 
      #ylim(c(0, 10)) + 
      labs(title = plot_title, y = "percentage", x= "length")
  
  return(this_plot)
}


plot_1cell = plot_stage_length_distribution(mouse_stage_length_dist, "1cell", BURNT_ORANGE, plot_title = "1cell")
plot_2cell = plot_stage_length_distribution(mouse_stage_length_dist, "2cell", BURNT_ORANGE, plot_title = "2cell")
plot_4cell = plot_stage_length_distribution(mouse_stage_length_dist, "4cell", BURNT_ORANGE, plot_title = "4cell")
plot_8cell = plot_stage_length_distribution(mouse_stage_length_dist, "8cell", BURNT_ORANGE, plot_title = "8cell")
plot_MII   = plot_stage_length_distribution(mouse_stage_length_dist, "MII", BURNT_ORANGE, plot_title   = "MII")
plot_GV    = plot_stage_length_distribution(mouse_stage_length_dist, "GV", BURNT_ORANGE, plot_title = "GV")


# title_main = ggdraw() + 
#   draw_label(
#     "Length distribution of ribosoe protected footprints",
#     fontface   = 'plain',
#     size       = FONT_TITLE_SIZE ,
#     fontfamily = 'helvetica',
#   ) +
#   theme(
#     # add margin on the left of the drawing canvas,
#     # so title is aligned with left edge of first plot
#     plot.margin = margin(0, 0, 0, 0)
#   )

raw_plot = plot_grid(plot_GV, plot_MII,
                     plot_1cell, plot_2cell, 
                     plot_4cell, plot_8cell,
                     ncol = 2)


# main_plot = plot_grid(title_main, raw_plot,
#                       ncol = 1,
#                       rel_heights = c(0.1, 1))

main_plot_mouse = raw_plot

main_plot_mouse


mouse_gv_mii = plot_grid( plot_GV, plot_MII, nrow = 1)
mouse_gv_mii

#######################################################################################

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


save_plot_pdf("length_distribution_mouse_main.pdf", main_plot_mouse, width = 3.54, height = 3.54 * 1.5)
#save_plot_pdf("length_distribution_mouse_mii_gv.pdf", mouse_mii_gv, width = 3.54, height = 3.54/2)

plot_length = 2.25

save_plot_pdf("length_distribution_1cell.pdf", plot_1cell, width = plot_length, height = plot_length)
save_plot_pdf("length_distribution_2cell.pdf", plot_2cell, width = plot_length, height = plot_length)
save_plot_pdf("length_distribution_4cell.pdf", plot_4cell, width = plot_length, height = plot_length)
save_plot_pdf("length_distribution_8cell.pdf", plot_8cell, width = plot_length, height = plot_length)
save_plot_pdf("length_distribution_MII.pdf", plot_MII, width = plot_length, height = plot_length)
save_plot_pdf("length_distribution_GV.pdf", plot_GV, width = plot_length, height = plot_length)

################################################################################

human_length_distribution = 
  get_length_distribution(ribo.object = human_ribo,
                          region      = "CDS",
                          compact     = FALSE
  )

# Rename the experiments
human_length_distribution$experiment = human_rename_experiments( human_length_distribution$experiment  )

# Add experiment group column "stage"
human_length_distribution = 
  human_length_distribution %>%
  mutate( stage = unlist(lapply ( strsplit(  as.vector(experiment), split = "-" ), "[[", 1) )  ) 


# Sum lengths of experiments and find the percentages
human_stage_length_dist = 
  human_length_distribution %>%
  group_by(experiment) %>%
  mutate( experiment_sum = sum(count) ) %>%
  mutate( experiment_percentage = (count / experiment_sum) * 100 ) 

## Add replicate number
human_stage_length_dist = 
  human_stage_length_dist %>%
  group_by(stage) %>%
    mutate(replicate_count = length(unique(experiment)) )

human_stage_length_dist =
  human_stage_length_dist%>%
  group_by(stage, length) %>%
  summarise( average_percentage = mean(experiment_percentage) , 
             experiment_se_high = average_percentage + 
                                  (  sd(experiment_percentage) / sqrt(replicate_count)  ) ,
             experiment_se_low = average_percentage - 
                                 ( sd(experiment_percentage) /sqrt(replicate_count)  )  )


plot_100 = plot_stage_length_distribution(human_stage_length_dist, "100", BURNT_ORANGE, ymax = 12, plot_title = "100 cells" ) +
              scale_y_continuous(breaks = c(0, 6 ,12), limits = c(0,12), expand = c(0,0.2) )
plot_10M = plot_stage_length_distribution(human_stage_length_dist, "10M", BURNT_ORANGE, ymax = 12, plot_title = "10M cells"  ) + 
             scale_y_continuous(breaks = c(0, 6 ,12), limits = c(0,12), expand = c(0,0.2) )


save_plot_pdf("length_distribution_100.pdf", plot_100, width = plot_length, height = plot_length)
save_plot_pdf("length_distribution_10M.pdf", plot_10M, width = plot_length, height = plot_length)

human_length_distribution_combined = 
  plot_grid(  plot_100, plot_10M, nrow = 1 )

human_length_distribution_combined
save_plot_pdf("length_distribution_human_combined.pdf", 
              human_length_distribution_combined, 
              width = 3.54, height = 3.54 / 2)
