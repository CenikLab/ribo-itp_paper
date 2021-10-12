#REGION COUNTS

mouse_ribo_file = "../../../mouse-itp_v5.ribo"
human_ribo_file = "../../../../itp/human-itp_v4.ribo"

GSE101018_file = "../../../external_data/GSE101018.ribo"
GSE78634_file = "../../../external_data/GSE78634.ribo"

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

GSE101018_ribo = Ribo(GSE101018_file, rename = rename_default)
GSE78634_ribo  = Ribo(GSE78634_file,  rename = rename_default)

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

CDS_GREEN  = "#00BA38"
UTR5_BLUE  = "#619CFF"
UTR3_GREEN = "619CFF"

START_SITE_COLOR = rgb(124,203,162, max = 255)
STOP_SITE_COLOR  = rgb(8,144,153, max = 255)



################################################################################
#########                 F O N T   S I Z E S                          #########

FONT_LABEL_SIZE = 8
FONT_TITLE_SIZE = 9

PDF_resolution = 600
FIGURE_FONT = "helvetica"

################################################################################
#########                 OTHER FIGURE SETTINGS                         ########

LINE_THICKNESS = 0.5

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


human_metagene_start = 
get_metagene(human_ribo, 
             site = "start",
             range.lower = MOUSE_MIN_LENGTH,
             range.upper = MOUSE_MAX_LENGTH,
             length = FALSE,
             transcript = TRUE,
             alias = TRUE,
             compact = FALSE,
             experiment = human_ribo@experiments)

human_metagene_stop = 
  get_metagene(human_ribo, 
               site = "stop",
               range.lower = MOUSE_MIN_LENGTH,
               range.upper = MOUSE_MAX_LENGTH,
               length = FALSE,
               transcript = TRUE,
               alias = TRUE,
               compact = FALSE,
               experiment = human_ribo@experiments)

human_metagene_start$experiment = human_rename_experiments(human_metagene_start$experiment)
human_metagene_stop$experiment = human_rename_experiments(human_metagene_stop$experiment)

mouse_metagene_start = 
  get_metagene(mouse_ribo, 
               site = "start",
               range.lower = MOUSE_MIN_LENGTH,
               range.upper = MOUSE_MAX_LENGTH,
               length = FALSE,
               transcript = TRUE,
               alias = TRUE,
               compact = FALSE,
               experiment = mouse_ribo@experiments)

mouse_metagene_stop = 
  get_metagene(mouse_ribo, 
               site = "stop",
               range.lower = MOUSE_MIN_LENGTH,
               range.upper = MOUSE_MAX_LENGTH,
               length = FALSE,
               transcript = TRUE,
               alias = TRUE,
               compact = FALSE,
               experiment = mouse_ribo@experiments)


mouse_metagene_start$experiment = rename_experiments( mouse_metagene_start$experiment  )
mouse_metagene_stop$experiment = rename_experiments( mouse_metagene_stop$experiment  )

get_experiment_coverage_at_given_len = function(df, this_exp, len){
  result_df = 
    df %>%
      filter( experiment == this_exp & length == len )
  
  return(result_df)
}



plot_individual_len = function(df, experiment){
  offsets = c()
  
  for( i in seq(MOUSE_MIN_LENGTH, MOUSE_MAX_LENGTH)  ){
    coverage_raw = get_experiment_coverage_at_given_len(df, experiment, i)
    coverage     = t( coverage_raw[-c(1,2)] )
    c            = data.frame("position" = as.integer(row.names(coverage) ), "count" = as.vector(coverage[,1] ))
    
    this_plot = ggplot(data = c, aes(x = position ,
                                     y = count )) + 
                geom_line(  ) + 
                labs(title = i)
    
    #print(i)
    print(which.max(c$count) )
    
    print(this_plot)
  }
  return(offsets)
}

################################################################################

plot_individual_len(mouse_metagene_stop, "8cell-1")
plot_individual_len(mouse_metagene_stop, "8cell-2")

plot_individual_len(mouse_metagene_stop, "2cell-1")
plot_individual_len(mouse_metagene_stop, "2cell-2")
plot_individual_len(mouse_metagene_stop, "2cell-3")

plot_individual_len(mouse_metagene_stop, "4cell-1")
plot_individual_len(mouse_metagene_stop, "4cell-2")
plot_individual_len(mouse_metagene_stop, "4cell-3")

plot_individual_len(mouse_metagene_stop, "1cell-1")
plot_individual_len(mouse_metagene_stop, "1cell-2")
plot_individual_len(mouse_metagene_stop, "1cell-3")
plot_individual_len(mouse_metagene_stop, "1cell-4")

plot_individual_len(mouse_metagene_start, "1cell-1")
plot_individual_len(mouse_metagene_start, "1cell-2")
plot_individual_len(mouse_metagene_start, "1cell-3")
plot_individual_len(mouse_metagene_start, "1cell-4")

plot_individual_len(mouse_metagene_start, "8cell-1")

plot_individual_len(mouse_metagene_start, "8cell-2")


plot_individual_len(mouse_metagene_start, "8cell-3")
plot_individual_len(mouse_metagene_start, "8cell-4")

plot_individual_len(mouse_metagene_start, "4cell-1")
plot_individual_len(mouse_metagene_start, "4cell-2")
plot_individual_len(mouse_metagene_start, "4cell-3")

plot_individual_len(mouse_metagene_start, "2cell-1")
plot_individual_len(mouse_metagene_start, "2cell-2")

################################################################################

plot_individual_len(human_metagene_start, "100-1")
plot_individual_len(human_metagene_start, "100-2")
plot_individual_len(human_metagene_start, "100-3")

plot_individual_len(human_metagene_start, "10M-1")
plot_individual_len(human_metagene_start, "10M-2")
plot_individual_len(human_metagene_start, "10M-3")


plot_individual_len(human_metagene_stop, "100-1")
plot_individual_len(human_metagene_stop, "100-2")
plot_individual_len(human_metagene_stop, "100-3")

plot_individual_len(human_metagene_stop, "10M-1")
plot_individual_len(human_metagene_stop, "10M-2")
plot_individual_len(human_metagene_stop, "10M-3")

z= list()

z[["3"]] = 7


plot_metagene(mouse_ribo,
              site        = "stop",
              experiment  = c("20210513-ITP-1cell-cross-50-A"),
              range.lower =MOUSE_MIN_LENGTH,
              range.upper = MOUSE_MAX_LENGTH)


plot_metagene(mouse_ribo,
              site        = "start",
              experiment  = c("20210513-ITP-1cell-cross-50-C"),
              range.lower =MOUSE_MIN_LENGTH,
              range.upper = MOUSE_MAX_LENGTH)

plot_metagene(mouse_ribo,
              site        = "stop",
              experiment  = c("20210513-ITP-1cell-cross-50-C"),
              range.lower =MOUSE_MIN_LENGTH,
              range.upper = MOUSE_MAX_LENGTH)


plot_metagene(mouse_ribo,
              site        = "stop",
              experiment  = c("20210513-ITP-1cell-cross-50-A"),
              range.lower =25,
              range.upper = 40)

################################################################################

### Take 3' end of the reads on metagene

push_vector_to_right = function(coverage, units){
  main_vector = coverage[ 1:(length(coverage) - units) ]
  left_pad = rep(0, units)
  return( c(left_pad, main_vector) )
}


get_metagene_three_p_ends = function(this_df, experiment){
  
  metagene_coverage = as.integer( 
                          get_experiment_coverage_at_given_len(df = this_df, this_exp = experiment, len = MOUSE_MIN_LENGTH)[c(-1, -2)] )
  metagene_coverage = push_vector_to_right(metagene_coverage, MOUSE_MIN_LENGTH - 1)
  #return(metagene_coverage)
  
  for(this_length in  ( (MOUSE_MIN_LENGTH + 1): MOUSE_MAX_LENGTH) ){
    this_coverage     = as.integer(get_experiment_coverage_at_given_len(df = this_df, this_exp = experiment, len = this_length)[c(-1, -2)] )
    this_coverage = push_vector_to_right(this_coverage, this_length - 1)
    metagene_coverage = metagene_coverage + this_coverage
                        
  }
  
  return(metagene_coverage)
}


get_metagene_five_p_ends = function(this_df, experiment){
  
  metagene_coverage = as.integer( 
    get_experiment_coverage_at_given_len(df = this_df, this_exp = experiment, len = MOUSE_MIN_LENGTH)[c(-1, -2)] )

  #return(metagene_coverage)
  
  for(this_length in  ( (MOUSE_MIN_LENGTH + 1): MOUSE_MAX_LENGTH) ){
    this_coverage     = as.integer(get_experiment_coverage_at_given_len(df = this_df, this_exp = experiment, len = this_length)[c(-1, -2)] )
    metagene_coverage = metagene_coverage + this_coverage
    
  }
  
  return(metagene_coverage)
}



plot_metagene_three_p_ends = function(this_df, experiment){
  metagene_coverage = get_metagene_three_p_ends(this_df, experiment)
  
  this_df = data.frame(position = -50:50, count = metagene_coverage)
  
  this_plot = ggplot(data = this_df, aes(x = position ,
                                   y = count )) + 
    geom_line(  ) + 
    labs(title = experiment)
  
  return(this_plot)
}

plot_metagene_five_p_ends = function(this_df, experiment){
  metagene_coverage = get_metagene_five_p_ends(this_df, experiment)
  
  this_df = data.frame(position = -50:50, count = metagene_coverage)
  
  this_plot = ggplot(data = this_df, aes(x = position ,
                                         y = count )) + 
    geom_line(  ) + 
    labs(title = experiment)
  
  return(this_plot)
}

combine_three_and_five_p = function(this_df, experiment){
  three_p_plot  = plot_metagene_three_p_ends(  this_df, experiment   )
  five_p_plot   = plot_metagene_five_p_ends(  this_df, experiment   )
  
  combined_plot = plot_grid(three_p_plot, five_p_plot, ncol = 1)
  
  return(combined_plot)
}


two_cell_plot = combine_three_and_five_p(  mouse_metagene_stop, "1cell-3"   ) 


eight_cell_plot = combine_three_and_five_p(  mouse_metagene_stop, "8cell-2"   ) 

#ggsave("~/scratch/2cell_comparison.pdf", two_cell_plot)
#ggsave("~/scratch/8cell_comparison.pdf", eight_cell_plot)


combine_three_and_five_p(  human_metagene_stop, "100-1"   ) 

combine_three_and_five_p(  mouse_metagene_stop, "1cell-1"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "1cell-2"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "1cell-3"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "1cell-4"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "MII-1"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "MII-2"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "MII-3"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "MII-4"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "MII-5"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "GV-1"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "GV-2"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "GV-3"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "GV-4"   ) 
combine_three_and_five_p(  mouse_metagene_stop, "GV-5"   ) 



################################################################################
###########  P L O T S   F O R   T H E    M A N U S C R I P T    ###############

# As representative replicates for the metagene plots, we picked thefollowing
# 1cell-5 for mouse
#  100-1 

combine_three_and_five_p(  mouse_metagene_stop, "1cell-3"   ) 
combine_three_and_five_p(  mouse_metagene_start, "1cell-3"   ) 

combine_three_and_five_p(  mouse_metagene_stop, "8cell-1"   ) 
combine_three_and_five_p(  mouse_metagene_start, "8cell-1"   ) 

combine_three_and_five_p(  human_metagene_stop, "100-1"   ) 
combine_three_and_five_p(  human_metagene_start, "100-1"   ) 


combine_three_and_five_p(  human_metagene_stop, "10M-3"   ) 
combine_three_and_five_p(  human_metagene_start, "10M-3"   ) 

# After plotting the stop site for each read length
# We see that the offset is 14 nucleotidea fo all read lengths
plot_individual_len(mouse_metagene_stop, "1cell-2")
plot_individual_len(human_metagene_stop, "100-1")

shift_by_offsets = function(this_df, experiment, offsets ){
  result = get_experiment_coverage_at_given_len(df = this_df, this_exp = experiment, len = MOUSE_MIN_LENGTH)[c(-1, -2)] 
  result = as.integer( push_vector_to_right(result, offsets[1]) )
  
  for(i in 2:(MOUSE_MAX_LENGTH - MOUSE_MIN_LENGTH + 1)  ){
       this_coverage = get_experiment_coverage_at_given_len(df = this_df, this_exp = experiment, len = MOUSE_MIN_LENGTH + i - 1)[c(-1, -2)]
       this_coverage = as.integer(push_vector_to_right(this_coverage, offsets[i]) )
       result = result + this_coverage
  }
  
  return(result)
}

mouse_1cell_5_adjusted_stop_site_coverage = shift_by_offsets(mouse_metagene_stop, "1cell-3", c(15, rep(15, 6 )) )

this_coverage = mouse_1cell_5_adjusted_stop_site_coverage

this_df = data.frame(position = -50:50, count = this_coverage)

this_plot = ggplot(data = this_df, aes(x = position ,
                                       y = count )) + 
  geom_line(  ) + 
  labs(title = "1cell-1")

this_plot


mouse_1cell_5_adjusted_start_site_coverage = shift_by_offsets(mouse_metagene_start, "1cell-3",  rep(15, 7 ) )

this_coverage = mouse_1cell_5_adjusted_start_site_coverage

this_df = data.frame(position = -50:50, count = this_coverage)

this_plot = ggplot(data = this_df, aes(x = position ,
                                       y = count )) + 
  geom_line(  ) + 
  labs(title = "1cell-1")

this_plot



mouse_1cell_5_adjusted_stop_site_coverage = shift_by_offsets(mouse_metagene_stop, "8cell-1", c(15, rep(15, 6 )) )

this_coverage = mouse_1cell_5_adjusted_stop_site_coverage

this_df = data.frame(position = -50:50, count = this_coverage)

this_plot = ggplot(data = this_df, aes(x = position ,
                                       y = count )) + 
  geom_line(  ) + 
  labs(title = "1cell-1")

this_plot


mouse_1cell_5_adjusted_start_site_coverage = shift_by_offsets(mouse_metagene_start, "8cell-1",  rep(15, 7 ) )

this_coverage = mouse_1cell_5_adjusted_start_site_coverage

this_df = data.frame(position = -50:50, count = this_coverage)

this_plot = ggplot(data = this_df, aes(x = position ,
                                       y = count )) + 
  geom_line(  ) + 
  labs(title = "1cell-1")

this_plot



human_100_1_adjusted_stop_site_coverage = shift_by_offsets(human_metagene_stop, "100-1", rep(15, 7 ) )

this_coverage = human_100_1_adjusted_stop_site_coverage

this_df = data.frame(position = -50:50, count = this_coverage)

this_plot = ggplot(data = this_df, aes(x = position ,
                                       y = count )) + 
  geom_line(  ) + 
  labs(title = "1cell-1")

this_plot


human_100_1_adjusted_start_site_coverage = shift_by_offsets(human_metagene_start, "100-1", rep(15, 7 ) )

this_coverage = human_100_1_adjusted_start_site_coverage

this_df = data.frame(position = -35:35, count = this_coverage[c( -1:-15, -86:-100  )  ])

this_title = "Stop Site"

this_plot = ggplot(data = this_df, aes(x = position ,
                                       y = count)) + 
  geom_line( color = BURNT_ORANGE , size = 1) + 
  labs(title = this_title, x = "", y = "frequency" ) + 
  theme(plot.title   =  element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid   = element_blank(),
        axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
        axis.line         = element_line(colour = "black", size = 0.35), ) +
  scale_x_continuous( name= "position", breaks =  c( -30, -15, 0, 15, 30)  )

this_plot

##################################################################################################
##################################################################################################


plot_metagene_unit = function(this_df, 
                              experiment, 
                              site, 
                              plot_title = "", 
                              offsets    = rep(15, 7 ),
                              y_upper    = 0 ){
  
  if(is.null(offsets) ){
    # This is for the external data
    this_coverage = get_experiment_coverage_at_given_len(
                                             df       = this_df, 
                                             this_exp = experiment, 
                                             len      = EXTERNAL_MIN_LENGTH)[c(-1, -2)] 
    this_coverage = as.integer( this_coverage )
    
    this_df = data.frame(position = -35:35, count = this_coverage)
  }
  else{
    this_coverage   = shift_by_offsets(this_df, experiment[1], offsets )

    if(length(experiment) >1){
      for(i in seq(2, length(experiment))){
        this_coverage = this_coverage + shift_by_offsets(this_df, experiment[i], offsets )
      }
    }

    
    this_df         = data.frame(position = -35:35, count = this_coverage[c( -1:-15, -86:-100  )  ])
  }
  
  
  if(site == "start"){
    this_color = START_SITE_COLOR
  }
  else{
    this_color = STOP_SITE_COLOR
  }
  
  if(y_upper == 0){
    y_upper = max(this_df$count)
  }
    
  if(y_upper > 75 & y_upper < 100){
    y_upper = 100
  }
  
  if(y_upper > 750 & y_upper < 1000){
    y_upper = 1000
  }
  
  expansion_factor = 0.02
  
  this_plot = ggplot(data = this_df, aes(x = position ,
                                         y = count)) + 
    geom_line( color = this_color , size = LINE_THICKNESS) + 
    labs(title = plot_title, x = "", y = "") +  
    theme(plot.title        =  element_text(hjust = 0.5, family = FIGURE_FONT, face = "plain", size = FONT_TITLE_SIZE, color = this_color),
          panel.border      = element_blank(),
          panel.background  = element_blank(),
          panel.grid        = element_blank(),
          axis.text.y       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.text.x       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.title.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
          axis.line         = element_line(colour = "black", size = 0.35), ) +
    scale_x_continuous( name= "", breaks =  c( -30, -15, 0, 15, 30)  )  +
    scale_y_continuous(limits = c(0, y_upper ), expand = expansion(mult = c(expansion_factor, 0)) )

  
  if(site == "stop" ){
    this_plot = this_plot + scale_y_continuous(position = "right", limits = c(0, y_upper ), expand = expansion(mult = c(expansion_factor, 0)))
  }
  
  return(this_plot)
  
}

##################################################################################################


combine_start_and_stop = function(start_plot, stop_plot, title){
  
  
  plot_title =  ggdraw() + 
    draw_label(
      title,
      fontface = 'plain',
      fontfamily = FIGURE_FONT,
      size = FONT_TITLE_SIZE, 
      angle = 90
    ) 
  
  this_plot = plot_grid( start_plot, stop_plot, plot_title, nrow = 1, rel_widths = c(1, 1, 0.02)  )
  
  #this_plot = plot_grid(plot_title, base_plot, nrow = 3, rel_widths = c(1, 1, 0.02))
  
  return(this_plot)
} 



combine_start_and_stop = function(start_plot, stop_plot, plot_title = ""){
  
  
  plot_title =  ggdraw() + 
    draw_label(
      plot_title,
      fontface = 'bold',
      fontfamily = FIGURE_FONT,
      size = FONT_TITLE_SIZE
    ) 
  
  base_plot = plot_grid( start_plot, stop_plot, nrow = 1   )
  
  this_plot = plot_grid(plot_title, base_plot, ncol = 1, rel_heights =  c(0.1, 1))
  
  return(this_plot)
} 



combine_plots_main = function(plot_upper, plot_lower){
  x_label =  ggdraw() + 
    draw_label(
      "position",
      fontface = 'plain',
      fontfamily = FIGURE_FONT,
      size = FONT_LABEL_SIZE
    ) 
  
  plot_title =  ggdraw() + 
    draw_label(
      "Metagene Coverage",
      fontface = 'bold',
      fontfamily = FIGURE_FONT,
      size = FONT_TITLE_SIZE
    ) 
  
  y_label =  ggdraw() + 
    draw_label(
      "frequency",
      fontface = 'plain',
      fontfamily = FIGURE_FONT,
      size = FONT_LABEL_SIZE,
      angle = 90
    ) 
  
  
  
  base_plot = plot_grid( plot_upper, plot_lower, ncol = 2)
  
  this_plot = plot_grid(y_label, base_plot, nrow = 1, rel_widths = c(0.05, 1)  )
  
  return(base_plot)
}


hundred_cell_start_plot = plot_metagene_unit(this_df    = human_metagene_start, 
                                             experiment = "100-1", 
                                             site       = "start", 
                                             plot_title = "Start Site",
                                             y_upper    = 1500 )

hundred_cell_stop_plot  = plot_metagene_unit(this_df    = human_metagene_stop, 
                                             experiment = "100-1", 
                                             site       = "stop", 
                                             plot_title = "Stop Site",
                                             y_upper = 4000)

hundred_cell_plot       = combine_start_and_stop(hundred_cell_start_plot, hundred_cell_stop_plot, "100 Cells")
hundred_cell_plot

million_cell_start_plot = plot_metagene_unit(this_df    = human_metagene_start, 
                                             experiment = "10M-3", 
                                             site       = "start" )

million_cell_stop_plot  = plot_metagene_unit(this_df    = human_metagene_stop, 
                                             experiment = "10M-3", 
                                             site       = "stop",
                                             y_upper    = 2500)

million_cell_plot       = combine_start_and_stop( million_cell_start_plot, 
                                                  million_cell_stop_plot, 
                                                  "10M Cells" ) 
million_cell_plot

human_metagene_plot = combine_plots_main(hundred_cell_plot, million_cell_plot)
human_metagene_plot

all_mouse_experiments       = unique(mouse_metagene_start$experiment)
one_cell_experiment_indices = grep("1cell", all_mouse_experiments )
one_cell_experiments        = all_mouse_experiments[one_cell_experiment_indices]

two_cell_experiment_indices = grep("2cell", all_mouse_experiments )
two_cell_experiments        = all_mouse_experiments[two_cell_experiment_indices]

four_cell_experiment_indices = grep("4cell", all_mouse_experiments )
four_cell_experiments        = all_mouse_experiments[four_cell_experiment_indices]

eight_cell_experiment_indices = grep("8cell", all_mouse_experiments )
eight_cell_experiments        = all_mouse_experiments[eight_cell_experiment_indices]


one_cell_combined_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                         experiment = one_cell_experiments, 
                                         site       = "start", 
                                         plot_title = "Start Site",
                                         y_upper    = 600)

one_cell_combined_stop_plot = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                                  experiment = one_cell_experiments, 
                                                  site       = "stop", 
                                                  plot_title = "Stop Site",
                                                  y_upper    = 4000 )


one_cell_combined_plot       = combine_start_and_stop( one_cell_combined_start_plot, 
                                                       one_cell_combined_stop_plot, 
                                              "1cell"  ) 

one_cell_combined_plot


two_cell_combined_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                                  experiment = two_cell_experiments, 
                                                  site       = "start", 
                                                  plot_title = "Start Site",
                                                  y_upper    = 500)



two_cell_combined_stop_plot = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                                 experiment = two_cell_experiments, 
                                                 site       = "stop", 
                                                 plot_title = "Stop Site",
                                                 y_upper    = 1500 )

two_cell_combined_plot       = combine_start_and_stop( two_cell_combined_start_plot, 
                                                       two_cell_combined_stop_plot, 
                                                       "2cell"  ) 

two_cell_combined_plot



four_cell_combined_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                                  experiment = four_cell_experiments, 
                                                  site       = "start", 
                                                  plot_title = "Start Site",
                                                  y_upper    = 400)



four_cell_combined_stop_plot = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                                 experiment = four_cell_experiments, 
                                                 site       = "stop", 
                                                 plot_title = "Stop Site",
                                                 y_upper    = 1200 )

four_cell_combined_plot       = combine_start_and_stop( four_cell_combined_start_plot, 
                                                       four_cell_combined_stop_plot, 
                                                       "4cell"  ) 

four_cell_combined_plot



eight_cell_combined_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                                   experiment = eight_cell_experiments, 
                                                   site       = "start", 
                                                   plot_title = "Start Site",
                                                   y_upper    = 1200)



eight_cell_combined_stop_plot = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                                  experiment = eight_cell_experiments, 
                                                  site       = "stop", 
                                                  plot_title = "Stop Site",
                                                  y_upper    = 4000 )

eight_cell_combined_plot       = combine_start_and_stop( eight_cell_combined_start_plot, 
                                                        eight_cell_combined_stop_plot, 
                                                        "8cell"  ) 

eight_cell_combined_plot


mouse_metagene_combined_plot = combine_plots_main(one_cell_combined_plot, eight_cell_combined_plot)
mouse_metagene_combined_plot


mouse_metagene_2_4_combined_plot = combine_plots_main(two_cell_combined_plot, four_cell_combined_plot)
mouse_metagene_2_4_combined_plot

################################################################################

one_cell_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                         experiment = "1cell-3", 
                                         site       = "start", 
                                         plot_title = "Start Site" )

one_cell_stop_plot  = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                         experiment = "1cell-3", 
                                         site       = "stop", 
                                         plot_title = "Stop Site",
                                         y_upper    = 600 )

one_cell_plot       = combine_start_and_stop( one_cell_start_plot, 
                                              one_cell_stop_plot, 
                                              "1cell"  ) 
one_cell_plot

eight_cell_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                           experiment = "8cell-1", 
                                           site       = "start",
                                           y_upper    = 200)

eight_cell_stop_plot  = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                           experiment = "8cell-1", 
                                           site       = "stop" )

eight_cell_plot       = combine_start_and_stop(eight_cell_start_plot, 
                                               eight_cell_stop_plot, "8cell")

mouse_metagene_plot = combine_plots_main(one_cell_plot, eight_cell_plot)

mouse_metagene_plot



two_cell_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                         experiment = "2cell-3", 
                                         site       = "start", 
                                         plot_title = "Start Site",
                                         y_upper    = 200 )

two_cell_stop_plot = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                         experiment = "2cell-3", 
                                         site       = "stop", 
                                         plot_title = "Stop Site" )

two_cell_plot = combine_start_and_stop(two_cell_start_plot,
                                       two_cell_stop_plot,
                                       "2cell")

four_cell_start_plot = plot_metagene_unit(this_df    = mouse_metagene_start, 
                                         experiment = "4cell-3", 
                                         site       = "start", 
                                         plot_title = "Start Site",
                                         y_upper    = 200 )

four_cell_stop_plot = plot_metagene_unit(this_df    = mouse_metagene_stop, 
                                        experiment = "4cell-3", 
                                        site       = "stop", 
                                        plot_title = "Stop Site",
                                        y_upper    = 500 )

four_cell_plot = combine_start_and_stop(four_cell_start_plot,
                                        four_cell_stop_plot,
                                        "4cell")

four_cell_plot

mouse_metagene_plot_2_4_cells = combine_plots_main(two_cell_plot, four_cell_plot)

mouse_metagene_plot_2_4_cells

################################################################################
######## E X T E R N A L    R N A - S E Q    D A T A         ###################

EXTERNAL_MIN_LENGTH = 35
EXTERNAL_MAX_LENGTH = 35

GSE101018_ribo_metagene_start = 
  get_metagene(GSE101018_ribo, 
               site        = "start",
               range.lower = EXTERNAL_MIN_LENGTH,
               range.upper = EXTERNAL_MAX_LENGTH,
               length      = FALSE,
               transcript  = TRUE,
               alias       = TRUE,
               compact     = FALSE,
               experiment  = GSE101018_ribo@experiments)

GSE101018_ribo_metagene_stop = 
  get_metagene(GSE101018_ribo, 
               site        = "stop",
               range.lower = EXTERNAL_MIN_LENGTH,
               range.upper = EXTERNAL_MAX_LENGTH,
               length      = FALSE,
               transcript  = TRUE,
               alias       = TRUE,
               compact     = FALSE,
               experiment  = GSE101018_ribo@experiments)


GSE78634_ribo_metagene_start = 
  get_metagene(GSE78634_ribo, 
               site        = "start",
               range.lower = EXTERNAL_MIN_LENGTH,
               range.upper = EXTERNAL_MAX_LENGTH,
               length      = FALSE,
               transcript  = TRUE,
               alias       = TRUE,
               compact     = FALSE,
               experiment  = GSE78634_ribo@experiments)

GSE78634_ribo_metagene_stop = 
  get_metagene(GSE78634_ribo, 
               site        = "stop",
               range.lower = EXTERNAL_MIN_LENGTH,
               range.upper = EXTERNAL_MAX_LENGTH,
               length      = FALSE,
               transcript  = TRUE,
               alias       = TRUE,
               compact     = FALSE,
               experiment  = GSE78634_ribo@experiments)




GSE101018_1_start_plot = plot_metagene_unit(this_df    = GSE101018_ribo_metagene_start, 
                                            experiment = GSE101018_ribo@experiments[1], 
                                            site       = "start", 
                                            plot_title = "Start Site",
                                            offsets    = c(),
                                            y_upper    = 60000)



GSE101018_1_stop_plot = plot_metagene_unit(this_df    = GSE101018_ribo_metagene_stop, 
                                            experiment = GSE101018_ribo@experiments[1], 
                                            site       = "stop", 
                                            plot_title = "Stop Site",
                                            offsets    = c(),
                                            y_upper    = 25000)




GSE101018_1_metagene_plot = combine_start_and_stop(GSE101018_1_start_plot, 
                                                   GSE101018_1_stop_plot, 
                                                   GSE101018_ribo@experiments[1])


#GSE101018_1_metagene_plot

GSE101018_2_start_plot = plot_metagene_unit(this_df    = GSE101018_ribo_metagene_start, 
                                            experiment = GSE101018_ribo@experiments[2], 
                                            site       = "start", 
                                            plot_title = "",
                                            offsets    = c(),
                                            y_upper    = 30000)



GSE101018_2_stop_plot = plot_metagene_unit(this_df    = GSE101018_ribo_metagene_stop, 
                                           experiment = GSE101018_ribo@experiments[2], 
                                           site       = "stop", 
                                           plot_title = "",
                                           offsets    = c(),
                                           y_upper    = 10000)


GSE101018_2_metagene_plot = combine_start_and_stop(GSE101018_2_start_plot, 
                                                   GSE101018_2_stop_plot, 
                                                   GSE101018_ribo@experiments[2])

#GSE101018_2_metagene_plot


GSE101018_metagene_combined_plot = 
  combine_plots_main(GSE101018_1_metagene_plot, GSE101018_2_metagene_plot)

GSE101018_metagene_combined_plot


GSE78634_1_start_plot = plot_metagene_unit(this_df    = GSE78634_ribo_metagene_start, 
                                            experiment = GSE78634_ribo@experiments[1], 
                                            site       = "start", 
                                            plot_title = "Start Site",
                                            offsets    = c(),
                                            y_upper    = 30000)

GSE78634_1_stop_plot = plot_metagene_unit(this_df    = GSE78634_ribo_metagene_stop, 
                                           experiment = GSE78634_ribo@experiments[1], 
                                           site       = "stop", 
                                           plot_title = "Stop Site",
                                           offsets    = c(),
                                           y_upper    = 10000)

GSE78634_1_metagene_plot = combine_start_and_stop(GSE78634_1_start_plot,
                                                  GSE78634_1_stop_plot,
                                                  GSE78634_ribo@experiments[1])

#GSE78634_1_metagene_plot


GSE78634_2_start_plot = plot_metagene_unit(this_df    = GSE78634_ribo_metagene_start, 
                                           experiment = GSE78634_ribo@experiments[2], 
                                           site       = "start", 
                                           plot_title = "",
                                           offsets    = c(),
                                           y_upper    = 0)

GSE78634_2_stop_plot = plot_metagene_unit(this_df    = GSE78634_ribo_metagene_stop, 
                                          experiment = GSE78634_ribo@experiments[2], 
                                          site       = "stop", 
                                          plot_title = "",
                                          offsets    = c(),
                                          y_upper    = 4000)

GSE78634_2_metagene_plot = combine_start_and_stop(GSE78634_2_start_plot,
                                                  GSE78634_2_stop_plot,
                                                  GSE78634_ribo@experiments[2])

#GSE78634_2_metagene_plot

GSE78634_combined_metagene_plot = combine_plots_main(GSE78634_1_metagene_plot,
                                                     GSE78634_2_metagene_plot)

GSE78634_combined_metagene_plot

################################################################################


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

## Also see https://www.sciencemag.org/sites/default/files/Figure_prep_guide.pdf
##
## Science articles are 3 columns
## One column figure width is 2.25.  So 2 inches would be safe
## Two column is 4.75 so 4.5 inches would be safe

# Figure and legend should fit in  frame of 7.25 x 8.9 inches
# If we consider the page divided into two columns, then, half of it is 3.54 inches
# Two columns together is 7.25 inches (there is some space between the columns)

PDF_WIDTH  = 3.54
PDF_HEIGHT = PDF_WIDTH



save_plot_pdf("mouse_combined_metagene.pdf", mouse_metagene_combined_plot , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)

save_plot_pdf("mouse_metagene_2_4_combined_plot.pdf", mouse_metagene_2_4_combined_plot , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)

save_plot_pdf("mouse_metagene.pdf", mouse_metagene_plot , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)

save_plot_pdf("mouse_metagene_2_4_cells.pdf", mouse_metagene_plot_2_4_cells , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)

save_plot_pdf("GSE78634_metagene.pdf", GSE78634_combined_metagene_plot , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)

save_plot_pdf("GSE101018_metagene.pdf", GSE101018_metagene_combined_plot , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)


save_plot_pdf("human_metagene.pdf", human_metagene_plot , 
              width  = PDF_WIDTH, 
              height = PDF_HEIGHT)
