
source("./volcano_plots.R")
#source("./differential_analysis.R")
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


### To be FIXED!!!!

point_plot_legend = get_legend(Ppig_point_plot)

point_plot_grid = plot_grid(
  Ppig_point_plot + theme(legend.position = "None"),
  Ppig_point_plot + theme(legend.position = "None"),
  Ppig_point_plot + theme(legend.position = "None"),
  point_plot_legend,
  nrow = 1,
  rel_widths = c(1,1,1,0.8)
)

point_plot_grid

################################################################################


separator = generate_blank_plot()

figure_left = plot_grid(distinct_genes_heatmap_from_normalization_averaged [[4]] , 
                        separator, 
                        rel_heights = c(1, 0.5),
                        nrow = 2  )

figure_right = plot_grid(volcano_main_plot, 
                         point_plot_grid,
                         separator,
                         ncol = 1,
                         rel_heights = c(1, 0.6, 2))


figure_layout = plot_grid(figure_left,
                          figure_right,
                          ncol = 2,
                          rel_widths = c(1, 2.8))

save_plot_pdf("diff_figure_layout.pdf", figure_layout, width = unit(7.20, "in"), height = unit(9, "in"))


figure_right = plot_grid(volcano_main_plot, 
                         separator,
                         separator,
                         ncol = 1,
                         rel_heights = c(1, 0.6, 2))

figure_layout_volcano = plot_grid(separator,
                          figure_right,
                          ncol = 2,
                          rel_widths = c(1, 2.8))

save_plot_pdf("diff_figure_layout_volcano.pdf", figure_layout_volcano, width = unit(7.20, "in"), height = unit(9, "in"))
