
#### QC FIGURE LAYOUT

source("./correlation.R")

source('./metagene.R')

source('./region_counts.R')


#####################################
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


## Scatter plots from human data
squish_factor_1 = 0.58
page_top      = plot_grid(sp_grid, separator, rel_widths = c(1,squish_factor_1))
mouse_sp_grid_modified = plot_grid(mouse_sp_grid, separator, rel_widths = c(1,squish_factor_1))


mouse_metagene_plot_modified = plot_grid(mouse_metagene_plot, separator, rel_widths = c(1,0.37))

mouse_region_counts_comperative_plot_modified = plot_grid( mouse_region_counts_comperative_plot, separator, rel_widths = c(1,0.34)  )

#bottom_legend = get_legend(ppig_stages_plot)

page_bottom = plot_grid(genes_accross_stages_main_figure_bar, 
                        separator,
                        rel_widths = c(1,0.25))

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


#metagene_and_region_counts = plot_grid(mouse_metagene_plot ,  mouse_region_counts_comperative_plot, nrow = 1 )

figure_layout = plot_grid(page_top,
                          schematic_label,
                          mouse_metagene_plot_modified,
                          mouse_region_counts_comperative_plot_modified,
                          mouse_sp_grid_modified, 
                          page_bottom,
                          ncol = 1,
                          rel_heights = c(0.9, 1, 1, 0.9, 1),
                          #labels = c("A", "B", "C", "D", "E"),
                          labels = "AUTO",
                          label_size = 10,
                          label_fontfamily = "helvetica",
                          align = "hv")

#figure_layout



save_plot_pdf("qc_figure_layout.pdf", figure_layout, width = unit(7.20, "in"), height = unit(9, "in"))


ggsave("pdf/figure_layout.png", 
       plot   = figure_layout, 
       device = png(), 
       width  = unit(7.20, "in"),
       height = unit(9, "in"),
       dpi    = 600 )


figure_layout_bottom_only = 
  plot_grid(separator,
            schematic_label,
            separator,
            separator,
            separator, 
            page_bottom,
            ncol = 1,
            rel_heights = c(0.9, 1, 1, 0.9, 1),
            #labels = c("A", "B", "C", "D", "E"),
            labels = "AUTO",
            label_size = 10,
            label_fontfamily = "helvetica",
            align = "hv")

save_plot_pdf("qc_figure_layout_bottom.pdf", 
              figure_layout_bottom_only, 
              width = unit(7.20, "in"), 
              height = unit(9, "in"))