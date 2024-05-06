###########################
# libraries required here #
###########################
library( ggrepel )
library( tidyverse )
library( RColorBrewer )
library( pheatmap )

###############################
# color blind friendly colors #
###############################
colorBlindPalette <- c( "#999999",
                        "#E69F00",
                        "#56B4E9",
                        "#009E73",
                        "#F0E442",
                        "#0072B2",
                        "#D55E00",
                        "#CC79A7" )

colorBlindPalette2 <- c( "#000000",
                         "#999999",
                         "#E69F00",
                         "#56B4E9",
                         "#009E73",
                         "#F0E442",
                         "#0072B2",
                         "#D55E00",
                         "#CC79A7" )

blues5Palette <- c( '#ffffcc',
                    '#a1dab4',
                    '#41b6c4',
                    '#2c7fb8',
                    '#253494' )

greens5Palette <- c( '#ffffcc',
                     '#c2e699',
                     '#78c679',
                     '#31a354',
                     '#006837' )

purples5Palette <- c( '#feebe2',
                      '#fbb4b9',
                      '#f768a1',
                      '#c51b8a',
                      '#7a0177' )

reds5Palette <- c( '#ffffb2',
                   '#f3cc5c',
                   '#fd8d3c',
                   '#f03b20',
                   '#bd0026' )

reds4Palette <- c( '#fef0d9',
                   '#fdcc8a',
                   '#fc8d59',
                   '#d7301f' )

reds3Palette <- c( '#fee0d2',
                   '#fc9272',
                   '#de2d26' )

####################
# ggplot modifiers #
####################
raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text = element_text(size = 16),
  legend.position = "none",
  plot.title = element_text(size = 20, hjust = 0.5),
  #panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(
    colour = 'black',
    #size=0.5,
    linewidth = 0.5,
    linetype = 'solid'
  ),
  axis.line.y = element_line(
    colour = 'black',
    #size=0.5,
    linewidth = 0.5,
    linetype = 'solid'
  )
)

gg_bigger_texts = theme(
  axis.title = element_text(size = 22),
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 15),
  plot.title = element_text(size = 22),
  strip.text.x = element_text(size = 17,
                              margin = margin(b = 5, t = 5)),
  strip.text.y = element_text(size = 15)
)

gg_multiplot_texts = theme(
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 18),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 13),
  plot.title = element_text(size = 20),
  strip.text.x = element_text(size = 16,
                              margin = margin(b = 5, t = 5)),
  strip.text.y = element_text(size = 15)
)

gg_quadplot_smaller_text = theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 9),
  plot.title = element_text(size = 15)
)

gg_reduce_pathway_text = theme(
  axis.title = element_text(size = 14),
  axis.text.y = element_text(size = 8),
  axis.text.x = element_text(size = 10),
  plot.title = element_text(size = 15)
)

gg_no_legend = theme(legend.position = 'none')

gg_no_grid = theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())

gg_no_x_grid = theme(panel.grid.major.x = element_blank())

gg_no_y_grid = theme(panel.grid.major.y = element_blank())

gg_center_title = theme(plot.title = element_text(hjust = 0.5))

gg_no_x_label = theme(axis.title.x = element_blank())

gg_no_y_label = theme(axis.title.y = element_blank())

gg_angled_x_text = theme (axis.text.x = element_text(
  angle = 45,
  vjust = 1,
  hjust = 1,
  color = 'black'
))

#################
# volcano plots #
#################
make_ggplot_volcano <- function( deg_dataframe,
                                 case_name,
                                 control_name,
                                 axis_steps = 2,
                                 fold_change_cutoff = 1.5,
                                 qvalue_cutoff = 0.05,
                                 max_label = 30,
                                 truncate_fc = Inf )
{
  ##############################
  # set significance threshold #
  ##############################
  deg_dataframe <- deg_dataframe %>%
    mutate( Significant = case_when(
      qval < qvalue_cutoff & abs( FoldChange ) >= fold_change_cutoff ~ "Large",
      qval < qvalue_cutoff ~ "Modest",
      TRUE ~ "Not" ) ) %>%
    mutate( Significant = factor( Significant,
                                  levels=c( "Not", "Modest", "Large" ) ) ) %>%
    mutate( log2FoldChange = case_when(
      log2FoldChange > truncate_fc ~ Inf,
      log2FoldChange < (-1 * truncate_fc) ~ -Inf,
      TRUE ~ log2FoldChange ) )

  ################################
  # set values for square x axis #
  ################################
  x_volcano_value <- ( abs( deg_dataframe$log2FoldChange[ is.finite( deg_dataframe$log2FoldChange ) ] ) + 0.051 ) %>%
    max( . ) %>%
    round( ., 1 )

  if ( x_volcano_value < 1.0 ) {
    x_volcano_value = 1.0
  }

  x_num_for_limits <- round( x_volcano_value, 0 )

  x_volcano_low <- x_volcano_value * -1
  x_volcano_high <- x_volcano_value

  x_break_list <- seq( -1 * x_num_for_limits, x_num_for_limits, by = axis_steps )

  ##############
  # plot lines #
  ##############
  horizontal_line <- log10( qvalue_cutoff ) * -1
  vertical_line_1 <- log2( fold_change_cutoff )
  vertical_line_2 <- vertical_line_1 * -1

  ###################################
  # actually make the volcano plots #
  ###################################
  plot_volcano <- ggplot( deg_dataframe,
                          aes( x=log2FoldChange,
                               y=-log10( qval ),
                               colour=Significant ) ) +
    scale_colour_manual( values = c( "darkgray",
                                     "blue",
                                     "red" ) ) +
    scale_x_continuous( limits = c( x_volcano_low, x_volcano_high ),
                        breaks = x_break_list ) +
    theme_bw() +
    gg_bigger_texts +
    gg_no_legend +
    gg_no_grid +
    gg_center_title +
    geom_point( size=1.2 ) +
    geom_hline( yintercept = horizontal_line,
                linetype=2 ) +
    geom_vline( xintercept = c( vertical_line_1, vertical_line_2 ),
                linetype = 2 ) +
    geom_text_repel( data = subset( deg_dataframe, Significant == "Large" )[c(1:max_label),],
                     colour = "black",
                     aes( label = symbol ),
                     size = 3 ) +
    xlab( parse( text = paste0( "log[2]~(", case_name, "/", control_name, ")" ) ) ) +
    ylab( parse( text = paste0( "-log[10]~(Adj.~p-value)" ) ) )

  return( plot_volcano )
}

############
# make PCA #
############
make_ggplot_pca <- function( input_data,
                             sample_info,
                             join_name = 'sample_id',
                             first_plot_pc = 1,
                             second_plot_pc = 2 ) {
  pca_out <- prcomp( x = t( input_data ),
                     scale. = FALSE )

  percent_variance <- round( pca_out$sdev^2 / sum( pca_out$sdev^2 ) * 100.0, 2 )

  pca_coords <- pca_out$x %>%
    as.data.frame( . ) %>%
    rownames_to_column( join_name ) %>%
    full_join( y = sample_info,
               by = join_name )

  # Plot dimension names
  x_dim_name <- paste0( "PC",
                        first_plot_pc,
                        " (",
                        percent_variance[first_plot_pc],
                        "%)" )

  y_dim_name <- paste0( "PC",
                        second_plot_pc,
                        " (",
                        percent_variance[second_plot_pc],
                        "%)" )

  ggplot( data = pca_coords,
          mapping = aes( x = PC1,
                         y = PC2,
                         shape = status,
                         colour = status ) ) +
    theme_bw() +
    geom_point( size = 5 ) +
    geom_text_repel( colour = "black",
                     aes( label = sample_id ),
                     size = 3 ) +
    gg_bigger_texts +
    gg_center_title +
    gg_no_grid +
    xlab( x_dim_name ) +
    ylab( y_dim_name ) +
    theme( legend.position = 'top' ) %>%
    return( . )
}

################
# make heatmap #
################
make_heatmap <- function( input_data,
                          sample_info,
                          showrownames = TRUE,
                          showcolnames = TRUE,
                          cutreerow = 2,
                          cutreecol = 2,
                          color_annotation = FALSE,
                          custom_palette = FALSE,
                          color_ramp_number = 20,
                          cluster_columns = TRUE,
                          column_gaps = NULL ) {

  heatmap_info <- sample_info %>%
    select( sample_id, status ) %>%
    distinct( . ) %>%
    column_to_rownames( 'sample_id' ) %>%
    mutate( status = factor( status, levels = c( 'Untreated', 'Treated' ) ) )

  heatmap_info = heatmap_info[ colnames( input_data ), , drop = FALSE ]

  if ( color_annotation == FALSE ) {
    col_ann_colors <- list(
      status = c( Untreated = "#f0f0f0", Treated = "#bdbdbd" )
    )
  } else {
    col_ann_colors <- list(
      status = c( Untreated = "#00BFC4", Treated = "#F8766D" )
    )
  }

  if ( custom_palette == TRUE ){
    this_palette <- colorRampPalette( c( "#4373b4", "#fdf8b5", "#d52a22" ) )(color_ramp_number)
  } else {
    this_palette <- colorRampPalette( rev( brewer.pal( n=10, "RdBu" ) ) )(color_ramp_number)
  }

  pheatmap( mat = input_data,
            color = this_palette,
            cluster_cols = cluster_columns,
            gaps_col = column_gaps,
            show_rownames = showrownames,
            show_colnames = showcolnames,
            annotation_col = heatmap_info,
            annotation_colors = col_ann_colors,
            scale = 'row',
            cutree_rows = cutreerow,
            cutree_cols = cutreecol,
            fontsize_row = 6 )
}
