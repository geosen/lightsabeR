#'Plot Signature strength over replicate correlation based on a Cmap2 results table
#'
#'Dependencies are:
#' library(cmapR)
#' library(tidyverse)
#' library(ggrepel)
#' 
#' 
#'date 22/04/2024
#'@param cmap_melted A melted data.frame containing Connectivity Map 2 results
#'@param n_top_highlight Number of top results to label
#'@param cell_types Cell types to keep
#'@param ids_to_highlight Particular perturbation ids to label. Overrides n_top_highlight if used.
#'@param max.overlaps geom_label_repel max overlaps parameter
#'@param min.segment.length geom_label_repel minimum segment length parameter
#'@param label_alpha Transparency parameter for labels
#'@param label_size Size parameter for labels
#'@param tas_color Color of the tas lines
#'@param tas_linetype Linetype of the tas lines
#'
#'
#'
#'@return A ggplot
#'@export

cmap2_cc_ss_plot <- function(cmap_melted, n_top_highlight = NULL, cell_types = NULL,
                             ids_to_highlight = NULL, max.overlaps = 40,
                             min.segment.length = unit(0.001,'mm'),
                             label_alpha = 0.8,
                             label_size = 3,
                             tas_color = 'blue',
                             tas_linetype = 'dashed') {
  
  if(!requireNamespace("cmapR")) {
    stop('Please load cmapR to continue')
  }
  if(!requireNamespace("tidyverse")) {
    stop('Please load tidyverse to continue')
  }
  if(!requireNamespace("purrr")) {
    stop('Please load purrr to continue')
  }
#tas_function <- function(x,y) {sqrt(x*y/978)}

#annotating signatures
df1 <- df1 %>% 
  mutate(Signature = case_when(cc_q75 >= 0.2 & ss_ngene <200 ~ 'Subtle & reproducible',
                                      cc_q75 >= 0.2 & ss_ngene >=200 ~ 'Strong & reproducible',
                                      cc_q75 < 0.2 & ss_ngene <200 ~ 'Inert',
                                      cc_q75 < 0.2 & ss_ngene >=200 ~ 'Noisy')) 

if(!is.null(cell_types)) {
df2 <- df1 %>%
  filter(cell_iname %in% cell_types)
} else {
  df2 <- df
}

if(!is.null(n_top_highlight) & is.null(ids_to_highlight)) {
ids_to_highlight <- df2 %>%
  slice_max(tas,n = n_top_highlight) %>%
  pull(id)
} 

geom_tas <- map(
  seq(0.2, 0.8, 0.2), 
  ~ geom_function(fun = function(x) (978 * .x^2) / x, 
                  color = tas_color, 
                  linetype = tas_linetype
  )
)





g <- ggplot(df2, aes(x = cc_q75,y = ss_ngene, color = Signature)) + 
  geom_hline(yintercept = 200) + 
  geom_vline(xintercept = 0.2)+
  geom_rect(xmin = 0, xmax = 0.2, ymin = 0, ymax = 200,fill = 'lavenderblush', 
            color = 'grey80', alpha = 0.4) + 
  geom_rect(xmin = 0.2, xmax = 1, ymin = 0, ymax = 200,fill = 'darkseagreen1', 
            color = 'grey80', alpha = 0.4) + 
  geom_rect(xmin = 0, xmax = 0.2, ymin = 200, ymax = 1000,fill = 'aliceblue', 
            color = 'grey80', alpha = 0.4) + 
  geom_rect(xmin = 0.2, xmax = 1, ymin = 200, ymax = 1000,fill = 'lemonchiffon', 
            color = 'grey80', alpha = 0.4) + 
  geom_point() + 
  geom_tas + 
  labs(y = 'Signature strength',
       x = 'Replicate correlation') + 
  scale_color_manual(values = c('thistle3','cadetblue3','goldenrod3','darkolivegreen3')) + 
  scale_y_continuous(
    limits = c(0, 1000), 
    sec.axis = sec_axis(~ sqrt(. / 978), "TAS", 
                        breaks = seq(0.2, 0.8, 0.2))) + 
  theme_minimal()



if(!is.null(ids_to_highlight)) {
  if(!requireNamespace("ggrepel")) {
    stop('Please load ggrepel to continue')
  }
g <- g + 
    geom_label_repel(data = df2 %>%
                       filter(id %in% ids_to_highlight),
                     aes(label = pert_iname),
                     max.overlaps = max.overlaps,
                     min.segment.length = min.segment.length, 
                     show.legend = F, 
                     alpha = label_alpha, 
                     size = label_size)
}
g

}
