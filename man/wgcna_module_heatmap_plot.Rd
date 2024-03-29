% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wgcna_module_heatmap_plot.R
\name{wgcna_module_heatmap_plot}
\alias{wgcna_module_heatmap_plot}
\title{Plot Module - Trait correlation using the output of WGCNA}
\usage{
wgcna_module_heatmap_plot(
  wgcna_trait_res,
  row_labels = NULL,
  col_labels = NULL,
  row_label_size = 14,
  col_label_size = 7,
  colf = NULL,
  row_angle = 0,
  column_km = 1,
  row_km = 1,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  title = "WGCNA Module - Trait Correlation",
  legend_title = "Correlation",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  border_linetype = 1,
  border_linewidth = 1,
  cell_color = "white",
  cell_linetype = 1,
  cell_linewidth = 1
)
}
\arguments{
\item{wgcna_trait_res}{Input wgcna results data.frame}

\item{row_labels}{A vector to name the row labels}

\item{col_labels}{A vector to name the column labels}

\item{row_label_size}{The size of the row labels}

\item{col_label_size}{The size of the column labels}

\item{colf}{A function for the heatmap colors}

\item{row_angle}{The angle of the x axis labels}

\item{column_km}{How many clusters should the columns be divided into}

\item{row_km}{How many clusters should the rows be divided into}

\item{show_row_dend}{Whether to show the row dendrogram}

\item{show_column_dend}{Whether to show the column dendrogram}

\item{title}{The title of the plot}

\item{legend_title}{The title of the legend}

\item{cluster_rows}{Whether to cluster rows or not}

\item{cluster_cols}{Whether to cluster columns or not}

\item{border_color}{Outer rectangle line color}

\item{border_linetype}{Outer rectangle linetype}

\item{border_linewidth}{Outer rectangle linewidth}

\item{cell_color}{Inner rectangle line color}

\item{cell_linetype}{Inner rectangle linetype}

\item{cell_linewidth}{Inner rectangle linewidth}
}
\value{
A ComplexHeatmap
}
\description{
Dependencies are:
library(ComplexHeatmap)
library(circlize)
library(stringr)
}
\details{
date 11/07/2023
}
