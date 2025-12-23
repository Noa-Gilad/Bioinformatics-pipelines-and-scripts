# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Function to clean pathway names
clean_pathway_name <- function(term, category) {
  if (grepl("KEGG", category)) {
    # Remove KEGG pathway ID (e.g., "hsa00030:")
    cleaned <- str_replace(term, "^hsa\\d+:", "")
  } else if (grepl("GOTERM", category)) {
    # Remove GO ID (e.g., "GO:0014829~")
    cleaned <- str_replace(term, "^GO:\\d+~", "")
  } else {
    cleaned <- term
  }
  
  # Clean up and capitalize first letter
  cleaned <- str_trim(cleaned)
  if (nchar(cleaned) > 0) {
    cleaned <- paste0(toupper(substr(cleaned, 1, 1)), substr(cleaned, 2, nchar(cleaned)))
  }
  
  return(cleaned)
}

# Function to create horizontal bar plot matching the exact style
create_horizontal_bar_plot <- function(df, title = "Pathway Enrichment", max_terms = 20, padj = NULL) {
  
  legend_name <- "p.adjust"  # Default legend title
  
  # Filter by padj if specified
  if (!is.null(padj)) {
    # First try Benjamini
    keep_adj <- !is.na(df$Benjamini) & df$Benjamini <= padj
    
    if (any(keep_adj)) {
      df <- df[keep_adj, , drop = FALSE]
      legend_name <- "p.adjust"
    } else {
      keep_raw <- !is.na(df$PValue) & df$PValue <= padj
      df <- df[keep_raw, , drop = FALSE]
      legend_name <- "PValue"
      cat("Note: No adjusted p-values <= cutoff; filtering by raw p-values instead\n")
    }
  }
  
  # Sort by p-value and take top terms
  df_sorted <- df %>%
    arrange(PValue) %>%
    head(max_terms)
  
  # Clean pathway names
  df_sorted$CleanTerm <- mapply(clean_pathway_name, df_sorted$Term, df_sorted$Category)
  
  # Assign p.adjust or PValue for coloring
  if (legend_name == "PValue") {
    df_sorted$plot_value <- df_sorted$PValue
    cat("Note: Using raw p-values for coloring\n")
  } else {
    df_sorted$plot_value <- df_sorted$Benjamini
    # If all Benjamini values are 1.0, switch to PValue
    if (length(unique(df_sorted$plot_value)) == 1 && unique(df_sorted$plot_value) == 1.0) {
      df_sorted$plot_value <- df_sorted$PValue
      legend_name <- "PValue"
      cat("Note: Using raw p-values as all Benjamini values are 1.0\n")
    }
  }
  
  # Reorder terms for plotting
  df_sorted$CleanTerm <- factor(df_sorted$CleanTerm, 
                                levels = rev(df_sorted$CleanTerm[order(df_sorted$PValue)]))
  
  # Plot
  p <- ggplot(df_sorted, aes(x = Count, y = CleanTerm, fill = plot_value)) +
    geom_col(width = 0.7, color = "white", size = 0.1) +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      name = legend_name,
      labels = function(x) {
        ifelse(x >= 0.01, sprintf("%.2f", x), 
               ifelse(x >= 0.001, sprintf("%.3f", x), 
                      formatC(x, format = "e", digits = 1)))
      },
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 0.8,
        barheight = 6,
        frame.colour = "black",
        frame.linewidth = 0.5,
        ticks.colour = "black",
        ticks.linewidth = 0.5
      )
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, max(df_sorted$Count) * 1.05),
      breaks = pretty(c(0, max(df_sorted$Count)), n = 6)
    ) +
    labs(
      title = title,
      x = "Count",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black", hjust = 1),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.margin = margin(l = 15),
      plot.margin = margin(15, 15, 15, 15),
      panel.border = element_blank()
    )
  
  return(p)
}


# Function to load and plot DAVID data
load_and_plot_david_data <- function(filepath, title_suffix = "", max_terms = 20, padj = NULL) {
  
  # Load the data
  df <- read.delim(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Determine if it's KEGG or GO data for title
  if (grepl("KEGG", df$Category[1])) {
    plot_title <- paste("KEGG Pathway Enrichment", title_suffix)
  } else if (grepl("GO", df$Category[1])) {
    plot_title <- paste("GO Term Enrichment", title_suffix)
  } else {
    plot_title <- paste("Pathway Enrichment", title_suffix)
  }
  
  # Create the plot
  p <- create_horizontal_bar_plot(df, title = plot_title, max_terms = max_terms, padj = padj)
  
  return(list(plot = p, data = df))
}

