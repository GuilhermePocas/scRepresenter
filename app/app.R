
library(shiny)
library(plotly)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(networkD3)
library(anndataR)
library(RColorBrewer)
library(htmlwidgets)

options(shiny.maxRequestSize = 10 * 1024^3)  # 10 GB

# =========================
# Helper functions
# =========================

safe_char <- function(x) {
  if (is.factor(x)) return(as.character(x))
  if (is.list(x)) return(vapply(x, function(z) paste(z, collapse = ";"), character(1)))
  as.character(x)
}

get_meta_df <- function(sce) {
  meta <- as.data.frame(colData(sce))
  meta$cell_id <- colnames(sce)
  meta
}

get_valid_grouping_cols <- function(meta_df) {
  cols <- colnames(meta_df)
  valid <- cols[sapply(cols, function(cl) {
    x <- meta_df[[cl]]
    if (cl == "cell_id") return(FALSE)
    if (is.list(x)) return(FALSE)
    ux <- unique(safe_char(x))
    ux <- ux[!is.na(ux) & ux != ""]
    length(ux) >= 2
  })]
  valid
}

detect_prediction_cols <- function(meta_df) {
  cols <- colnames(meta_df)
  preds <- cols[grepl("prediction$", cols, ignore.case = TRUE)]
  preds
}

make_discrete_palette <- function(values) {
  vals <- sort(unique(values[!is.na(values) & values != ""]))
  n <- length(vals)
  if (n == 0) return(NULL)
  
  base_cols <- c(
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Set3"),
    brewer.pal(9, "Set1")
  )
  
  if (n > length(base_cols)) {
    cols <- grDevices::colorRampPalette(base_cols)(n)
  } else {
    cols <- base_cols[seq_len(n)]
  }
  
  names(cols) <- vals
  cols
}

downsample_cells <- function(df, mode = "All", n = 3000, stratify_col = NULL) {
  if (mode == "All" || nrow(df) <= n) return(df)
  
  set.seed(123)
  
  if (mode == "Random") {
    keep <- sample(seq_len(nrow(df)), n)
    return(df[keep, , drop = FALSE])
  }
  
  if (mode == "Stratified" && !is.null(stratify_col) && stratify_col %in% colnames(df)) {
    df2 <- df %>%
      mutate(.grp = safe_char(.data[[stratify_col]])) %>%
      filter(!is.na(.grp), .grp != "")
    
    if (nrow(df2) == 0) return(df)
    
    props <- table(df2$.grp) / nrow(df2)
    target_n <- pmax(1, round(props * n))
    
    sampled <- lapply(names(target_n), function(g) {
      idx <- which(df2$.grp == g)
      take <- min(length(idx), target_n[g])
      df2[sample(idx, take), , drop = FALSE]
    })
    
    out <- bind_rows(sampled)
    out$.grp <- NULL
    return(out)
  }
  
  df
}

make_embedding_df <- function(sce, reduction_name, color_by = NULL) {
  rd <- reducedDim(sce, reduction_name)
  
  if (is.null(rd)) stop(paste("Reduction not found:", reduction_name))
  if (ncol(rd) < 2) stop(paste("Reduction", reduction_name, "has fewer than 2 dimensions."))
  
  rd <- as.data.frame(rd)
  colnames(rd)[1:2] <- c("Dim1", "Dim2")
  rd$cell_id <- colnames(sce)
  
  meta <- get_meta_df(sce)
  df <- left_join(rd, meta, by = "cell_id")
  
  if (!is.null(color_by) && color_by %in% colnames(df)) {
    df$color_value <- df[[color_by]]
  } else {
    df$color_value <- "Cells"
  }
  
  df
}

make_embedding_plot <- function(df, title_text, color_by, pt_size = 0.8, pt_alpha = 0.85) {
  hover_text <- paste0("Cell: ", df$cell_id)
  if (!is.null(color_by) && color_by %in% colnames(df)) {
    hover_text <- paste0(
      hover_text,
      "<br>", color_by, ": ", safe_char(df[[color_by]])
    )
  }
  
  p <- ggplot(df, aes(x = Dim1, y = Dim2, text = hover_text))
  
  if (is.numeric(df$color_value)) {
    p <- p +
      geom_point(aes(color = color_value), size = pt_size, alpha = pt_alpha) +
      scale_color_viridis_c()
  } else {
    df$color_value <- safe_char(df$color_value)
    pal <- make_discrete_palette(df$color_value)
    p <- p +
      geom_point(aes(color = color_value), size = pt_size, alpha = pt_alpha) +
      scale_color_manual(values = pal, drop = FALSE)
  }
  
  p +
    labs(
      title = title_text,
      x = "Dim 1",
      y = "Dim 2",
      color = color_by
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.title = element_text(face = "bold")
    )
}

make_sankey_debug_df <- function(meta_df, labelA, labelB,
                                 filter_col = NULL, filter_vals = NULL,
                                 min_freq = 1, top_n = 100) {
  req(labelA %in% colnames(meta_df), labelB %in% colnames(meta_df))
  
  df <- meta_df
  
  if (!is.null(filter_col) &&
      filter_col != "None" &&
      filter_col %in% colnames(df) &&
      !is.null(filter_vals) &&
      length(filter_vals) > 0) {
    df <- df[df[[filter_col]] %in% filter_vals, , drop = FALSE]
  }
  
  A <- unname(safe_char(df[[labelA]]))
  B <- unname(safe_char(df[[labelB]]))
  
  keep <- !is.na(A) & !is.na(B) & A != "" & B != ""
  A <- A[keep]
  B <- B[keep]
  
  if (length(A) == 0 || length(B) == 0) {
    return(data.frame(A = character(0), B = character(0), Freq = numeric(0)))
  }
  
  tab <- as.data.frame(table(A, B), stringsAsFactors = FALSE)
  colnames(tab) <- c("A", "B", "Freq")
  tab <- tab[tab$Freq >= min_freq, , drop = FALSE]
  tab <- tab[order(-tab$Freq), , drop = FALSE]
  if (nrow(tab) > top_n) tab <- tab[seq_len(top_n), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

make_sankey_network_data <- function(sankey_df) {
  if (is.null(sankey_df) || nrow(sankey_df) == 0) return(NULL)
  
  sankey_df$A <- as.character(sankey_df$A)
  sankey_df$B <- as.character(sankey_df$B)
  sankey_df$Freq <- as.numeric(sankey_df$Freq)
  
  left_nodes  <- unique(sankey_df$A)
  right_nodes <- unique(sankey_df$B)
  
  nodes <- data.frame(
    name = c(paste0("Left: ", left_nodes), paste0("Right: ", right_nodes)),
    stringsAsFactors = FALSE
  )
  
  source_lookup <- setNames(seq_along(left_nodes) - 1L, left_nodes)
  target_lookup <- setNames(length(left_nodes) + seq_along(right_nodes) - 1L, right_nodes)
  
  links <- data.frame(
    source = as.integer(unname(source_lookup[sankey_df$A])),
    target = as.integer(unname(target_lookup[sankey_df$B])),
    value  = as.numeric(sankey_df$Freq)
  )
  
  list(nodes = nodes, links = links)
}

make_confusion_objects <- function(meta_df, truth_col, pred_col, mode = "Row-normalized") {
  req(truth_col %in% colnames(meta_df), pred_col %in% colnames(meta_df))
  
  df <- meta_df[, c(truth_col, pred_col), drop = FALSE]
  colnames(df) <- c("truth", "pred")
  
  df$truth <- safe_char(df$truth)
  df$pred  <- safe_char(df$pred)
  df <- df[!is.na(df$truth) & !is.na(df$pred) & df$truth != "" & df$pred != "", , drop = FALSE]
  
  if (nrow(df) == 0) {
    return(list(mat = NULL, df = NULL))
  }
  
  conf_counts <- table(df$truth, df$pred)
  
  conf_mat <- switch(
    mode,
    "Counts" = conf_counts,
    "Row-normalized" = prop.table(conf_counts, margin = 1),
    "Column-normalized" = prop.table(conf_counts, margin = 2)
  )
  
  conf_df <- as.data.frame(as.table(conf_mat))
  colnames(conf_df) <- c("TrueLabel", "PredLabel", "Value")
  
  list(
    mat = conf_mat,
    counts = conf_counts,
    df = conf_df,
    raw_df = df
  )
}

make_confusion_heatmap_plot <- function(conf_df, title_text, mode = "Row-normalized") {
  if (is.null(conf_df) || nrow(conf_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No valid data available") +
        theme_void()
    )
  }
  
  label_fmt <- if (mode == "Counts") "%.0f" else "%.2f"
  
  ggplot(conf_df, aes(x = PredLabel, y = TrueLabel, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf(label_fmt, Value)), size = 3) +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(
      title = title_text,
      x = "Predicted cell type",
      y = "True cell type",
      fill = if (mode == "Counts") "Count" else "Value"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

align_confusion_matrices <- function(mat1, mat2) {
  if (is.null(mat1) || is.null(mat2)) return(NULL)
  
  all_rows <- union(rownames(mat1), rownames(mat2))
  all_cols <- union(colnames(mat1), colnames(mat2))
  
  m1 <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
               dimnames = list(all_rows, all_cols))
  m2 <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
               dimnames = list(all_rows, all_cols))
  
  m1[rownames(mat1), colnames(mat1)] <- mat1
  m2[rownames(mat2), colnames(mat2)] <- mat2
  
  list(m1 = m1, m2 = m2)
}

make_difference_heatmap_plot <- function(mat1, mat2, title_text = "Difference heatmap") {
  aligned <- align_confusion_matrices(mat1, mat2)
  
  if (is.null(aligned)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No difference data available") +
        theme_void()
    )
  }
  
  diff_mat <- aligned$m2 - aligned$m1
  diff_df <- as.data.frame(as.table(diff_mat))
  colnames(diff_df) <- c("TrueLabel", "PredLabel", "Value")
  
  ggplot(diff_df, aes(x = PredLabel, y = TrueLabel, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Value)), size = 3) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
    labs(
      title = title_text,
      x = "Predicted cell type",
      y = "True cell type",
      fill = "Delta"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

compute_metrics_table <- function(meta_df, truth_col, pred_col) {
  df <- meta_df[, c(truth_col, pred_col), drop = FALSE]
  colnames(df) <- c("truth", "pred")
  df$truth <- safe_char(df$truth)
  df$pred  <- safe_char(df$pred)
  df <- df[!is.na(df$truth) & !is.na(df$pred) & df$truth != "" & df$pred != "", , drop = FALSE]
  
  if (nrow(df) == 0) return(data.frame())
  
  labels <- sort(unique(c(df$truth, df$pred)))
  
  out <- lapply(labels, function(lbl) {
    tp <- sum(df$truth == lbl & df$pred == lbl)
    fn <- sum(df$truth == lbl & df$pred != lbl)
    fp <- sum(df$truth != lbl & df$pred == lbl)
    
    recall <- if ((tp + fn) == 0) NA else tp / (tp + fn)
    precision <- if ((tp + fp) == 0) NA else tp / (tp + fp)
    f1 <- if (is.na(recall) || is.na(precision) || (recall + precision) == 0) NA else 2 * recall * precision / (recall + precision)
    support <- sum(df$truth == lbl)
    
    data.frame(
      CellType = lbl,
      Support = support,
      Recall = recall,
      Precision = precision,
      F1 = f1,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(out)
}

compute_overall_metrics <- function(meta_df, truth_col, pred_col) {
  df <- meta_df[, c(truth_col, pred_col), drop = FALSE]
  colnames(df) <- c("truth", "pred")
  df$truth <- safe_char(df$truth)
  df$pred  <- safe_char(df$pred)
  df <- df[!is.na(df$truth) & !is.na(df$pred) & df$truth != "" & df$pred != "", , drop = FALSE]
  
  if (nrow(df) == 0) return(data.frame())
  
  class_metrics <- compute_metrics_table(meta_df, truth_col, pred_col)
  overall_acc <- mean(df$truth == df$pred)
  macro_f1 <- mean(class_metrics$F1, na.rm = TRUE)
  balanced_acc <- mean(class_metrics$Recall, na.rm = TRUE)
  
  data.frame(
    Metric = c("Overall accuracy", "Balanced accuracy", "Macro F1"),
    Value = c(overall_acc, balanced_acc, macro_f1)
  )
}

top_misclassifications <- function(meta_df, truth_col, pred_col, top_n = 15, normalize = TRUE) {
  df <- meta_df[, c(truth_col, pred_col), drop = FALSE]
  colnames(df) <- c("truth", "pred")
  df$truth <- safe_char(df$truth)
  df$pred  <- safe_char(df$pred)
  df <- df[!is.na(df$truth) & !is.na(df$pred) & df$truth != "" & df$pred != "", , drop = FALSE]
  
  if (nrow(df) == 0) return(data.frame())
  
  tab <- as.data.frame(table(df$truth, df$pred), stringsAsFactors = FALSE)
  colnames(tab) <- c("TrueLabel", "PredLabel", "Count")
  tab <- tab[tab$TrueLabel != tab$PredLabel, , drop = FALSE]
  
  if (normalize) {
    truth_sizes <- as.data.frame(table(df$truth), stringsAsFactors = FALSE)
    colnames(truth_sizes) <- c("TrueLabel", "TruthTotal")
    tab <- left_join(tab, truth_sizes, by = "TrueLabel")
    tab$Fraction <- tab$Count / tab$TruthTotal
    tab <- tab[order(-tab$Fraction, -tab$Count), ]
  } else {
    tab <- tab[order(-tab$Count), ]
  }
  
  head(tab, top_n)
}

save_plot_file <- function(plot_obj, file, width = 9, height = 7) {
  ext <- tools::file_ext(file)
  if (tolower(ext) == "pdf") {
    ggsave(file, plot = plot_obj, width = width, height = height, device = grDevices::cairo_pdf)
  } else {
    ggsave(file, plot = plot_obj, width = width, height = height, dpi = 300)
  }
}

# =========================
# UI
# =========================

ui <- navbarPage(
  "scRepresenter Explorer",
  
  tabPanel(
    "Input",
    fluidPage(
      br(),
      fluidRow(
        column(
          8,
          wellPanel(
            h3("scRepresenter Explorer"),
            tags$p(
              "This application is designed to explore and compare cell representations, predicted labels, and annotation agreement across different embedding and classification methods."
            ),
            tags$p(
              "Upload a ",
              tags$strong(".h5ad"),
              " object containing cell metadata and low-dimensional embeddings. After loading, you can compare embeddings, inspect classification heatmaps, review performance summaries, and explore label transitions."
            ),
            tags$ul(
              tags$li("Embeddings: visual comparison of multiple low-dimensional representations."),
              tags$li("Classification Heatmaps: compare true cell labels with predictions from different methods."),
              tags$li("Performance Summary: review overall metrics, per-cell-type performance, and top misclassifications."),
              tags$li("Label Flow: visualize how annotations transition between two metadata columns.")
            ),
            fileInput("h5ad_file", "Upload h5ad object", accept = c(".h5ad")),
            actionButton("load_obj", "Load object")
          )
        ),
        column(
          4,
          wellPanel(
            h4("Expected input"),
            tags$p("Your .h5ad object should ideally contain:"),
            tags$ul(
              tags$li("Cell metadata in obs, such as true labels and prediction columns."),
              tags$li("Low-dimensional embeddings in obsm."),
              tags$li("Prediction columns like scgpt_prediction, avg_prediction, conc_prediction, or similar.")
            ),
            tags$p(
              tags$strong("Example embeddings: "),
              "FM, Knowledge_GOcomm, FM_GOcomm_Avg, FM_GOcomm_Con"
            )
          )
        )
      ),
      br(),
      conditionalPanel(
        condition = "output.object_loaded == true",
        fluidRow(
          column(
            6,
            h4("Object summary"),
            verbatimTextOutput("object_summary")
          )
        ),
        br(),
        fluidRow(
          column(6, h4("Available embeddings"), DTOutput("available_embeddings")),
          column(6, h4("Available metadata columns"), DTOutput("available_metadata"))
        )
      )
    )
  ),
  
  tabPanel(
    "Explore",
    fluidPage(
      br(),
      
      fluidRow(
        column(
          12,
          wellPanel(
            h4("Global settings"),
            fluidRow(
              column(3, selectInput("downsample_mode", "Downsampling", choices = c("All", "Random", "Stratified"), selected = "All")),
              column(3, numericInput("downsample_n", "Cells to display", value = 3000, min = 100, step = 100)),
              column(3, selectInput("downsample_group", "Stratify by", choices = NULL)),
              column(3, sliderInput("global_alpha", "Point alpha", min = 0.1, max = 1, value = 0.85, step = 0.05))
            )
          )
        )
      ),
      
      br(),
      
      tabsetPanel(
        id = "explore_tabs",
        
        tabPanel(
          "Embeddings",
          br(),
          fluidRow(
            column(
              6,
              wellPanel(
                h4("Embedding panel 1"),
                selectInput("umap1_reduction", "Reduction", choices = NULL),
                selectInput("umap1_color", "Color by", choices = NULL),
                sliderInput("umap1_pt_size", "Point size", min = 0.2, max = 3, value = 0.8, step = 0.1),
                fluidRow(
                  column(6, downloadButton("download_umap1_png", "Export PNG")),
                  column(6, downloadButton("download_umap1_pdf", "Export PDF"))
                ),
                br(),
                plotlyOutput("umap1_plot", height = "450px")
              )
            ),
            column(
              6,
              wellPanel(
                h4("Embedding panel 2"),
                selectInput("umap2_reduction", "Reduction", choices = NULL),
                selectInput("umap2_color", "Color by", choices = NULL),
                sliderInput("umap2_pt_size", "Point size", min = 0.2, max = 3, value = 0.8, step = 0.1),
                fluidRow(
                  column(6, downloadButton("download_umap2_png", "Export PNG")),
                  column(6, downloadButton("download_umap2_pdf", "Export PDF"))
                ),
                br(),
                plotlyOutput("umap2_plot", height = "450px")
              )
            )
          )
        ),
        
        tabPanel(
          "Classification Heatmaps",
          br(),
          fluidRow(
            column(
              12,
              wellPanel(
                h4("Classification heatmap settings"),
                fluidRow(
                  column(3, selectInput("cm_truth", "True labels", choices = NULL)),
                  column(3, selectInput("cm_pred1", "Method 1", choices = NULL)),
                  column(3, selectInput("cm_pred2", "Method 2", choices = NULL)),
                  column(3, selectInput("cm_mode", "Display mode", choices = c("Counts", "Row-normalized", "Column-normalized"), selected = "Row-normalized"))
                )
              )
            )
          ),
          br(),
          fluidRow(
            column(
              6,
              wellPanel(
                h4("Classification heatmap 1"),
                fluidRow(
                  column(6, downloadButton("download_cm1_png", "Export PNG")),
                  column(6, downloadButton("download_cm1_pdf", "Export PDF"))
                ),
                br(),
                plotOutput("cm1_plot", height = "500px")
              )
            ),
            column(
              6,
              wellPanel(
                h4("Classification heatmap 2"),
                fluidRow(
                  column(6, downloadButton("download_cm2_png", "Export PNG")),
                  column(6, downloadButton("download_cm2_pdf", "Export PDF"))
                ),
                br(),
                plotOutput("cm2_plot", height = "500px")
              )
            )
          ),
          br(),
          fluidRow(
            column(
              12,
              wellPanel(
                h4("Method difference heatmap (Method 2 - Method 1)"),
                fluidRow(
                  column(3, downloadButton("download_cmdiff_png", "Export PNG")),
                  column(3, downloadButton("download_cmdiff_pdf", "Export PDF"))
                ),
                br(),
                plotOutput("cmdiff_plot", height = "550px")
              )
            )
          ),
          br(),
          fluidRow(
            column(
              12,
              wellPanel(
                h3("Performance Summary"),
                fluidRow(
                  column(
                    6,
                    wellPanel(
                      h4("Overall metrics"),
                      DTOutput("overall_metrics_table")
                    )
                  ),
                  column(
                    6,
                    wellPanel(
                      h4("Per-cell-type metrics"),
                      DTOutput("class_metrics_table")
                    )
                  )
                ),
                br(),
                fluidRow(
                  column(
                    12,
                    wellPanel(
                      h4("Top misclassifications"),
                      DTOutput("top_misclass_table")
                    )
                  )
                )
              )
            )
          )
        ),
        
        tabPanel(
          "Label Flow",
          br(),
          fluidRow(
            column(
              12,
              wellPanel(
                h4("Label flow settings"),
                fluidRow(
                  column(2, selectInput("sankey_label1", "Left labels", choices = NULL)),
                  column(2, selectInput("sankey_label2", "Right labels", choices = NULL)),
                  column(2, selectInput("sankey_filter_col", "Filter column", choices = NULL)),
                  column(2, uiOutput("sankey_filter_values_ui")),
                  column(2, numericInput("sankey_min_freq", "Min flow", value = 5, min = 1, step = 1)),
                  column(2, numericInput("sankey_top_n", "Top N flows", value = 50, min = 5, step = 5))
                ),
                br(),
                fluidRow(
                  column(3, downloadButton("download_sankey_html", "Export Sankey HTML")),
                  column(3, downloadButton("download_sankey_csv", "Export Sankey CSV"))
                ),
                br(),
                sankeyNetworkOutput("sankey_plot", height = "500px")
              )
            )
          )
        ),
        
        tabPanel(
          "Help / Interpretation",
          br(),
          fluidRow(
            column(
              12,
              wellPanel(
                h4("How to use this page"),
                tags$p(
                  tags$strong("Embeddings: "),
                  "Compare low-dimensional representations across methods. Use the same metadata column in both panels to visually assess whether cell groups remain compact and well separated."
                ),
                tags$p(
                  tags$strong("Classification Heatmaps: "),
                  "These compare true labels against predicted labels. A strong diagonal indicates better agreement. The method difference heatmap highlights where Method 2 improves or worsens relative to Method 1."
                ),
                tags$p(
                  tags$strong("Performance Summary: "),
                  "Overall metrics give a quick summary, while per-cell-type metrics and top misclassifications identify which cell types are easiest or hardest to predict."
                ),
                tags$p(
                  tags$strong("Label Flow: "),
                  "This shows how labels correspond between two metadata columns. It is useful for exploring agreement, relabelling, and dominant transitions between annotation schemes."
                )
              )
            )
          )
        )
      )
    )
  )
)

# =========================
# Server
# =========================

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    sce = NULL,
    meta_df = NULL,
    rd_names = NULL,
    meta_cols = NULL,
    grouping_cols = NULL,
    prediction_cols = NULL
  )
  
  output$object_loaded <- reactive({
    !is.null(rv$sce)
  })
  outputOptions(output, "object_loaded", suspendWhenHidden = FALSE)
  
  observeEvent(input$load_obj, {
    req(input$h5ad_file)
    
    withProgress(message = "Loading h5ad object...", value = 0, {
      shiny::incProgress(0.1, detail = "Reading file...")
      obj <- anndataR::read_h5ad(input$h5ad_file$datapath, as = "SingleCellExperiment")
      shiny::incProgress(0.1, detail = "Processing metadata...")
      incProgress(0.5)
      
      validate(
        need(inherits(obj, "SingleCellExperiment"),
             "Uploaded object must be a valid .h5ad file readable as SingleCellExperiment")
      )
      
      rd_names <- reducedDimNames(obj)
      validate(
        need(length(rd_names) > 0, "No reduced dimensions found in uploaded object")
      )
      
      meta_df <- get_meta_df(obj)
      meta_cols <- colnames(meta_df)
      grouping_cols <- get_valid_grouping_cols(meta_df)
      prediction_cols <- detect_prediction_cols(meta_df)
      
      validate(
        need(length(grouping_cols) > 0, "No valid grouping metadata columns found")
      )
      
      rv$sce <- obj
      rv$meta_df <- meta_df
      rv$meta_cols <- meta_cols
      rv$rd_names <- rd_names
      rv$grouping_cols <- grouping_cols
      rv$prediction_cols <- prediction_cols
      
      default_meta <- if ("celltype" %in% meta_cols) "celltype" else meta_cols[1]
      default_group <- if ("celltype" %in% grouping_cols) "celltype" else grouping_cols[1]
      pred1_default <- if ("conc_prediction" %in% prediction_cols) "conc_prediction" else if (length(prediction_cols) >= 1) prediction_cols[1] else default_group
      pred2_default <- if ("avg_prediction" %in% prediction_cols) "avg_prediction" else if (length(prediction_cols) >= 2) prediction_cols[2] else pred1_default
      
      updateSelectInput(session, "downsample_group", choices = grouping_cols, selected = default_group)
      
      updateSelectInput(session, "umap1_reduction", choices = rd_names, selected = rd_names[1])
      updateSelectInput(session, "umap2_reduction", choices = rd_names, selected = rd_names[min(2, length(rd_names))])
      updateSelectInput(session, "umap1_color", choices = meta_cols, selected = default_meta)
      updateSelectInput(session, "umap2_color", choices = meta_cols, selected = default_meta)
      
      updateSelectInput(session, "cm_truth", choices = grouping_cols, selected = default_group)
      updateSelectInput(session, "cm_pred1", choices = grouping_cols, selected = pred1_default)
      updateSelectInput(session, "cm_pred2", choices = grouping_cols, selected = pred2_default)
      
      updateSelectInput(session, "sankey_label1", choices = grouping_cols, selected = default_group)
      updateSelectInput(session, "sankey_label2", choices = grouping_cols, selected = pred1_default)
      updateSelectInput(session, "sankey_filter_col", choices = c("None", grouping_cols), selected = "None")
      
      incProgress(0.3)
    })
  })
  
  output$object_summary <- renderText({
    req(rv$sce)
    paste0(
      "Cells: ", ncol(rv$sce), "\n",
      "Features: ", nrow(rv$sce), "\n",
      "Embeddings: ", paste(rv$rd_names, collapse = ", "), "\n",
      "Metadata columns: ", paste(rv$meta_cols, collapse = ", "), "\n",
      "Grouping columns: ", paste(rv$grouping_cols, collapse = ", "), "\n",
      "Prediction columns: ", paste(rv$prediction_cols, collapse = ", ")
    )
  })
  
  output$available_embeddings <- renderDT({
    req(rv$rd_names)
    datatable(data.frame(Embeddings = rv$rd_names), options = list(dom = "t"))
  })
  
  output$available_metadata <- renderDT({
    req(rv$meta_cols)
    datatable(data.frame(Metadata = rv$meta_cols), options = list(dom = "t"))
  })
  
  output$sankey_filter_values_ui <- renderUI({
    req(rv$meta_df, input$sankey_filter_col)
    if (input$sankey_filter_col == "None") return(NULL)
    
    vals <- unique(safe_char(rv$meta_df[[input$sankey_filter_col]]))
    vals <- vals[!is.na(vals) & vals != ""]
    vals <- sort(vals)
    
    selectizeInput(
      "sankey_filter_values",
      "Filter values",
      choices = vals,
      selected = NULL,
      multiple = TRUE
    )
  })
  
  emb1_df <- reactive({
    req(rv$sce, input$umap1_reduction, input$umap1_color)
    df <- make_embedding_df(rv$sce, input$umap1_reduction, input$umap1_color)
    downsample_cells(df, input$downsample_mode, input$downsample_n, input$downsample_group)
  })
  
  emb2_df <- reactive({
    req(rv$sce, input$umap2_reduction, input$umap2_color)
    df <- make_embedding_df(rv$sce, input$umap2_reduction, input$umap2_color)
    downsample_cells(df, input$downsample_mode, input$downsample_n, input$downsample_group)
  })
  
  emb1_plot_obj <- reactive({
    make_embedding_plot(
      emb1_df(),
      paste0("Embedding: ", input$umap1_reduction),
      input$umap1_color,
      input$umap1_pt_size,
      input$global_alpha
    )
  })
  
  emb2_plot_obj <- reactive({
    make_embedding_plot(
      emb2_df(),
      paste0("Embedding: ", input$umap2_reduction),
      input$umap2_color,
      input$umap2_pt_size,
      input$global_alpha
    )
  })
  
  output$umap1_plot <- renderPlotly({
    ggplotly(emb1_plot_obj(), tooltip = "text")
  })
  
  output$umap2_plot <- renderPlotly({
    ggplotly(emb2_plot_obj(), tooltip = "text")
  })
  
  cm1_obj <- reactive({
    req(rv$meta_df, input$cm_truth, input$cm_pred1, input$cm_mode)
    make_confusion_objects(rv$meta_df, input$cm_truth, input$cm_pred1, input$cm_mode)
  })
  
  cm2_obj <- reactive({
    req(rv$meta_df, input$cm_truth, input$cm_pred2, input$cm_mode)
    make_confusion_objects(rv$meta_df, input$cm_truth, input$cm_pred2, input$cm_mode)
  })
  
  cm1_plot_obj <- reactive({
    make_confusion_heatmap_plot(
      cm1_obj()$df,
      paste0("Classification heatmap: ", input$cm_pred1),
      input$cm_mode
    )
  })
  
  cm2_plot_obj <- reactive({
    make_confusion_heatmap_plot(
      cm2_obj()$df,
      paste0("Classification heatmap: ", input$cm_pred2),
      input$cm_mode
    )
  })
  
  diff_plot_obj <- reactive({
    make_difference_heatmap_plot(
      cm1_obj()$mat,
      cm2_obj()$mat,
      title_text = paste0(input$cm_pred2, " - ", input$cm_pred1)
    )
  })
  
  output$cm1_plot <- renderPlot({
    print(cm1_plot_obj())
  })
  
  output$cm2_plot <- renderPlot({
    print(cm2_plot_obj())
  })
  
  output$cmdiff_plot <- renderPlot({
    print(diff_plot_obj())
  })
  
  output$overall_metrics_table <- renderDT({
    req(rv$meta_df, input$cm_truth, input$cm_pred1)
    datatable(
      compute_overall_metrics(rv$meta_df, input$cm_truth, input$cm_pred1),
      options = list(dom = "t", pageLength = 5)
    )
  })
  
  output$class_metrics_table <- renderDT({
    req(rv$meta_df, input$cm_truth, input$cm_pred1)
    datatable(
      compute_metrics_table(rv$meta_df, input$cm_truth, input$cm_pred1),
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  output$top_misclass_table <- renderDT({
    req(rv$meta_df, input$cm_truth, input$cm_pred1)
    datatable(
      top_misclassifications(rv$meta_df, input$cm_truth, input$cm_pred1, top_n = 15, normalize = TRUE),
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  sankey_debug_df <- reactive({
    req(rv$meta_df, input$sankey_label1, input$sankey_label2)
    make_sankey_debug_df(
      meta_df = rv$meta_df,
      labelA = input$sankey_label1,
      labelB = input$sankey_label2,
      filter_col = input$sankey_filter_col,
      filter_vals = input$sankey_filter_values,
      min_freq = input$sankey_min_freq,
      top_n = input$sankey_top_n
    )
  })
  
  sankey_widget <- reactive({
    req(sankey_debug_df())
    sank <- make_sankey_network_data(sankey_debug_df())
    validate(need(!is.null(sank), "No Sankey data available for the selected inputs."))
    
    sankeyNetwork(
      Links = sank$links,
      Nodes = sank$nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      fontSize = 12,
      nodeWidth = 30,
      sinksRight = FALSE
    )
  })
  
  output$sankey_plot <- renderSankeyNetwork({
    sankey_widget()
  })
  
  # =========================
  # Downloads
  # =========================
  
  output$download_umap1_png <- downloadHandler(
    filename = function() "embedding_panel1.png",
    content = function(file) save_plot_file(emb1_plot_obj(), file)
  )
  
  output$download_umap1_pdf <- downloadHandler(
    filename = function() "embedding_panel1.pdf",
    content = function(file) save_plot_file(emb1_plot_obj(), file)
  )
  
  output$download_umap2_png <- downloadHandler(
    filename = function() "embedding_panel2.png",
    content = function(file) save_plot_file(emb2_plot_obj(), file)
  )
  
  output$download_umap2_pdf <- downloadHandler(
    filename = function() "embedding_panel2.pdf",
    content = function(file) save_plot_file(emb2_plot_obj(), file)
  )
  
  output$download_cm1_png <- downloadHandler(
    filename = function() "classification_heatmap_1.png",
    content = function(file) save_plot_file(cm1_plot_obj(), file, width = 10, height = 8)
  )
  
  output$download_cm1_pdf <- downloadHandler(
    filename = function() "classification_heatmap_1.pdf",
    content = function(file) save_plot_file(cm1_plot_obj(), file, width = 10, height = 8)
  )
  
  output$download_cm2_png <- downloadHandler(
    filename = function() "classification_heatmap_2.png",
    content = function(file) save_plot_file(cm2_plot_obj(), file, width = 10, height = 8)
  )
  
  output$download_cm2_pdf <- downloadHandler(
    filename = function() "classification_heatmap_2.pdf",
    content = function(file) save_plot_file(cm2_plot_obj(), file, width = 10, height = 8)
  )
  
  output$download_cmdiff_png <- downloadHandler(
    filename = function() "method_difference_heatmap.png",
    content = function(file) save_plot_file(diff_plot_obj(), file, width = 10, height = 8)
  )
  
  output$download_cmdiff_pdf <- downloadHandler(
    filename = function() "method_difference_heatmap.pdf",
    content = function(file) save_plot_file(diff_plot_obj(), file, width = 10, height = 8)
  )
  
  output$download_sankey_html <- downloadHandler(
    filename = function() "label_flow.html",
    content = function(file) saveWidget(sankey_widget(), file = file, selfcontained = TRUE)
  )
  
  output$download_sankey_csv <- downloadHandler(
    filename = function() "label_flow.csv",
    content = function(file) write.csv(sankey_debug_df(), file, row.names = FALSE)
  )
}

shinyApp(ui, server, options = list(port = 3838))