# Modified scripts from lefser package to return p-values

## Kruskal-Wallis Rank Sum Test for the classes
filterKruskal <- function(relab, class, p.value, method = method) {
  # applies "kruskal.test.alt" function to each row (feature) of relab
  # to detect differential abundance between classes, 0 and 1
  kw.res <- apply(relab, 1L, function(x) {
    kruskal.test(x ~ class)[["p.value"]]
  })
  # TRUE for p-values less than or equal to kw.threshold
  kw.res <- stats::p.adjust(kw.res, method = method)
  kw.sub <- kw.res < p.value
  
  # NAs are FALSE
  kw.sub[is.na(kw.sub)] <- FALSE
  
  kw.res <<- kw.res[kw.sub, drop = FALSE]
  # extracts features with statistically significant differential abundance
  # from "relab" matrix
  relab[kw.sub,, drop = FALSE]
}

my_dropFeatures <- function(x) {
  se <- x
  row_data <- as.data.frame(SummarizedExperiment::rowData(se))
  nrow1 <- nrow(row_data)
  for (i in seq_along(row_data)) {
    sumNA <- sum(is.na(row_data[[i]]))
    if (sumNA > 0) {
      message(
        sumNA, " features don't have ", colnames(row_data)[i],
        " information."
      )
    }
  }
  # row_data <- tidyr::drop_na(row_data)
  nrow2 <- nrow(row_data)
  if (nrow1 > nrow2) {
    message(
      "Dropped ", nrow1 - nrow2,
      " features without full taxonomy information."
    )
  }
  se <- se[rownames(row_data),]
  pathStrings <- sort(unique(lefser:::.rowData2PathStrings(se)))
  l1 <- length(pathStrings)
  tips <- stringr::str_extract(pathStrings, "(|\\w__\\w+)?$")
  dupTips <- tips[which(duplicated(tips))]
  pathStrings <- pathStrings[!tips %in% dupTips]
  l2 <- length(pathStrings)
  if (l1 > l2) {
    message(
      "Dropped ", l1 - l2,
      " features with duplicated tip names."
    )
  }
  list(se = se, pathStrings = pathStrings)
}

## lefser function
my_lefser <-  function(relab,
           kruskal.threshold = 0.05,
           wilcox.threshold = 0.05,
           lda.threshold = 2.0,
           classCol = "CLASS",
           subclassCol = NULL,
           assay = 1L,
           trim.names = FALSE,
           checkAbundances = TRUE,
           method = "none",
           ...
  ) {
    relab_data <- assay(relab, i = assay)
    
    ## Check whether relative abundance is provided or not
    if (checkAbundances && !identical(all.equal(colSums(relab_data),
                                                rep(1e6, ncol(relab_data)),
                                                check.attributes = FALSE),
                                      TRUE)) {
      warning("Convert counts to relative abundances with 'relativeAb()'")
    }
    
    ## Extract the class/subclass information
    classf <- colData(relab)[[classCol]]
    classf <- as.factor(classf)
    lclassf <- levels(classf)
    if (is.null(classf) || !identical(length(lclassf), 2L)) {
      msg <- "'classCol' must refer to a valid dichotomous (two-level) variable"
      stop(msg) # ensure the class has only two levels
    }
    message(
      "The outcome variable is specified as '", classCol,
      "' and the reference category is '", lclassf[1],
      "'.\n See `?factor` or `?relevel` to change the reference category."
    )
    
    ## Kruskal-Wallis Rank Sum Test for the classes
    relab_sub <- filterKruskal(relab = relab_data,
                               class = classf,
                               p.value = kruskal.threshold,
                               method = method)
    
    ## Wilcoxon Rank-Sum Test for the sub-classes
    if (!is.null(subclassCol)) {
      subclass <- as.factor(colData(relab)[[subclassCol]])
      subclass <- droplevels(subclass)
      ## z-statistics result of the features passing the significance cut-off
      relab_sub <- fillPmatZmat(class = classf,
                                subclass = subclass,
                                relab_sub = relab_sub,
                                p.threshold = wilcox.threshold,
                                method = method)
    }
    
    ## Return an empty data table if there is no significant features
    if(nrow(relab_sub) == 0L){
      return(lefser:::.return_no_results())
    }
    
    ## Transposed relative abundance matrix with the 'class' column
    relab_sub_t <- t(relab_sub)
    relab_sub_t_df <- as.data.frame(relab_sub_t)
    
    # relab_sub_t_df <- createUniqueValues(df = relab_sub_t_df, class = classf)
    relab_sub_t_df <- cbind(relab_sub_t_df, class = classf)
    
    ## LDA model
    warn <- testthat::capture_warnings(
      raw_lda_scores <- lefser:::ldaFunction(relab_sub_t_df, lclassf)
    )
    
    ## Warning collinearity and recommend `get_terminal_nodes`
    if (length(warn) && nzchar(warn)) {
      msg <- "Variables in the input are collinear. Try only with the terminal nodes using `get_terminal_nodes` function"
      warning(msg)
    }
    ## Processing LDA scores
    processed_scores <-
      sign(raw_lda_scores) * log((1 + abs(raw_lda_scores)), 10)
    processed_sorted_scores <- sort(processed_scores)   # sorting of scores
    scores_df <- data.frame(features = names(processed_sorted_scores),
                            scores = as.vector(processed_sorted_scores),
                            stringsAsFactors = FALSE)
    scores_df <- lefser:::.trunc(scores_df, trim.names)   # short-form of taxa name
    
    ## Filter with LDA threshold
    threshold_scores <- abs(scores_df$scores) >= lda.threshold
    res_scores <- scores_df[threshold_scores, , drop = FALSE]
    class(res_scores) <- c("lefser_df", class(res_scores))
    attr(res_scores, "classes") <- lclassf
    if (nrow(res_scores) == 0L) {
      return(lefser:::.return_no_results())
    }
    
    res_scores$kw_p.value = kw.res[res_scores$features]
    
    ## Add attributes with argument values
    ## This is used for plottin functions
    attr(res_scores, "inputSE") <- relab
    attr(res_scores, "kth") <-  kruskal.threshold
    attr(res_scores, "wth") <- wilcox.threshold
    attr(res_scores, "ldath") <- lda.threshold
    attr(res_scores, "class_arg") <- classCol
    attr(res_scores, "subclass_arg") <- subclassCol
    attr(res_scores, "method") <- method
    # attr(res_scores, "lgroupf") <- lgroupf[1]
    # attr(res_scores, "case") <- lgroupf[2]
    
    ## Some more attributes to create the cladogram.
    # pathStrings <- .selectPathStrings(relab, res_scores)
    # # attr(res_scores, "pathStrings") <- pathStrings
    # attr(res_scores, "tree") <- .toTree(pathStrings)
    attr(res_scores, "lclassf") <- lclassf[1]
    attr(res_scores, "case") <- lclassf[2]
    
    res_scores
  }

## lefserClades
my_lefserClades <- function(relab, ...) {
  se <- lefser:::.selectTaxRanks(relab)
  se <- lefser:::.appendRankLetter(se)
  l <- my_dropFeatures(se)
  se <- l[["se"]]
  pathStrings <- l[["pathStrings"]]
  seL <- as.list(mia::splitByRanks(se))
  ## Kingdom would not be informative
  # seL <- seL[!names(seL) %in% "kingdom"]
  msgRanks <- paste(names(seL), collapse = ", ")
  msgRanks <- msgRanks[length(msgRanks):1]
  message(
    "lefser will be run at the ", msgRanks, " level."
  )
  seL <- purrr::imap(seL, function(x, idx) {
    seVar <- x
    row_data <- as.data.frame(SummarizedExperiment::rowData(seVar))
    row_data <- row_data[,1:which(colnames(row_data) == idx), drop = FALSE]
    # row_data <- purrr::discard(row_data, ~ all(is.na(.x)))
    SummarizedExperiment::rowData(seVar) <- S4Vectors::DataFrame(row_data)
    newRowNames <- lefser:::.rowData2PathStrings(seVar)
    BiocGenerics::rownames(seVar) <- newRowNames
    seVar
  })
  resL <- purrr::imap(seL, function(x, idx, ...) {
    message(
      "\n>>>> Running lefser at the ", idx, " level.",
      " <<<<"
    )
    withCallingHandlers(
      my_lefser(relab = x,...),
      warning = function(w) {
        if (grepl("relativeAb", w$message)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  },
  ...
  )
  controlVar <- resL |>
    purrr::map(~ attr(.x, "lclassf")) |>
    unlist(use.names = FALSE) |>
    unique()
  caseVar <- resL |>
    purrr::map(~ attr(.x, "case")) |>
    unlist(use.names = FALSE) |>
    unique()
  classArgVar <- resL |>
    purrr::map(~ attr(.x, "class_arg")) |>
    unlist(use.names = FALSE) |>
    unique()
  subclassArgVar <- resL |>
    purrr::map(~ attr(.x, "subclass_arg")) |>
    unlist(use.names = FALSE) |>
    unique()
  names(resL) <- names(seL)
  res <- dplyr::bind_rows(resL, .id = "Rank") |>
    dplyr::relocate(.data$Rank, .after = tidyselect::last_col())
  class(res) <- c("lefser_df_clades", class(res))
  attr(res, "pathStrings") <- pathStrings
  attr(res, "tree") <- lefser:::.toTree(pathStrings)
  attr(res, "lclassf") <- controlVar
  attr(res, "case") <- caseVar
  attr(res, "inputSE") <- seL
  attr(res, "class_arg") <- classArgVar
  attr(res, "subclass_arg") <- subclassArgVar
  return(res)
}

my_lefserPlot <- function(df,
                          colors = "c",
                          trim.names = TRUE,
                          title = "",
                          label.font.size = 3,
                          other.font.size = 20) {
  
  df <- lefser:::.trunc(df, trim.names)
  lclassf <- attr(df, "lclassf")
  colors <- lefser:::.selectPalette(colors)
  
  ## Create the `class` column
  if (!is.null(lclassf)) {
    class <- ifelse(df$scores > 0,  attr(df, "case"), attr(df, "lclassf"))
    df$class <- factor(class, levels = c(attr(df, "case"), attr(df, "lclassf")))
  } else {
    class <- ifelse(df$scores > 0, 1, 0)
    df$class <- as.factor(class)
  }
  
  ## Add the `order` column based on the scores
  ## To make duplicated features behave independently
  df <- df %>%
    arrange(scores) %>%
    dplyr::mutate(order = row_number())
  # mutate(order = seq_len(nrow(.)))
  
  plt <-
    ggplot(df, aes(factor(order), scores, width = 0.75)) + # Plot same x-axis values separately
    ylab("LDA SCORE (log 10)") +
    ggtitle(title) +
    theme_void() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = other.font.size),
      axis.text.y  = element_blank(),
      axis.text.x  = element_text(vjust = 0.7, size = other.font.size)) +
    geom_point(
      stat = "identity", aes(color = class, size = -log10(kw_p.value))) +
    labs(color = NULL, size = "-log10(p-value)") +
    theme(    # Legends
      legend.position = "top",
      #legend.title = element_blank(),
      legend.text = element_text(size = other.font.size),
      legend.key.height = unit(0.07, 'cm'),
      legend.key.width = unit(0.6, 'cm')) +
    scale_fill_manual(values = colors) +
    geom_text(    # Feature labeling
      aes(y = 0, label = features),
      hjust = ifelse(df$scores < 0, 0, 1),
      nudge_y = ifelse(df$scores < 0, 0.1, -0.1),
      color = "black",
      size = label.font.size) +
    theme(    # Guide lines
      panel.grid.major.x = element_line(
        color = "grey", linewidth = 0.5, linetype = "dotted"),
      panel.grid.minor.x = element_line(
        color = "grey", linewidth = 0.5, linetype = "dotted")) +
    coord_flip()
  
  return(plt)
}


