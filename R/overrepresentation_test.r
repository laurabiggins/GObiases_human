#' overrep_test
#'
#' Overrepresentation test - Functional (usually gene ontology) analysis
#'
#' @param categories list of named character vectors containing the functional groups.
#' Each vector should contain gene names or IDs. The name of each vector should
#' be the functional category.
#' @param query_genes character vector of gene names
#' @param bg_genes character vector of gene names to use as the background set
#' containing genes and associated information including location of genes
#' @param min_query minimum number of query genes in the category for it be tested,
#' no point having a category with one query gene in it
#' @param pval_threshold p value threshold. Only results with p-value/corrected p-value
#' less than this thrreshold will be returned.
#' @param mult_test apply multiple testing correction (Benjamini-Hochberg FDR is used)
#' and use the corrected value to filter for significant results.
#' This should usually be set to TRUE (default). If set to false, the correction is
#' still applied but the uncorrected pvalue is used to filter by.
#' @param super_strict stricter pvalue correction where it takes the number of
#' tests as being the total number of functional categories. By default the number
#' of tests corrected for is only the number of functional categories that contain
#' > min_query genes
#' @return results of functional overrepresentation test. If no categories have a
#' p-value <= pval_threshold a NULL object will be returned.
#' @examples
#' go_results <- overrep_test(all_go_categories, genes_1, bg_genes)
#' head(go_results)


# list of categories is all the functional categories from the gmt file
overrep_test <- function(categories, query_genes, bg_genes, min_query = 3,
                         pval_threshold = 0.05, ease = TRUE, sig_digits = 4,
                         mult_test = TRUE, super_strict = FALSE,
                         unique_genes_in_categories = NULL) {
  #  browser()
  query_genes <- clean_text(query_genes)
  bg_genes    <- clean_text(bg_genes)

  if(is.null(unique_genes_in_categories)){
      unique_genes_in_categories <- unique(as.vector(unlist(categories)))
  }

  query_genes <- query_genes[query_genes %in% unique_genes_in_categories]
  #bg_genes    <- bg_genes[bg_genes %in% unique_genes_in_categories]

  matched_categories <- categories[sapply(categories, function(x) {
    sum(!is.na(fastmatch::fmatch(query_genes, x))) >= min_query
  })]

  if (length(matched_categories) < 1) {
    warning("The set of query genes were not found in the functional categories")
    return(NULL)
  }

  df <- data.frame(
    query_in_category = sapply(matched_categories, function(x) {
      sum(query_genes %in% x)
    }),
    bg_in_category = sapply(matched_categories, function(x) {
      sum(!is.na(fastmatch::fmatch(x, bg_genes)))
    }),
    category_length = sapply(matched_categories, length)
  )

  df$enrichment <- (df$query_in_category/length(query_genes))/(df$bg_in_category/length(bg_genes))

  query_not_in_category <- length(query_genes) - df$query_in_category
  bg_not_in_category    <- length(bg_genes)    - df$bg_in_category

  ifelse(
    ease == TRUE,
    query_count <- df$query_in_category - 1,
    query_count <- df$query_in_category
  )

  contingency_values <- cbind(
    query_count,
    df$bg_in_category,
    query_not_in_category,
    bg_not_in_category
  )

  df$pval <- apply(
    X = contingency_values,
    MARGIN = 1,
    FUN = function(x) {
      fisher.test(matrix(x,nrow = 2), alternative = "greater")$p.value
    }
  )

  ifelse(
    super_strict,
    n_tests <- length(categories),
    n_tests <- length(matched_categories)
  )

  df$adj_pval <- p.adjust(df$pval, method = "BH", n = n_tests)

  if (mult_test) {

    if (sum(df$adj_pval <= pval_threshold) == 0) {
      return(NULL)
    } else {

      df <- df[df$adj_pval <= pval_threshold, ]
    }
  }else{

    if (sum(df$pval <= pval_threshold) == 0) {
      return(NULL)
    }else {
      df <- df[df$pval <= pval_threshold, ]
    }
  }

  df[, 5:6] <- apply(df[, 5:6], 2, signif, digits = sig_digits)

  df[order(df$adj_pval), ]

}

#' Clean text
#'
#' Remove erroneous text characters.
#'
#' @param text character vector
#' @param chars character vector of characters to remove from text
#' @param remove_empty remove empty elements from the vector.
#' @param remove_dup remove duplicated elements from the vector.
#' @param upper_case convert to upper case
#' @return A character vector
#' @examples
#' clean_text("\thello_world")
#' clean_text("\thello_world", upper_case = FALSE)

clean_text <- function(text,
                       chars = c("\t", "\r", "\n", "\\", "\ ", "\""),
                       remove_empty = TRUE,
                       remove_dup = TRUE,
                       upper_case = TRUE
) {
    char_pattern <- paste(chars, collapse = "|")
    char_pattern <- paste0("[", char_pattern, "]")
    cleaned <- gsub(pattern = char_pattern, replacement = "", x = as.character(text))
    if (remove_dup) cleaned <- remove_duplicates(cleaned)
    if (upper_case) cleaned <- toupper(cleaned)
    ifelse(remove_empty, return(cleaned[sapply(cleaned, nchar) > 0]), return(cleaned))
}


#' remove_duplicates
#'
#' remove duplicates from vector or data frame
#'
#' @keywords internal
#'
#' @param x vector or data frame to remove duplicates from
#' @param column If a data frame is passed in, the column that is used
#' for deduplicating (default = 1)
#'
#' @return A character vector or data frame.
#' @examples
#' hello_string <- c("h", "e", "l", "l", "o")
#' hello_string
#' remove_duplicates(hello_string)
#'
#' df <- data.frame(hello_string, 1:5)
#' df
#' remove_duplicates(df)
remove_duplicates <- function(x, column = 1) {
    if (is.data.frame(x) | is.matrix(x)) {
        return(x[!duplicated(x[,column]), ])
    } else if (is.list(x)) {
        stop("x must be a vector, data frame or matrix")
    } else if (is.vector(x)) {
        return(unique(x))
    } else{
        warning("data must be a vector, data frame or matrix")
    }
}
