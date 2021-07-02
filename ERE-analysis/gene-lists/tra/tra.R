#' Binarize expression values based on max gap between values
#' @export
find_tra_tissues <- function(tissue_expression) {
    nrtissues <- length(tissue_expression)
    expression_sort <- sort(tissue_expression)
    diff_expression <- diff(expression_sort)
    max_gap <- which.max(diff_expression) + 1
    thr_exp <- as.numeric(expression_sort[max_gap])
    sapply(tissue_expression, function(x) ifelse(x < thr_exp, 0, 1))
}

#' Maximum expression across genes adjusting for NA

#' fmax takes a vector of gene expression for a single gene across tissues;
#' tissue with unknown expression for this gene indicated as NA. If all values
#' are NA, NA is returned, else NAs are ignored and max estimated on available
#' expression values.
#' Function from Mostacci et al (2016) Bioinformatics, Supplement.
#'
#' @param tissue_expression [numeric vector] Vector of gene expression
#' measurements for a single gene across length(tissue_expression) tissues
#' @return maximum expression value across tissues [double]
#'
#' @export
fmax <- function(tissue_expression) {
    if (!all(is.na(tissue_expression))) {
        res <- max(tissue_expression, na.rm = TRUE)
    } else {
        res <- NA
    }
    return(res)
}

#' Calculate tissue specificity measure tau
#'
#'
#' Tau takes a vector of gene expression for a single gene across tissues;
#' tissue with unknown expression for this gene indicated as NA; normalisation
#' of tau to number of tissues - 1 is computed across all tissues, ie taking
#' NA tissues into account.
#' Function adapted from Mostacci et al (2016) Bioinformatics, Supplement.
#'
#' @param tissue_expression [numeric vector] Vector of gene expression
#' measurements for a single gene across length(tissue_expression) tissues
#' @return tau [double]
#' @export
tau <- function(tissue_expression) {
    if (all(!is.na(tissue_expression))) {
        nrtissue <- length(tissue_expression)
        tissue_expression <- tissue_expression[!is.na(tissue_expression)]
        if (min(tissue_expression) >= 0) {
            if (max(tissue_expression) != 0) {
                tissue_expression <- tissue_expression / max(tissue_expression)
                res <- sum((1-tissue_expression))
                res <- res / (nrtissue - 1)
            } else {
                res <- 0
            }
        } else {
            res <- NA
        }
    } else {
        res <- NA
    }
    return(res)
}
