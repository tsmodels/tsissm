check_weights <- function(w, n)
{
    if (length(w) != n) stop(paste0("\nweights must be a vector of length ", n))
    w <- as.numeric(w)
    return(w)
}

#' Ensembling of Models and Predictions
#'
#' @description Ensembles estimated, predicted and simulated objects arising from the automatic selection method. 
#' @param object an object of class \dQuote{tsissm.selection}, \dQuote{tsissm.selection_predict} or
#' \dQuote{tsissm.selection_simulate}.
#' @param weights a vector of weights equal to the number of models to be ensembled.
#' @param start the index for the state decomposition (when all lambda equal). If 1,
#' will return the predicted states else will return the lagged predicted states (which
#' can be summed to return the predictive distribution).
#' @param ... not used.
#' @returns For the estimated object, a list with the weighted fit and errors, whilst
#' for the predicted and simulated objects a list with the ensemble predictions
#' as a \dQuote{tsmodel.predict} object. If the models were all estimated with
#' the same Box Cox lambda, then the weighted state decomposition is also returned
#' inside the \dQuote{tsmodel.predict} object.
#' Additionally, for the predicted object, the ensemble analytic mean is also returned.
#' @note
#' Only the size of the weights is checked (should be equal to number of selected models),
#' but not checks are performed on the values of the weights or whether they sum to 1. This
#' is left to the user.
#' When lambda is 1, the sum of component will be off by 1 versus the weighted distribution
#' since the Box Cox transform contains an offset (so it is to be expected). To replicate the
#' value of the predicted distribution by summing the decomposed components, the argument
#' \dQuote{start} should be set to 0 to return the lagged predicted components.
#' @aliases tsensemble
#' @method tsensemble tsissm.selection
#' @rdname tsensemble
#' @export
#'
tsensemble.tsissm.selection <- function(object, weights = NULL, start = 1, ...)
{
    weights <- check_weights(weights, n = length(object$models))
    tab <- object$table
    m <- length(object$models)
    if (length(unique(tab$lambda)) == 1) {
        S <- .ensemble_estimated_states(object, weights, start)
        f <- do.call(cbind, lapply(1:m, function(i) coredata(fitted(object$models[[i]]))))
        f <- f %*% weights
        f <- xts(f, object$models[[1]]$spec$target$index)
        colnames(f) <- "Fitted"
        original_series <- xts(object$models[[1]]$spec$target$y_orig, index(f))
        L <- list(fitted = f, original_series = original_series, decomposition = S)
    } else {
        f <- do.call(cbind, lapply(1:m, function(i) coredata(fitted(object$models[[i]]))))
        f <- f %*% weights
        f <- xts(f, object$models[[1]]$spec$target$index)
        colnames(f) <- "Fitted"
        original_series <- xts(object$models[[1]]$spec$target$y_orig, index(f))
        L <- list(fitted = f, original_series = original_series)
    }
    return(L)
}

#' @method tsensemble tsissm.selection_predict
#' @rdname tsensemble
#' @export
#'
tsensemble.tsissm.selection_predict <- function(object, weights = NULL, start = 1, ...)
{
    # check is all lambda is ths same
    # if yes, then ensemble states and decompose   
    # return tsmodel.predict with decomposition
    weights <- check_weights(weights, n = length(object$models))
    nsim <- NROW(object$models[[1]]$distribution)
    h <- NCOL(object$models[[1]]$distribution)
    tab <- object$table
    m <- length(object$models)
    freq <- sort(unique(sapply(object$models, function(x) x$spec$seasonal$seasonal_frequency)))
    if (length(unique(tab$lambda)) == 1) {
        S <- .ensemble_predicted_states(object, weights, start)
        list_of_distributions <- lapply(1:m, function(i) object$models[[i]]$distribution)
        distribution <- Reduce(`+`, Map(`*`, list_of_distributions, weights))
        original_series <- object$models[[1]]$original_series
        analytic_mean <- do.call(cbind, lapply(1:m, function(i) coredata(object$models[[i]]$analytic_mean)), quote = T)
        analytic_mean <- analytic_mean %*% weights
        L <- list(distribution = distribution, original_series = original_series, analytic_mean = analytic_mean, frequency = freq,
                  decomposition = S)
        class(L) <- c("tsissm.predict", "tsmodel.predict")
    } else {
        list_of_distributions <- lapply(1:m, function(i) object$models[[i]]$distribution)
        distribution <- Reduce(`+`, Map(`*`, list_of_distributions, weights))
        original_series <- object$models[[1]]$original_series
        analytic_mean <- do.call(cbind, lapply(1:m, function(i) coredata(object$models[[i]]$analytic_mean)), quote = T)
        analytic_mean <- analytic_mean %*% weights
        L <- list(distribution = distribution, original_series = original_series, analytic_mean = analytic_mean, frequency = freq)
        class(L) <- c("tsissm.predict", "tsmodel.predict")
    }
    return(L)
}

#' @method tsensemble tsissm.selection_simulate
#' @rdname tsensemble
#' @export
#'
tsensemble.tsissm.selection_simulate <- function(object, weights = NULL, start = 1, ...)
{
    # check is all lambda is ths same
    # if yes, then ensemble states and decompose   
    # return tsmodel.predict with decomposition
    weights <- check_weights(weights, n = length(object$models))
    nsim <- NROW(object$models[[1]]$distribution)
    h <- NCOL(object$models[[1]]$distribution)
    tab <- object$table
    m <- length(object$models)
    if (length(unique(tab$lambda)) == 1) {
        S <- .ensemble_simulated_states(object, weights, start)
        list_of_distributions <- lapply(1:m, function(i) object$models[[i]]$distribution)
        distribution <- Reduce(`+`, Map(`*`, list_of_distributions, weights))
        L <- list(distribution = distribution, decomposition = S)
    } else {
        list_of_distributions <- lapply(1:m, function(i) object$models[[i]]$distribution)
        distribution <- Reduce(`+`, Map(`*`, list_of_distributions, weights))
        L <- list(distribution = distribution)
    }
    return(L)
}

.ensemble_decompose <- function(object, start = 1)
{
    n_cores <- nbrOfWorkers()
    if (n_cores <= 1) {
        out <- lapply(1:length(object$models), function(i){
            tsdecompose(object$models[[i]], start = start)
        })
    } else {
        out <- future_lapply(1:length(object$models), function(i){
            tsdecompose(object$models[[i]], start = start)
        }, future.seed = TRUE, future.packages = "tsissm")
        out <- eval(out)
    }
    return(out)
}


.ensemble_estimated_states <- function(object, weights, start = 1)
{
    S <- .ensemble_decompose(object, start = start)
    m <- length(S)
    d_names <- unique(unlist(sapply(1:m, function(i) colnames(S[[i]]))))
    n <- length(d_names)
    L <- vector(mode = "list", length = n)
    zero_series <- coredata(S[[1]][,1] * 0)
    original_index <- index(S[[1]])
    for (j in 1:n) {
        ors <- do.call(cbind, lapply(1:m, function(i) {
            snames <- colnames(S[[i]])
            if (any(grepl(paste0("^",d_names[j]), snames))) {
                ix <- which(grepl(paste0("^",d_names[j]), snames))
                return(coredata(S[[i]][,ix]))
            } else {
                return(zero_series)
            }
        }))
        ors <- xts(ors %*% weights, original_index)
        L[[j]] <- ors
    }
    L <- do.call(cbind, L)
    names(L) <- d_names
    reindex <- which(colnames(L) == "Irregular")
    # place Irregular component last
    L <- cbind(L[,-reindex], L[,reindex])
    return(L)
}

.ensemble_predicted_states <- function(object, weights, start = 1)
{
    S <- .ensemble_decompose(object, start = start)
    m <- length(S)
    d_names <- unique(unlist(sapply(1:m, function(i) names(S[[i]]))))
    n <- length(d_names)
    L <- vector(mode = "list", length = n)
    zero_matrix <- S[[1]][[1]]$distribution * 0
    zero_series <- coredata(S[[1]][[1]]$original_series * 0)
    original_index <- index(S[[1]][[1]]$original_series)
    for (j in 1:n) {
        tmp <- lapply(1:m, function(i) {
            snames <- names(S[[i]])
            if (any(grepl(paste0("^",d_names[j]),snames))) {
                ix <- which(grepl(paste0("^",d_names[j]), snames))
                return(S[[i]][[ix]]$distribution)
            } else {
                return(zero_matrix)
            }
        })
        tmp <- Reduce(`+`, Map(`*`, tmp, weights))
        
        ors <- do.call(cbind, lapply(1:m, function(i) {
            snames <- names(S[[i]])
            if (any(grepl(paste0("^",d_names[j]),snames))) {
                ix <- which(grepl(paste0("^",d_names[j]), snames))
                return(coredata(S[[i]][[ix]]$original_series))
            } else {
                return(zero_series)
            }
        }))
        ors <- xts(ors %*% weights, original_index)
        colnames(ors) <- d_names[j]
        Z <- list(original_series = ors, distribution = tmp)
        class(Z) <- "tsmodel.predict"
        L[[j]] <- Z
    }
    names(L) <- d_names
    return(L)
}

.ensemble_simulated_states <- function(object, weights, start = 1)
{
    S <- .ensemble_decompose(object, start = start)
    m <- length(S)
    d_names <- unique(unlist(sapply(1:m, function(i) names(S[[i]]))))
    n <- length(d_names)
    L <- vector(mode = "list", length = n)
    zero_matrix <- S[[1]][[1]]$distribution * 0
    for (j in 1:n) {
        tmp <- lapply(1:m, function(i) {
            snames <- names(S[[i]])
            if (any(grepl(paste0("^",d_names[j]),snames))) {
                ix <- which(grepl(paste0("^",d_names[j]), snames))
                return(S[[i]][[ix]]$distribution)
            } else {
                return(zero_matrix)
            }
        })
        tmp <- Reduce(`+`, Map(`*`, tmp, weights))
        L[[j]] <- tmp
    }
    names(L) <- d_names
    return(L)
}

