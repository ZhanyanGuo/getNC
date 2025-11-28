#' Sweep partner knockouts for a single target
#'
#' For each candidate partner gene \code{g} (excluding the target and already-knocked genes),
#' compute \eqn{\mathbb{E}[ \text{target} \mid (\text{knocked} \cup \{g\}) = 0 ]} and
#' \eqn{\mathrm{Var}[ \text{target} \mid (\text{knocked} \cup \{g\}) = 0 ]}.
#' Inputs mirror \code{\link{predict_conditional_knockout}}.
#'
#' @param Sigma Covariance matrix (features x features).
#' @param genes Character vector of feature names aligned with `Sigma`.
#' @param target Character or integer (must specify exactly one target gene).
#' @param knocked Character or integer vector of genes already set to 0.
#' @param mu Optional mean vector aligned with `genes`; if NULL, zeros used.
#'
#' @return list with two named numeric vectors (same names/order):
#'   \itemize{
#'     \item \code{mean}: conditional means of the target under each partner knockout
#'     \item \code{var}:  conditional variances of the target under each partner knockout
#'   }
#'
#' @examples
#' # Minimal synthetic example (3 genes): target = G1; knock out G2 or G3
#' set.seed(1)
#' genes  <- c("G1","G2","G3")
#' R      <- matrix(c(1, 0.3, -0.2,
#'                    0.3, 1,   0.1,
#'                   -0.2, 0.1, 1), 3, 3, dimnames = list(genes, genes))
#' sdv    <- c(1.0, 0.8, 1.2)
#' Sigma  <- diag(sdv) %*% R %*% diag(sdv)  # covariance in original units
#' out    <- sweep_partner_knockouts(Sigma, genes, target = "G1", knocked = "G2", mu = c(0,0,0))
#' str(out)
#'
#' @seealso \code{\link{predict_conditional_knockout}}
#'
#' @references
#' Anderson, T.W. (2003). \emph{An Introduction to Multivariate Statistical Analysis} (3rd ed.).
#'   Wiley. (Conditional Normal formulae)
#' Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). \emph{Multivariate Analysis}. Academic Press.
#'
#' @export
sweep_partner_knockouts <- function(Sigma,
                                    genes,
                                    target,
                                    knocked = character(0),
                                    mu = NULL) {
  # basic checks
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("`Sigma` must be a square matrix.")
  }
  if (length(genes) != nrow(Sigma)) {
    stop("`genes` length must match nrow(Sigma).")
  }

  # resolve indices
  to_index <- function(sel) if (is.numeric(sel)) sel else match(sel, genes)

  t_idx <- to_index(target)
  if (length(t_idx) != 1L || is.na(t_idx)) {
    stop("`target` must refer to exactly one valid gene.")
  }

  b_idx <- integer(0)
  if (length(knocked)) {
    b_idx <- to_index(knocked)
    if (any(is.na(b_idx))) stop("Some `knocked` genes not found in `genes`.")
  }

  # candidate partners: all genes except target and already-knocked
  all_idx <- seq_along(genes)
  cand_idx <- setdiff(all_idx, c(t_idx, b_idx))
  if (length(cand_idx) == 0L) {
    warning("No candidate partner genes found (everything excluded).")
    return(list(mean = numeric(0), var = numeric(0)))
  }

  # pre-allocate
  partner_names <- genes[cand_idx]
  means <- numeric(length(cand_idx))
  vars  <- numeric(length(cand_idx))
  names(means) <- partner_names
  names(vars)  <- partner_names

  # loop over partners, add each to the knockout set, call conditioning
  base_knocked <- if (length(b_idx)) genes[b_idx] else character(0)
  for (k in seq_along(cand_idx)) {
    gname <- genes[cand_idx[k]]
    res <- predict_conditional_knockout(
      Sigma  = Sigma,
      genes  = genes,
      target = genes[t_idx],
      knocked = c(base_knocked, gname),
      mu     = mu
    )
    # extract scalar mean/var for the (single) target
    means[k] <- unname(res$mean_cond[1])
    vars[k]  <- unname(res$cov_cond[1, 1])
  }

  list(mean = means, var = vars)
}

#' 3D density plots from a glasso fit, ranked by mean and by variance (colored by magnitude)
#'
#' Produces two interactive 3D line plots (via plotly): top-K partners ranked by
#' \code{|conditional mean|} and top-K partners ranked by \code{conditional variance}.
#' Colors encode the magnitude used for ranking (|mean| or variance respectively).
#'
#' @param fit List from \code{auto_fit_glasso()}/\code{fit_glasso()} (needs \code{$sigma}, \code{$features};
#'   for original units also needs \code{$sd}, and optionally \code{$mu}).
#' @param target Single gene (name or index) to predict.
#' @param knocked Vector of already-knocked genes (names or indices).
#' @param k Integer, how many partners to plot for each ranking (default 15).
#' @param use_original_units TRUE to recover covariance via sd (\eqn{\Sigma = D R D}). If FALSE, use corr scale.
#' @param renormalize TRUE to \code{cov2cor(fit$sigma)} before scaling when recovering \eqn{\Sigma}.
#' @param xlim Optional length-2 numeric x-range for the pdf curves; NULL = auto.
#' @param n Number of x points per curve (default 200).
#'
#' @return list(mean_plot = plotly_object, var_plot = plotly_object)
#'
#' @examples
#' \donttest{
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   # Tiny synthetic "fit" object on 6 genes to keep examples fast:
#'   set.seed(2)
#'   genes <- paste0("G", 1:6)
#'   R     <- diag(6)
#'   R[1,2] <- R[2,1] <- 0.4
#'   R[1,3] <- R[3,1] <- -0.3
#'   sdv   <- runif(6, 0.8, 1.2)
#'   Sigma <- diag(sdv) %*% R %*% diag(sdv)
#'   fit   <- list(
#'     sigma    = stats::cov2cor(Sigma),  # pretend we estimated corr
#'     sd       = sdv,
#'     mu       = rep(0, 6),
#'     features = genes
#'   )
#'   plt <- plot_partner_knockout_densities_dual(
#'     fit, target = "G1", knocked = "G2", k = 3,
#'     use_original_units = TRUE
#'   )
#'   # In interactive R, print to view:
#'   # plt$mean_plot; plt$var_plot
#' }
#' }
#'
#' @seealso \code{\link{recover_covariance}}, \code{\link{predict_conditional_knockout}},
#'   \code{\link{auto_fit_glasso}}, \code{\link{fit_glasso}}
#'
#' @references
#' C. Sievert (2020). \emph{Interactive Web-Based Data Visualization with R, plotly, and shiny}.
#'   Chapman & Hall/CRC.
#' Anderson, T.W. (2003). \emph{An Introduction to Multivariate Statistical Analysis} (3rd ed.).
#'   Wiley.
#'
#' @importFrom stats dnorm cov2cor
#' @import plotly
#' @export
plot_partner_knockout_densities_dual <- function(
  fit, target, knocked = character(0),
  k = 15,
  use_original_units = TRUE,
  renormalize = TRUE,
  xlim = NULL,
  n = 200
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. install.packages('plotly')")
  }

  # Build Sigma, genes, mu from the fit
  genes <- fit$features
  if (isTRUE(use_original_units)) {
    if (is.null(fit$sd)) stop("fit$sd missing; set use_original_units=FALSE or refit storing sd.")
    Sigma <- recover_covariance(fit, renormalize = renormalize)
    mu    <- if (!is.null(fit$mu)) fit$mu else rep(0, length(genes))
  } else {
    Sigma <- stats::cov2cor(fit$sigma)
    mu    <- rep(0, length(genes))   # z-scored/corr scale => zero mean
  }

  # Get sweep results (means/vars per partner)
  sv <- sweep_partner_knockouts(Sigma, genes, target, knocked, mu)
  if (length(sv$mean) == 0L) stop("No candidate partners to plot.")

  partners_all <- names(sv$mean)
  m_all <- as.numeric(sv$mean)
  v_all <- as.numeric(sv$var)
  sd_all <- sqrt(pmax(v_all, 0))

  # The following part involves chat-GPT 5 for handleling colors and orgnizing

  # color mapper (0..1 -> hex)
  .mk_mapper <- function() {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      ramp <- grDevices::colorRamp(viridisLite::viridis(256))
    } else {
      ramp <- grDevices::colorRamp(c("#2166AC", "#F4A582", "#B2182B")) # blue→salmon→red
    }
    function(x01) {
      x01 <- pmin(pmax(x01, 0), 1)
      rgb <- ramp(x01)
      grDevices::rgb(rgb[,1], rgb[,2], rgb[,3], maxColorValue = 255)
    }
  }
  map_col <- .mk_mapper()

  # Plot builder with per-trace colors by a magnitude vector
  .build_plot <- function(partners, m, sd, mag, title_suffix) {
    # normalize magnitudes for color mapping
    mag01 <- if (length(unique(mag)) > 1) {
      (mag - min(mag, na.rm = TRUE)) / (max(mag, na.rm = TRUE) - min(mag, na.rm = TRUE))
    } else {
      rep(0.5, length(mag))
    }
    cols <- map_col(mag01)

    # x-range
    if (is.null(xlim)) {
      x_min <- min(m - 4 * sd); x_max <- max(m + 4 * sd)
      if (!is.finite(x_min) || !is.finite(x_max) || x_min == x_max) {
        x_min <- -3; x_max <- 3
      }
      xrng <- c(x_min, x_max)
    } else {
      xrng <- xlim
    }
    xs <- seq(xrng[1], xrng[2], length.out = n)

    p <- plotly::plot_ly()
    for (i in seq_along(partners)) {
      dens <- stats::dnorm(xs, mean = m[i], sd = ifelse(sd[i] > 0, sd[i], 1e-8))
      p <- plotly::add_trace(
        p,
        x = xs,
        y = rep(i, length(xs)),   # stack by index; label with partner names
        z = dens,
        type = "scatter3d",
        mode = "lines",
        name = sprintf("%s | mean=%.3f sd=%.3f", partners[i], m[i], sd[i]),
        line = list(width = 4, color = cols[i])
      )
    }
    targ_name <- if (is.numeric(target)) genes[target] else target
    plotly::layout(
      p,
      title = paste("density of", targ_name, title_suffix),
      scene = list(
        xaxis = list(title = paste0("Target ", targ_name)),
        yaxis = list(title = "Partner",
                     tickmode = "array",
                     tickvals = seq_along(partners),
                     ticktext = partners),
        zaxis = list(title = "Density f(x)")
      ),
      legend = list(title = list(text = "Partner (color = magnitude)"))
    )
  }

  # Top-K by |mean|
  ord_mean <- order(abs(m_all), decreasing = TRUE)
  idx_mean <- ord_mean[seq_len(min(k, length(ord_mean)))]
  mean_plot <- .build_plot(
    partners = partners_all[idx_mean],
    m  = m_all[idx_mean],
    sd = sd_all[idx_mean],
    mag = abs(m_all[idx_mean]),                         # color by |mean|
    title_suffix = "(top-K means)"
  )

  # Top-K by variance
  ord_var <- order(v_all, decreasing = TRUE)
  idx_var <- ord_var[seq_len(min(k, length(ord_var)))]
  var_plot <- .build_plot(
    partners = partners_all[idx_var],
    m  = m_all[idx_var],
    sd = sd_all[idx_var],
    mag = v_all[idx_var],                               # color by variance
    title_suffix = "(top-K variance)"
  )

  list(mean_plot = mean_plot, var_plot = var_plot)
}

#' Get indices of the top-k correlated genes (including if exist target gene(s))
#'
#' @param Sigma Covariance (or correlation) matrix.
#' @param genes Character vector of gene names (length = ncol(Sigma)).
#' @param target Gene index(es) or name(s).
#' @param k Number of top genes to return (including the target).
#'
#' @return Integer vector of indices for the top-k genes.
#' @keywords internal
get_top_k_gene_indices <- function(Sigma, genes, target, k) {

  # --- safety checks ---
  stopifnot(is.matrix(Sigma), is.numeric(Sigma))
  stopifnot(length(genes) == ncol(Sigma))
  k <- min(max(0, k), length(genes))

  # Convert target to indices

  if (is.character(target)) {
    # names → indices
    if (!all(target %in% genes)) {
      stop("Some target genes not found in `genes`.")
    }
    target_idx <- match(target, genes)

  } else if (is.numeric(target)) {
    # numeric index
    if (!all(target >= 1 & target <= length(genes))) {
      stop("Numeric `target` contains invalid indices.")
    }
    target_idx <- as.integer(target)

  } else {
    stop("`target` must be character or numeric.")
  }

  target_idx <- unique(target_idx)

  # Compute association strength

  # Submatrix: correlations of all genes to target set
  v <- Sigma[, target_idx, drop = FALSE]

  # Compute association strength: strongest absolute correlation to ANY target gene
  assoc <- apply(v, 1, function(x) max(abs(x)))

  # Rank by absolute correlation (descending)
  ranked <- order(assoc, decreasing = TRUE)

  # Take top k
  top_k <- ranked[seq_len(min(k, length(ranked)))]

  unique(top_k)
}

#' Create a GLnode object
#'
#' @param index Integer. The gene index for this node.
#' @param neighbours Integer vector of neighbor gene indices.
#' @param knocked Logical. Whether the gene is knocked or not.
#'
#' @return An object of class "GLnode".
#' @export 
new_GLnode <- function(index, neighbours = integer(0), knocked = FALSE) {
  stopifnot(is.numeric(index), length(index) == 1)
  stopifnot(is.logical(knocked), length(knocked) == 1)

  node <- list(
    index      = as.integer(index),
    neighbours = as.integer(neighbours),
    knocked    = knocked
  )
  class(node) <- "GLnode"
  node
}

#' Create a GLgraph object
#'
#' @param fit A list from run_glasso()/fit_glasso()/fit_glasso_raw().
#'   Must contain: $sigma, $omega, $features.
#' @param k Integer. Number of genes to include.
#' @param target Gene name(s) or index(es).
#'
#' @return An object of class "GLgraph".
#' @export
new_GLgraph <- function(fit, k, target) {

  # basic input checks
  stopifnot(is.list(fit))
  stopifnot(all(c("sigma", "omega", "features") %in% names(fit)))
  stopifnot(is.numeric(k), length(k) == 1)

  # extract components
  genes <- fit$features
  Sigma <- fit$sigma
  Omega <- fit$omega

  # get the top-k gene indices based on correlation to target
  top_k_idx <- get_top_k_gene_indices(
    Sigma  = Sigma,
    genes  = genes,
    target = target,
    k      = k
  )

  # extract the k-by-k precision matrix
  Omega_sub <- Omega[top_k_idx, top_k_idx, drop = FALSE]

  # prepare a list to store node objects
  nodes <- vector("list", length(top_k_idx))

  # build each node
  for (i in seq(from = 1, to = length(top_k_idx))) {
    # global gene index
    global_idx <- top_k_idx[i]

    # row of precision matrix for this node
    row_i <- Omega_sub[i, ]

    # indices of neighbors in the submatrix
    neigh_local <- which(row_i != 0)

    # remove self
    neigh_local <- setdiff(neigh_local, i)

    # convert local indices to global indices
    neigh_global <- top_k_idx[neigh_local]

    # create node
    nodes[[i]] <- new_GLnode(
      index      = global_idx,
      neighbours = neigh_global
    )
  }

  g <- list(
    fit          = fit,
    k            = k,
    nodes        = nodes,
    selected_idx = top_k_idx
  )

  class(g) <- "GLgraph"
  return(g)
}

#' Convert GLgraph object to a visNetwork graph
#'
#' @param G A GLgraph object.
#' @param color_normal Color for unknocked nodes.
#' @param color_knocked Color for knocked nodes.
#'
#' @return A visNetwork object.
#' @export
visnetwork_from_GLgraph <- function(G,
                                    color_normal = "lightblue",
                                    color_knocked = "salmon") {

  # Ensure the input is a GLgraph
  stopifnot(inherits(G, "GLgraph"))

  # Build node table: id, label, knocked status
  idx_vec <- vapply(G$nodes, `[[`, integer(1), "index")
  nodes_df <- data.frame(
    id = idx_vec,
    label = G$fit$features[idx_vec],
    knocked = vapply(G$nodes, `[[`, logical(1), "knocked"),
    stringsAsFactors = FALSE
  )

  # Assign node colors
  nodes_df$color <- ifelse(nodes_df$knocked, color_knocked, color_normal)

  ### 1. Build edges using ONLY the upper triangle of Omega

  top_k_idx <- G$selected_idx
  Omega_sub <- G$fit$omega[top_k_idx, top_k_idx, drop = FALSE]

  edges_df <- data.frame(
    from  = integer(0),
    to    = integer(0),
    color = character(0),
    stringsAsFactors = FALSE
  )

  p <- length(top_k_idx)

  for (i in seq_len(p - 1)) {
    for (j in seq((i + 1), p)) {
      if (Omega_sub[i, j] != 0) {
        edges_df <- rbind(edges_df, data.frame(
          from  = top_k_idx[i],
          to    = top_k_idx[j],
          color = color_normal,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # locate knocked nodes
  node_idx_vec   <- vapply(G$nodes, `[[`, integer(1), "index")
  node_knocked   <- vapply(G$nodes, `[[`, logical(1), "knocked")
  knocked_global <- node_idx_vec[node_knocked]

  if (length(knocked_global) > 0 && nrow(edges_df) > 0) {
    for (kn in knocked_global) {
      hit <- (edges_df$from == kn | edges_df$to == kn)
      edges_df$color[hit] <- color_knocked
    }
  }

  # Construct and return the visNetwork graph
  visNetwork::visNetwork(nodes_df, edges_df) %>%
    visNetwork::visLayout(randomSeed = 1) %>%
    visNetwork::visOptions(
      highlightNearest = FALSE,
      nodesIdSelection = TRUE
    )
}

#' Toggle knock state of a node and update visNetwork graph
#'
#' @param G A GLgraph object.
#' @param clicked Integer ID of clicked node.
#' @param color_normal Color for unknocked nodes.
#' @param color_knocked Color for knocked nodes.
#'
#' @return Updated GLgraph object.
#' @export
visnetwork_toggle_knock <- function(G, clicked,
                                    color_normal = "lightblue",
                                    color_knocked = "salmon") {

  # Basic input checks: GLgraph and numeric node ID
  if (is.null(clicked)){
    return(visnetwork_from_GLgraph(G, color_normal, color_knocked))
  }
  stopifnot(inherits(G, "GLgraph"))
  stopifnot(is.numeric(clicked), length(clicked) == 1)

  # Identify which node matches the clicked ID
  idx <- which(vapply(G$nodes, `[[`, integer(1), "index") == clicked)

  # If not found, warn and return unchanged graph
  if (length(idx) == 0) {
    warning("Clicked node not found in GLgraph.")
    return(visnetwork_from_GLgraph(G, color_normal, color_knocked))
  }

  # Flip the knocked status (TRUE <-> FALSE)
  G$nodes[[idx]]$knocked <- !G$nodes[[idx]]$knocked

  # return the updated G
  G
}


# Gen Ai used for documentation and input verify