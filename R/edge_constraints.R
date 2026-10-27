# Edge constraints for graph exploration.
#
# forbidden: d x d logical, forbidden[i, j] = TRUE means the edge i -> j can never be in the graph.
# required:  d x d logical, required[i, j] = TRUE means the edge i -> j is always in the graph.
# Forcing the direction between a and b is `forbidden_edges = cbind(b, a)` (the edge may be absent
# but can never point b -> a); add `required_edges = cbind(a, b)` to also force the edge to exist.

# Turn a user supplied constraint into a d x d logical matrix.
# `x` can be NULL, a d x d 0/1 (or logical) matrix, or a two column matrix / data.frame / list of
# (from, to) pairs given as trait indices or trait names.
edge_constraint_matrix <- function(x, d, trait_names = NULL, arg = "forbidden_edges") {
    out <- matrix(FALSE, nrow = d, ncol = d)
    if (is.null(x) || length(x) == 0) {
        return(out)
    }
    if (is.list(x) && !is.data.frame(x)) {
        if (!all(lengths(x) == 2)) stop(arg, ": each pair must have length 2 (from, to)")
        x <- do.call(rbind, x)
    }
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x) && nrow(x) == d && ncol(x) == d && !is.character(x) && all(x %in% c(0, 1, TRUE, FALSE))) {
        # d x d adjacency style matrix. Ambiguous only when d == 2 and a single pair is given;
        # a 1 x 2 matrix is never d x d for d == 2, so this is safe.
        out <- x != 0
    } else {
        if (!is.matrix(x) || ncol(x) != 2) {
            stop(arg, " must be a d x d 0/1 matrix or a two column matrix of (from, to) pairs")
        }
        resolve <- function(v) {
            if (is.character(v)) {
                if (is.null(trait_names)) stop(arg, ": trait names given but the data has no column names")
                ix <- match(v, trait_names)
                if (anyNA(ix)) stop(arg, ": unknown trait(s): ", paste(v[is.na(ix)], collapse = ", "))
                ix
            } else {
                ix <- as.integer(v)
                if (anyNA(ix) || any(ix < 1 | ix > d)) stop(arg, ": trait index out of range 1..", d)
                ix
            }
        }
        out[cbind(resolve(x[, 1]), resolve(x[, 2]))] <- TRUE
    }
    if (any(diag(out))) stop(arg, " cannot contain self edges")
    out
}

# Check a flattened graph string or 0/1 matrix against the constraints
violates_edge_constraints <- function(B, forbidden, required) {
    B <- B != 0
    any(B[forbidden]) || !all(B[required])
}

# Make an adjacency matrix consistent with the constraints.
# Starts from the required edges and greedily adds the remaining edges (largest `weight` first)
# as long as the graph stays acyclic. Returns a 0/1 matrix.
enforce_edge_constraints <- function(adj, forbidden, required, weight = abs(adj)) {
    d <- nrow(adj)
    B <- matrix(0, d, d)
    B[required] <- 1
    if (!is_dag(B)) {
        stop("The required edges contain a cycle.")
    }
    cand <- (adj != 0) & !forbidden & !required
    diag(cand) <- FALSE
    ix <- which(cand)
    ix <- ix[order(weight[ix], decreasing = TRUE)]
    for (k in ix) {
        B[k] <- 1
        if (!is_dag(B)) B[k] <- 0
    }
    B
}

check_edge_constraints <- function(forbidden, required) {
    if (any(forbidden & required)) {
        stop("Some edges are both forbidden and required.")
    }
    if (any(required) && !is_dag(required * 1)) {
        stop("The required edges contain a cycle.")
    }
    invisible(TRUE)
}
