# Return those indices that are not rejected by IndBH_plus.
IndBH_plus_neg_check <- function(setup, rej_init_BH){
    BH_idx = setup$BH_idx
    alpha = setup$alpha
    m = setup$m
    use_pvals = setup$use_pvals
    R_UB = length(rej_init_BH) + 1
    return(which(use_pvals > alpha * R_UB / m))
}

IndBH_plus_pos_check <- function(check_type,
                                 rej_init_BH,
                                 method_to_calibrate,
                                 cached_results,
                                 rejectable){
    if(check_type == 1){
        rej = rej_init_BH
    }

    if(check_type == 2){
        rej = IndBH_plus_pos_check_from_mask(cached_results,
                                              method_to_calibrate,
                                              rejectable)
    }
    rejectable = setdiff(rejectable, rej)
    return(list(rej = rej,
                rejectable = rejectable))
}


# for all rejectable pvalues
IndBH_plus_pos_check_from_mask <- function(cached_results, method_to_calibrate,
                                           rejectable){
    setup = cached_results$setup

    alpha = setup$alpha
    m = setup$m
    BH_ig = setup$BH_ig
    membership = setup$membership
    use_pvals = setup$use_pvals

    # mask only the rejectable p-values
    mask = do.call(c, igraph::adjacent_vertices(BH_ig, rejectable))
    mask = unique(as.integer(mask))
    ks = unique(membership[rejectable])
    masked_cached_results =
        modify_cache_pvals(cached_results = cached_results,
                           new_mask = mask,
                           new_modified_comp = ks)
    R_LB = length(method_to_calibrate(cached_results = masked_cached_results,
                                      return_BH_idx = TRUE)) + 1
    return(which(use_pvals <= alpha * R_LB / m))
}


