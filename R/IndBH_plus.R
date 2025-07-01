IndBH_plus <- function(recurse = 1,
                        alpha, pvals, adjlist, block = NULL){
    # Build a cache
    cached_results = build_cache(alpha, pvals, adjlist, block)
    if(cached_results$setup$nc == 0){
        message("BH made no rejections.")
        return(integer(0))  
    } 
    # Run IndBH_plus based on cache
    return(.IndBH_plus(recurse = recurse, cached_results))
}


.IndBH_plus <- function(recurse = 1,
                       cached_results = list(setup = NULL,
                                             ivs_lst = NULL,
                                             Ns = NULL),
                       return_BH_idx = FALSE){

    # Deal with recursions:
    # If recurse > 1, then define the method to be calibrated as the result of
    # one_step_calibration, one level of recursion down
    if(recurse > 1){
        method_to_calibrate <- function(cached_results, return_BH_idx){
            .IndBH_plus(recurse = recurse - 1, cached_results = cached_results,
                       return_BH_idx = return_BH_idx)
        }
    }

    # Otherwise, the method to be calibrated is IndBH
    else{
        method_to_calibrate <- .IndBH
    }



    # Get the initial rejection set of the method to be calibrated (in BH subset)
    rej_prev_init_BH = method_to_calibrate(cached_results = cached_results,
                                   return_BH_idx = TRUE)




    # Main function
    block = cached_results$setup$block
    if(block){
        # Use the block dependence version of one-step calibration
        return(.block_OSC(rej_prev_init_BH,
                            method_to_calibrate,
                            cached_results,
                            return_BH_idx))
    }
    else{
        # Get a block-dependence based rejection
        block_cache = cached_results
        block_cache$setup$block = TRUE
        rej_block_init_BH = .block_OSC(rej_prev_init_BH,
                                  method_to_calibrate,
                                  block_cache,
                                  return_BH_idx = TRUE)

        return(.nonblock_OSC(rej_prev_init_BH,
                             rej_block_init_BH,
                          method_to_calibrate,
                          cached_results,
                          return_BH_idx))
    }
}


.block_OSC <- function(rej_init_BH,
                           method_to_calibrate,
                           cached_results,
                           return_BH_idx){
    # Define all required variables
    setup = cached_results$setup
    m = setup$m # number of hypotheses
    alpha = setup$alpha # FDR level
    pvals = setup$pvals # p-values
    BH_idx = setup$BH_idx # BH rejections
    BH_pvals = setup$BH_pvals # p-values corresponding to BH rejections
    BH_ig = setup$BH_ig # BH graph
    clu = setup$clu # igraph connected components object of BH graph
    membership = setup$membership # membership vector of BH graph
    comps = setup$comps # list of components of BH graph
    comps_noname = setup$comps_noname
    nc = setup$nc # number of components of BH graph
    use_pvals = setup$use_pvals

    # Get the indices of the minimum p-values (in the BH subset)
    idx_mins = get_idx_of_mins(use_pvals, comps_noname)
    #

    # get which components are confirmed rejected (Positive Check)
    comp_allrej = get_chcked_comps(comps_noname, rej_init_BH)

    # get which components are confirmed NOT rejected (Negative Check)
    nrej_BH = which(use_pvals > alpha * length(rej_init_BH) / m)
    comp_allnrej = get_chcked_comps(comps_noname, nrej_BH)

    # For each component, update the threshold,
    # taking advantage of the block dependence shortcut,
    # and add the new rejections to the threshold.

    # Add in all the rejections we have already.
    rej_BH = rej_init_BH

    for(k in 1:nc){
        # get the indices of the vertices in the kth component (in the BH
        # subset)
        cp_BH = as.integer(comps_noname[[k]])
        # check if we can skip this component
        if(comp_allrej[k] | comp_allnrej[k]){
            next
        }
        # Make a new cache object, setting ALL the pvalues in cp_BH to 1.
        # TODO this is a shortcut to the full rejection,
        # which we can do because...?
        masked_cached_results = modify_cache_pvals(cached_results = cached_results,
                                                 new_mask = cp_BH,
                                                 new_modified_comp = k)
        # Get the threshold to calibrate for this block (in the BH subset)
        rhat_k = method_to_calibrate(cached_results = masked_cached_results,
                                     return_BH_idx = TRUE)
        # calib_pvals = use_pvals[c(idx_mins[k], rhat_k)]
        # # This threshold is equivalent to rejecting p_i based on
        # # use_pvals[c(i, rhat_k)].
        # thr_k = alpha * rbh(alpha, calib_pvals, m = m) / m
        thr_k = alpha * (length(rhat_k) + 1) / m
        cp_BH_rejs = cp_BH[which(use_pvals[cp_BH] <= thr_k)]
        rej_BH = c(rej_BH, cp_BH_rejs)
    }
    # Get the final rejection set
    rej_BH = unique(rej_BH)
    return(rejections_to_return(rej_BH, return_BH_idx, BH_idx))
}

.nonblock_OSC <- function(rej_prev_init_BH,
                          rej_block_init_BH,
                         method_to_calibrate,
                         cached_results,
                         return_BH_idx){

    setup = cached_results$setup
    BH_idx = setup$BH_idx

    # Do checks.
    rejectable = 1:length(BH_idx)

    # First a negative check.
    neg_notrejected = IndBH_plus_neg_check(setup, rej_prev_init_BH)
    all_checked = neg_notrejected
    all_rej = integer(0)
    rejectable = setdiff(rejectable, neg_notrejected)


    # Next, do positive checks.
    check_type = 1
    still_checking = TRUE
    rej_init_BH = union(rej_prev_init_BH, rej_block_init_BH)

    while(still_checking){

        pos_res = IndBH_plus_pos_check(check_type = check_type,
                                       rej_init_BH,
                                       method_to_calibrate,
                                       cached_results,
                                       rejectable)
        pos_rejected = pos_res$rej
        old_rejectable = rejectable
        rejectable = pos_res$rejectable

        all_checked = union(all_checked, pos_rejected)
        all_rej = union(all_rej, pos_rejected)


        if(length(all_checked) == length(BH_idx)){
            # Alternately, length(rejectable) == 0.
            # We're done checking.
            return(rejections_to_return(all_rej, return_BH_idx, BH_idx))

        }
        if(check_type == 1){
            check_type = 2
        }
        # Give up if the positive checks stop helping
        if(length(old_rejectable) == length(rejectable)){
            still_checking = FALSE
        }
    }

    # Execute the full calculation for the remaining rejectable p-values
    membership = cached_results$setup$membership
    BH_ig = cached_results$setup$BH_ig
    use_pvals = cached_results$setup$use_pvals
    rej_full_calc = numeric(0)
    alpha = setup$alpha
    m = setup$m
    for(i in rejectable){
        k = membership[i]
        mask = unique(as.integer(igraph::neighbors(BH_ig, i)))
        # Drop values where the mask would cover i
        mask = mask[mask != i]
        masked_cached_results =
            modify_cache_pvals(cached_results = cached_results,
                           new_mask = mask,
                           new_modified_comp = k)
        rhat_i = method_to_calibrate(cached_results = masked_cached_results,
                                     return_BH_idx = TRUE)
        thr_i = alpha * (length(union(rhat_i, i))) / m
        if(use_pvals[i] <= thr_i){
            rej_full_calc = c(rej_full_calc, i)
        }
    }
    all_rej = union(all_rej, rej_full_calc)
    return(rejections_to_return(all_rej, return_BH_idx, BH_idx))
}
