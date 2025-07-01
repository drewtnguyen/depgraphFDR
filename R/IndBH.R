#' Run the Independent Set BH (IndBH) procedure
#'
#' @param alpha FDR level
#' @param pvals A vector of p-values for the null hypotheses
#' @param adjlist An adjacency list, encoding dependencies between p-values as
#' a graph
#' @param block If TRUE, use the block dependence version of the procedure
#' for faster computation
#' @param cached_results Cached results from previous runs of the procedure
#'
#' @return If return_cache = FALSE, a vector of integers, representing
#' the indices of which nulls to reject. If return_cache = TRUE, a list of
#' two elements: the first is the vector of indices of which nulls to reject,
#' and the second is a IndBH_cache object, which can be passed to future
#' calls of IndBH to speed up computation.
#'
#' @export
#'
#' @examples
IndBH <- function(alpha, pvals, adjlist, block = NULL){
    # Build a cache
    cached_results = build_cache(alpha, pvals, adjlist, block)
    if(cached_results$setup$nc == 0){
        message("BH made no rejections.")
        return(integer(0))  
    } 
    # Run IndBH based on cache
    return(.IndBH(cached_results))
}


.IndBH <- function(cached_results = list(setup = NULL,
                                        ivs_lst = NULL,
                                        Ns = NULL),
                           return_BH_idx = FALSE){
    block = cached_results$setup$block
    if(block){
        # use the block dependence version of IndBH, which drastically
        # simplifies the problem
        return(.block_IndBH(cached_results, return_BH_idx))
    }
    else{
        return(.nonblock_IndBH(cached_results, return_BH_idx))
    }
}



.block_IndBH <- function(cached_results, return_BH_idx = TRUE){
    setup = cached_results$setup

    m = setup$m # number of hypotheses
    alpha = setup$alpha # FDR level
    pvals = setup$pvals # p-values
    BH_idx = setup$BH_idx # BH rejections
    BH_comps = setup$comps_noname
    nc = setup$nc # number of components of BH graph
    modified_comps = setup$modified_comps
    pmins_orig = setup$pmins
    use_pvals = setup$use_pvals # p-values corresponding to BH rejections

    pmins = numeric(nc)
    checked_comp = logical(nc)

    for(k in setup$modified_comps){
        checked_comp[k] = TRUE
        # get the indices of the BH vertices in the kth component
        idx = BH_comps[[k]]
        # extract the minimum p-value and add to list
        p = min(use_pvals[idx])
        pmins[k] = p
    }

    for(k in 1:nc){
        if(!checked_comp[k]){
            pmins[k] = pmins_orig[k]
        }
    }

    # get the IndBH threshold
    Rmins = rbh(alpha = alpha, pvals = pmins, m = m)

    thr = alpha * Rmins / m
    rej_BH = which(use_pvals <= thr)
    if(return_BH_idx){
        return(rej_BH)
    }
    else{
        rej = sub_to_orig(rej_BH, BH_idx)
        return(rej)
    }
}


.nonblock_IndBH <- function(cached_results, return_BH_idx = TRUE){

    setup = cached_results$setup
    BH_idx = setup$BH_idx
    BH_comps = setup$comps_noname

    if(length(setup$modified_comps) > 0){
        cached_results$Ns = update_Ns(cached_results)
    }

    Ntot = cached_results$Ns$Ntot
    # Do the checks.
    # TODO These checks don't matter much, it seems, for these dense graphs
    # where the number of ivs sets is quite small. THe speed savings of the
    # checks are just 330 to 300 ms. 10%? It's dwarfed by the computation of
    # the Ns.

    # Set up variables
    check_type = 1
    skip_comps = numeric(0)

    # First a negative check.
    neg_notrejected = IndBH_neg_check(setup, Ntot)
    all_checked = neg_notrejected
    all_rej = integer(0)

    skip_comps = c(skip_comps,
                   get_chcked_comps(BH_comps, neg_notrejected))

    # Do the positive checks
    while(check_type <= 3){
        pos_res = IndBH_pos_check(check_type = check_type,
                                       cached_results,
                                       skip_comps)
        pos_rejected = pos_res$rej
        all_checked = union(all_checked, pos_rejected)
        all_rej = union(all_rej, pos_rejected)
        if(length(all_checked) == length(BH_idx)){
            # We're done checking
            return(rejections_to_return(all_rej, return_BH_idx, BH_idx))
        }
        skip_comps = c(skip_comps, pos_res$checked_comps)
        check_type = check_type + 1
    }

    # If we haven't returned by now, we're not done checking,
    # and can no longer check component by component, but
    # must now check each p-value individually.

    # Execute the full calculation for the remaining p-values
    rejectable = setdiff(1:length(BH_idx), all_checked)
    rej_full = IndBH_full_calculation(rejectable, cached_results)
    all_rej = union(all_rej, rej_full)
    return(all_rej)
}
