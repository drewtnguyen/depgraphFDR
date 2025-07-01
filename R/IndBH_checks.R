IndBH_pos_check <- function(check_type = 1, cached_results,
                            skip_comps = numeric(0)){
    setup = cached_results$setup
    ivs_lst = cached_results$ivs_lst
    Ns = cached_results$Ns
    if(check_type == 1){
        # Do a positive check based on block IndBH. It might be more powerful
        # sometimes than the coarse check.
        rej = .block_IndBH(cached_results)
    }
    if(check_type == 2){
        rej = IndBH_pos_check_coarse(setup, Ns$Ntot, Ns$Nplus)
    }
    if(check_type == 3){
        rej = IndBH_pos_check_fine(setup, ivs_lst, Ns,
                                    skip_comps = skip_comps)
    }

    checked_comps = get_chcked_comps(setup$comps_noname, rej)
    return(list(
        rej = rej,
        checked_comps = checked_comps))
}


# Return those indices that are not rejected by IndBH.
IndBH_neg_check <- function(setup, Ntot){
    BH_idx = setup$BH_idx
    alpha = setup$alpha
    m = setup$m
    use_pvals = setup$use_pvals

    thrs = 1:length(setup$BH_idx)

    if(all(thrs > Ntot)){
        R_UB = 0
    }
    else{
        R_UB = max(which(thrs <= Ntot), -Inf)
    }
    return(which(use_pvals > alpha * R_UB / m))
}

# Return those indices that are rejected by IndBH (coarse).
IndBH_pos_check_coarse <- function(setup, Ntot, Nplus){
    BH_idx = setup$BH_idx
    alpha = setup$alpha
    m = setup$m
    use_pvals = setup$use_pvals

    thrs = 1:length(BH_idx)

    R_LB = max(which(thrs <= Ntot - Nplus + 1), -Inf)
    return(which(use_pvals <= alpha * R_LB / m))
}

# Return those indices that are rejected by IndBH (fine).
IndBH_pos_check_fine <- function(setup, ivs_lst, Ns,
                            skip_comps = integer(0)){
    m = setup$m # number of hypotheses
    alpha = setup$alpha # FDR level
    comps = setup$comps # list of components of BH graph
    BH_idx = setup$BH_idx
    BH_comps = setup$comps_noname
    nc = setup$nc # number of components of BH graph
    use_pvals = setup$use_pvals # p-values corresponding to BH rejections

    Ntot = Ns$Ntot

    thrs = 1:length(BH_idx)

    comps2check = setdiff(1:nc, skip_comps)

    R_LB = integer(length(use_pvals))
    for(k in comps2check){
        # get component indices (in BH subset)
        comp = BH_comps[[k]]
        # Retrieve the N vector for this component.
        Nk = Ns$Nlist[[k]]
        # Subtract it from Ntot
        N_no_k = Ntot - Nk + 1
        # Set R_LB[comp] to the computed max.
        # The skipped ones will be zero, which is fine.
        R_LB[comp] = max(which(thrs <= N_no_k), -Inf)
    }
    return(which(use_pvals <= alpha * R_LB / m))
}

IndBH_full_calculation <- function(rejectable, cached_results){
    alpha = cached_results$setup$alpha
    m = cached_results$setup$m
    membership = cached_results$setup$membership
    use_pvals = cached_results$setup$use_pvals
    thrs = 1:length(cached_results$setup$BH_idx)
    BH_comps = cached_results$setup$comps_noname
    comps_connected = cached_results$setup$comps_connected
    ivs_lst = cached_results$ivs_lst
    Ns = cached_results$Ns
    Ntot = Ns$Ntot

    R_true = integer(length(use_pvals))
    for(i in rejectable){
        k = membership[i]
        comp = BH_comps[[k]]
        comp_ivs = ivs_lst[[k]]
        comp_con = comps_connected[k]

        N_i = compute_Ncomp(alpha, use_pvals = use_pvals,
                           thrs = thrs,
                           m = m,
                           comp = comp,
                           comp_ivs = comp_ivs,
                           comp_connected = comp_con,
                           inclusion_restriction = i)
        Nk = Ns$Nlist[[k]]
        # The skipped ones will be zero, which is fine.
        R_true[i] = max(which(thrs <= Ntot - Nk + N_i), -Inf)
    }
    return(which(use_pvals <= alpha * R_true / m))
}
