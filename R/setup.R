build_cache <- function(alpha, pvals, adjlist, block = NULL){
    cached_results = list(setup = NULL,
                          ivs_lst = NULL,
                          Ns = NULL)
    setup = graphical_procedure_setup(alpha, pvals, adjlist, block)
    cached_results$setup = setup
    block = setup$block

    # If not block, augment cached results with independent set calculations
    if(!block){
        cached_results$ivs_lst = compute_ivs_all(cached_results)
        cached_results$Ns = compute_Ns(cached_results)
    }
    return(cached_results)
}

# run boilerplate code required for every graphical procedure
graphical_procedure_setup <- function(alpha, pvals, adjlist, block = NULL){
    # define constants
    m = length(pvals)
    # Convert adjlist to igraph object.
        # I want this code to work:
        # ig = igraph::graph_from_adj_list(adjlist, mode = 'total')
        # But see issue
        # https://github.com/igraph/igraph/issues/2382
        # for why code is written this way instead.
    ig = igraph::as.undirected(
        igraph::graph_from_adj_list(adjlist, mode='out'),
        mode='mutual')

    # make a copy to stay unnamed
    ig_copy = ig
    # rename the vertices of the graph to remember the mapping back
    # to the p-value vector. very important for indexing later.
    igraph::V(ig)$name = as.character(1:m)

    # reduce problem to the BH rejection set, and extract mapping back
    # to the original problem
    BH_idx = BH(alpha, pvals)

    BH_pvals = pvals[BH_idx]
    # subset the ig as well
    BH_ig = igraph::induced_subgraph(ig, BH_idx)
    BH_ig_noname = igraph::induced_subgraph(ig_copy, BH_idx)
    # get the connected components of the graph
    clu = igraph::components(BH_ig)
    membership = clu$membership
    comps = igraph::groups(clu)

    # get components of the copied graph
    clu_noname = igraph::components(BH_ig_noname)
    comps_noname = igraph::groups(clu_noname)

    # get the number of connected components
    nc = clu$no

    # Check if they're fully connected
    is_complete <- function(g) {
        g = igraph::simplify(g)
        n <- igraph::vcount(g)
        m <- igraph::gsize(g)          # actual number of edges
        m == n*(n-1)/2
    }

    comps_connected = logical(nc)
    for(k in seq_len(nc)){
        comp_ig = igraph::induced_subgraph(BH_ig, comps_noname[[k]])
        comps_connected[k] = is_complete(comp_ig)
    }

    # Get the mininum p-values
    pmins = logical(nc)
    for(k in seq_len(nc)){
        # get the indices of the BH vertices in the kth component
        idx = comps_noname[[k]]
        # extract the minimum p-value and add to list
        p = min(BH_pvals[idx])
        pmins[k] = p
    }
    if(is.null(block)){
        block = FALSE
        if(all(comps_connected)){
            block = TRUE
            message("Problem instance is compatible with the cluster graph
                computational shortcut.")
        }
    }

    return(list(
        m = m,
        alpha = alpha,
        pvals = pvals,
        BH_idx = BH_idx,
        BH_pvals = BH_pvals,
        BH_ig = BH_ig,
        clu = clu,
        clu_noname = clu_noname,
        membership = membership,
        comps = comps,
        comps_noname = comps_noname,
        comps_connected = comps_connected,
        pmins = pmins,
        nc = nc,
        use_pvals = BH_pvals,
        masked_BH = integer(0),
        modified_comps = integer(),
        block = block
    ))
}

# just updates use_pvals to record the masking
modify_cache_pvals <- function(cached_results, new_mask, new_modified_comp){
    setup = cached_results$setup
    # record the masking
    masked_BH = union(setup$masked_BH, new_mask)
    modified_comps = union(setup$modified_comps, new_modified_comp)

    # mask the p-values
    masked_BH_pvals = setup$use_pvals
    masked_BH_pvals[masked_BH] = 1

    # define new setup
    new_setup = setup
    new_setup$use_pvals = masked_BH_pvals
    new_setup$masked_BH = masked_BH
    new_setup$modified_comps = modified_comps

    # modify the cache and return it
    new_cached_results = cached_results
    new_cached_results$setup = new_setup
    return(new_cached_results)
}

# updates the Ns object---but not Nplus.
update_Ns <- function(cached_results){
    old_Ns = cached_results$Ns
    setup = cached_results$setup
    ivs_lst = cached_results$ivs_lst

    m = setup$m # number of hypotheses
    alpha = setup$alpha # FDR level
    use_pvals = setup$use_pvals
    modified_comps = setup$modified_comps

    BH_idx = setup$BH_idx # BH rejections
    BH_ig = setup$BH_ig # BH graph
    BH_comps = setup$comps_noname
    comps_connected = setup$comps_connected
    nc = setup$nc # number of components of BH graph

    thrs = 1:length(BH_idx)
    Ntot = old_Ns$Ntot
    Nlist = old_Ns$Nlist
    for(k in modified_comps){
        # get component indices (in BH subset)
        comp = BH_comps[[k]]
        # get component independent set indices (in BH subset)
        comp_ivs = ivs_lst[[k]]
        comp_con = comps_connected[k]
        # Compute the N vector for this component
        Nk = compute_Ncomp(alpha, use_pvals = use_pvals,
                           thrs = thrs,
                           m = m,
                           comp = comp,
                           comp_ivs = comp_ivs,
                           comp_connected = comp_con)
        # Re-define Nlist[[k]] to be the new vector
        Nlist[[k]] = Nk
        Ntot = Ntot + Nk - old_Ns$Nlist[[k]]
    }
    return(list(Ntot = Ntot, Nplus = old_Ns$Nplus, Nlist = Nlist))
}

compute_ivs_all <- function(cached_results){
    comps = cached_results$setup$comps
    ig = cached_results$setup$BH_ig
    # for each component, get collection of all maximal independent sets
    ivs_lst = vector(mode = 'list', length = length(comps))
    for(i in 1:length(comps)){
        comp = comps[[i]]
        gcomp = igraph::induced_subgraph(ig, comp)
        comp_ivs = igraph::maximal_ivs(gcomp)
        ivs_lst[[i]] = comp_ivs
    }
    return(ivs_lst)
}


compute_Ncomp <- function(alpha, use_pvals,
                          thrs = 1:length(use_pvals), m, comp, comp_ivs,
                          comp_connected,
                          inclusion_restriction = integer(0)){
    # initialize contribution Ncomp
    Ncomp = numeric(length = length(thrs))
    # if comp, the connected component, is a singleton,
    # the contribution N of this component at thresh r
    # is 1 if p <= alpha*r/m and 0 otherwise
    if(length(comp) == 1){
        p = use_pvals[comp]
        Ncomp = as.numeric(p <= alpha*thrs/m)
    }
    else if(comp_connected){
        p = min(use_pvals[comp])
        Ncomp = as.numeric(p <= alpha*thrs/m)
    }
    else{
        for(iset in comp_ivs){
            if(length(inclusion_restriction) > 0){
                if(!(inclusion_restriction %in% iset)){
                    next
                }
            }
            # get pvals corresponding to iset
            p = use_pvals[comp[iset]]
            # calculate smallest r s.t. p <= alpha*r/m.
            r = ceiling(m*p/alpha)
            # Handling zero p-values: replace any 0's in r with 1's.
            r[r == 0] = 1
            # Handling 1 p-values: drop r that is just too large.
            r = r[r <= m]
            # make a table rawtab of how many occurrences of each r
            # there are, for each of the bins between each of the thresholds
            # also see stackoverflow link, though data.table gave me an error.
            #
            rawtab = tabulate(r)
            # Allocating memory can take a long time if length(thrs) is large
            tab = integer(length(thrs))
            if(length(r) > 0){
                # This part assumes that thrs = 1:length(use_pvals)
                tab[1:max(r)] = rawtab
            }
            # calculate how many pvalues are less than alpha*thr/m
            Nthis = cumsum(tab)
            # if this iset, for some r, has a larger number of p-values less than
            # alpha*r/m, overwrite Ncomp
            Ncomp = pmax(Nthis, Ncomp)
        }
    }
    return(Ncomp)
}


# Compute the Nplus and Ntot, before any p-values are masked
compute_Ns <- function(cached_results){
    setup = cached_results$setup
    ivs_lst = cached_results$ivs_lst

    m = setup$m # number of hypotheses
    alpha = setup$alpha # FDR level
    BH_pvals = setup$BH_pvals # p-values corresponding to BH rejections

    BH_idx = setup$BH_idx # BH rejections
    BH_ig = setup$BH_ig # BH graph
    BH_comps = setup$comps_noname
    comps_connected = setup$comps_connected
    nc = setup$nc # number of components of BH graph

    thrs = 1:length(BH_idx)
    Ntot = numeric(length(thrs))
    Nplus = numeric(length(thrs))
    Nlist = vector(mode = 'list', length = nc)

    for(k in 1:nc){
        # get component indices (in BH subset)
        comp = BH_comps[[k]]
        # get component independent set indices (in BH subset)
        comp_ivs = ivs_lst[[k]]
        comp_con = comps_connected[k]
        # Compute the N vector for this component
        Nk = compute_Ncomp(alpha, use_pvals = BH_pvals,
                           thrs = thrs,
                           m = m,
                           comp = comp,
                           comp_ivs = comp_ivs,
                           comp_connected = comp_con)
        Nlist[[k]] = Nk
        Ntot = Ntot + Nk
        Nplus = pmax(Nplus, Nk)
    }
    return(list(Ntot = Ntot, Nplus = Nplus, Nlist = Nlist))
}
