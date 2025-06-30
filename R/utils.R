rbh <- function(alpha, pvals, m = NULL){

    if(is.null(m)){
        m = length(pvals)
    }
    pvals_sort = sort(pvals)
    v = which(pvals_sort <= alpha*(1:length(pvals))/m)
    if(length(v) == 0){
        res = 0
    }
    else{
        res = max(v)
    }
    return(res)
}


BH <- function(alpha, pvals, m = NULL){
    if(is.null(m)){
        m = length(pvals)
    }
    which(pvals <= alpha * rbh(alpha, pvals, m)/m)
}



# sub_idx is a vector taking values in 1...m,
# representing a subset of 1...m.
#
# for an index / vector i in 1:length(sub_idx),
# representing a coordinate of sub_idx---as if sub_idx
# were relabeled from 1...length(sub_idx)---
# this function maps i back to the original index in 1...m.
sub_to_orig <- function(i, sub_idx){
    sub_idx[i]
}

orig_to_sub <- function(i, sub_idx){
    match(i, sub_idx)
}


# this function gets the index of the minimum p-value
# in each connected component.
get_idx_of_mins <- function(pvals, comps){
    nc = length(comps)
    idx_mins = integer(nc)
    for(k in 1:nc){
        # get the indices of the vertices in the kth component
        idx = as.integer(comps[[k]])
        # get the location of the minimum in this component
        loc_of_min = which.min(pvals[idx])
        # get the index of the minimum in the original p-value vector
        idx_of_min = idx[loc_of_min]
        # extract the minimum p-value and add to list
        p = pvals[idx_of_min]
        idx_mins[k] = idx_of_min
    }
    return(idx_mins)
}

`%fast_in%` <- function(a,b){
    fastmatch::`%fin%`(a,b)
}

# check if the component has been entirely checked
get_chcked_comps <- function(comps, chcked_idx){
    nc = length(comps)
    chcked_comps = logical(nc)
    for(k in 1:nc){
        cp = as.integer(comps[[k]])
        if(all(cp %fast_in% chcked_idx)){
            # if so, report it as checked
            chcked_comps[k] = TRUE
        }
    }
    return(chcked_comps)
}


rejections_to_return <- function(rej_BH, return_BH_idx, BH_idx){
    if(return_BH_idx){
        return(rej_BH)
    }
    else{
        rej = sub_to_orig(rej_BH, BH_idx)
        return(rej)
    }


}

