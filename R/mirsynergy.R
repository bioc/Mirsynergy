# Mirsynergy: 
# Solve miRNA synergistic regulatory network (MSRN) problem by
# clustering miRNA and genes with overlapping neighbourhood expansion
# The algorithm is similar to ClusterONE except the seed selection is
# restricted to miRNA only

# Input:
# W: N x M weight matrix of N mRNA and M miRNA
# H: N x N gene-gene interaction matrix including TFBS and PPI
# den: threshold for density score of each cluster
# ovp: threshold for overlap score between two clusters being merged

# Output: sparse binary matrix with (N+M) by C dimension, where C is
# the number of clusters

mirsynergy <- function(W, H, alpha=2, 
	merge.tol=0.8, density1.tol=1e-2, density2.tol=5e-3, verbose=FALSE) {
	
	# record running time
	if(verbose) ptm <- proc.time()
	
	# check weight matrix first
	W <- check_w(W)
	
	# Stage I get miRNA clusters
	if(verbose)
		message(sprintf("\n*Stage I. Cluster %s miRNA:\n", ncol(W)))
	
	# create miRNA-miRNA synergistic matrix
	S <- get_mirmir_w(W)
		
	V <- mircluster(S, alpha, merge.tol, density1.tol, verbose)
	
	if(verbose) message(sprintf("\nTime elapsed for Stage I: %.3f", 
		(proc.time() - ptm)[3]))
		
	# Stage II get gene clusters using miRNA clusters as seed	
	H <- as.matrix(H[rownames(W),rownames(W)])
		
	if(verbose)
		message(sprintf("\n*Stage II. Grow MRM with mRNA:\n", ncol(W)))
		
	V2 <- gencluster(V, W, H, alpha, density2.tol, verbose)
		
	if(verbose) {
		
		message("\nFinal MRM:")
		
		print_modules2(V2)
		
		message(sprintf("\nTotal time elapsed: %.3f (sec)", 
			(proc.time() - ptm)[3]))
	}
	
	V2
}

# check weight matrix
check_w <- function(W) {
		
	W.old <- W
	
	if(any(W>1) || any(W<0)) {
		
		stop(paste("Only probabilistic weights (0<=w<=1) are allowed!",	
		"Transform weights using normalize_weight.",sep="\n"))
	}
		
	W <- W[apply(W, 1, max) > 0, apply(W, 2, max) > 0]
		
	if(!identical(W.old, W))
		message("Discard row/column of all zeros in W")
		
	as.matrix(W)
}


# create miRNA-miRNA synergistic matrix
get_mirmir_w <- function(W) {
	
	# total outdegree weights per miRNA
	mir_w <- apply(W, 2, sum)
	
	# get connection weight b/w miRNA
	S <- t(W) %*% W
	
	# get min(sum_w_k, sum_w_j)
	min_w <- lapply(mir_w, function(w_sum_k) 
		apply(cbind(rep(w_sum_k, length(mir_w)),mir_w),1,min))
	
	min_w <- do.call("rbind", min_w)
	
	S / min_w	
}


# cluster miRNA based on miRNA-miRNA weight matrix
mircluster <- function(W, alpha, merge.tol, density.tol, verbose) {
	
	# Step1: Grow clusters
	# keep track of miRNA k that's being included in the modules
	V <- list()
	
	# index miRNA ID to enable numerical search
	# (more efficient than string search)
	original.id <- rownames(W)
	
	names(original.id) <- 1:nrow(W)
				
	dimnames(W) <- list(1:nrow(W), 1:nrow(W))
	
	# set self-connection weights to zero
	diag(W) <- 0
	
	# precompute weights for downstream speed-up
	w.sum <- as.numeric(apply(W,1,sum))
	
	# save neighbour nodes for each node as a list for downstream speed-up
	n.list <- lapply(split(W, as.numeric(rownames(W))), 	
		function(x) as.numeric(colnames(W)[which(x>0)]))		
	
	remain.k <- as.numeric(colnames(W))
	
	while(length(remain.k)>0) {
		
		seed.k <- remain.k[which.max(w.sum[remain.k])]		
				
		v.t <- create_v(seed.k, w.sum, n.list, alpha)
		
		v.t.old <- list()
		
		if(verbose) 
			message(sprintf("\nForming new module with seed: %s", 
				original.id[seed.k]))
		
		# grow module v.t
		while(!identical(v.t, v.t.old)) {
			
			v.t.old <- v.t
			
			v.t <- grow_v(v.t, W, w.sum, n.list, 
				alpha, verbose, original.id)
		}
		
		V <- c(V, list(v.t))
		
		remain.k <- binary_delete(remain.k, seed.k)
		
		for(v in v.t$v.in) {
			remain.k <- binary_delete(remain.k, v)
		}				
	}
	
	if(verbose) {
		message("\nMRM before merging:")
		print_modules(V, FALSE, FALSE, original.id)
	}
	
	# Step 2: Merging clusters
	if(length(V)>1)
		V <- merge_v(V, W, w.sum, merge.tol, alpha, verbose)

			
	# Step 3: Discard clusters with low density
	V[sapply(V, function(v) v$density < density.tol)] <- NULL
	
	
	if(verbose) {
		message("\nMRM after merging/filtering:")
		print_modules(V, TRUE, FALSE, original.id)
	}		
	
	return(V)
}



# create a new module with seed k
create_v <- function(seed.k, w.sum, n.list, alpha) {
	
	list(v.in=seed.k,
		w.in=0, w.bound=w.sum[seed.k],
		penalty=alpha, synergy=0, density=0, card=1,
		v.bound = n.list[[seed.k]])
}


# try add one node at a time
grow_v_helper <- function(x, v.c, W, w.sum, alpha, 
	action=c("add","remove")) {
			
	w.in.x <- sum(W[x, v.c$v.in])

	w.bound.x <- w.sum[x] - w.in.x	
	
	if(action=="add") {
		
		w.in.new <- v.c$w.in + w.in.x
		w.bound.new <- v.c$w.bound + w.bound.x - w.in.x
		penalty.new <- v.c$penalty + alpha
		
	} else if(action=="remove" & w.bound.x > 0){
		
		w.in.new <- v.c$w.in - w.in.x
		w.bound.new <- v.c$w.bound - w.bound.x + w.in.x	
		penalty.new <- v.c$penalty - alpha
		
	} else if(action=="remove" & w.bound.x == 0){
		
		# do not remove node that does not have
		# connection with any external nodes
		# setting w.in.new to zero will result in
		# zero synergy, so the node will not be removed
		w.in.new <- 0		
		w.bound.new <- 0
		penalty.new <- alpha
	}
			
	synergy.new <- w.in.new/(w.in.new + w.bound.new + penalty.new)
		
	list(node.x=x,
		action=action,
		w.in.new = w.in.new,
		w.bound.new = w.bound.new,
		penalty.new = penalty.new,				
		synergy.new = synergy.new
	)
}



# grow module v by deciding whether to add i or k based on synergy
grow_v <- function(v.c, W, w.sum, n.list, alpha,
	verbose, original.id, gencluster=FALSE, N) {
		
	# for gencluster (stage 2) only mRNA are 
	# allowed to be added/remove
	# check whether the smallest id of node.bound
	# is a mRNA first
	if(gencluster) {

		if(v.c$v.bound[1] <= N) {
			
			idx <- binary_search(v.c$v.bound, N)
			
			v.bound <- v.c$v.bound[1:idx]
		} else {
			v.bound <- NULL
		}
		
	} else {
		v.bound <- v.c$v.bound
	}
	
	# if(verbose) message(cat("consider adding: ", original.id[v.bound]))
	
	# try all boundary miRNA/mRNA to select 
	# the one that gives the highest synergy
	v.add.list <- lapply(v.bound, grow_v_helper, 
		v.c, W, w.sum, alpha, "add")
			
	if(length(v.add.list) > 0) {
		v.add <- v.add.list[[which.max(sapply(v.add.list, 
			function(x) x$synergy.new))]]
	} else {
		v.add <- list(synergy.new=0)
	}
	
	if(v.c$card > 1) {
		
		# check whether the smallest id of node.in
		# is a mRNA first
		if(gencluster) {# only for Stage II
			
			if(v.c$v.bound[1] <= N) {
				
				idx <- binary_search(v.c$v.in, N)
			
				v.in <- v.c$v.in[1:idx]
				
			} else {
				v.in <- NULL
			}
			
		} else {
			v.in <- v.c$v.in
		}
		
		v.rmv.list <- lapply(v.in, grow_v_helper, 
			v.c, W, w.sum, alpha, "remove")
		
		if(length(v.rmv.list) > 0) {
			v.rmv <- v.rmv.list[[which.max(sapply(v.rmv.list, 
				function(x) x$synergy.new))]]
		} else {
			v.rmv <- list(synergy.new=0)
		}
				
		if(v.add$synergy.new > v.rmv$synergy.new) 
			v.chg <- v.add
		else
			v.chg <- v.rmv
		
	} else if (v.add$synergy.new > 0){
		
		v.chg <- v.add
		
	} else {
		
		return(v.c)
	}
	
	
	if(v.chg$synergy.new > v.c$synergy) {		
				
		v.c$w.in <- v.chg$w.in.new
		
		v.c$w.bound <- v.chg$w.bound.new
		
		v.c$synergy <- v.chg$synergy.new
		
		v.c$penalty <- v.chg$penalty.new
		
		if(v.chg$action=="add") {			
			
			v.c$v.in <- binary_insert(v.c$v.in, v.chg$node.x)
			
			v.c$card <- v.c$card + 1
			
		} else if(v.chg$action=="remove") {						
			
			v.c$v.in <- binary_delete(v.c$v.in, v.chg$node.x)
			
			v.c$card <- v.c$card - 1
		}
		
		v.c <- update_boundary(v.c, v.chg$node.x, n.list, v.chg$action)			
		if(verbose) {
			message(sprintf(
			"%s: %s; w.in: %.3f; w.bound: %.3f; new synergy: %.3f", 
			v.chg$action, original.id[v.chg$node.x], 
			v.c$w.in, v.c$w.bound, v.chg$synergy.new))
		}				
	}	
	
	v.c
}


# update boundary nodes for the newly formed module
update_boundary <- function(v.c, x, n.list, action=c("add","remove"), 
	verbose=FALSE, original.id) {		
			
	x.n <- n.list[[x]]
	
	if(action=="add") {
		
		# binary insert
		v.c$v.in <- binary_insert(v.c$v.in, x)
		
		# binary remove
		v.c$v.bound <- binary_delete(v.c$v.bound, x)
		
		# the neighbour nodes of x that are not internal nodes
		# and not existing boundary nodes become new boundary nodes		
		for(n in x.n) {
			
			idx <- binary_search(v.c$v.in, n)
			
			if(idx == 0) {
				
				v.c$v.bound <- binary_insert(v.c$v.bound, n)
								
			} else if(v.c$v.in[idx] != n) {
				
				v.c$v.bound <- binary_insert(v.c$v.bound, n)
			}
		}
		
	} else if(action=="remove") {
		
		# binary remove
		v.c$v.in <- binary_delete(v.c$v.in, x)
		
		# binary insert
		v.c$v.bound <- binary_insert(v.c$v.bound, x)
		
		# boundary neighbours of x that is not connected to
		# any of the internal nodes should be removed from
		# boundary nodes		
		for(n in x.n) {
			
			for(v.in in v.c$v.in) {		
				
				v.in.n <- n.list[[v.in]]
				
				idx <- binary_search(v.in.n, n)
				
				if(idx==0) {
					
					v.c$v.bound <- binary_delete(v.c$v.bound, n)
					
				} else if(v.in.n[idx] != n) {
					
					v.c$v.bound <- binary_delete(v.c$v.bound, n)
				}				
			}
		}		
	}
	
	if(verbose & !missing(original.id)) {
		message("update boundary nodes:")
		cat(as.character(original.id[v.c$v.bound]),"\n")
	}
	
	v.c
}

# calculate synergy of module v.c
get_synergy <- function(v.c) {
	
	v.c$w.in/(v.c$w.in + v.c$w.bound + v.c$penalty)
}

# calculate density of module v.c
get_density <- function(v.c) {
		
	if(v.c$card > 1)
		2*v.c$w.in/(v.c$card * (v.c$card-1))
	else
		0
}


# get overlap score
get_omega <- function(v.a, v.b) {
	
	ovp <- sum(v.a$v.in %in% v.b$v.in)
	
	ovp^2/(v.a$card * v.b$card)
}


# merge two clusters if their overlap score (omega) above the tol
merge_v <- function(V, W, w.sum, merge.tol, alpha, verbose) {
					
	# C x C overlapping graph for C clusters
	O <- lapply(V, function(v) sapply(V, get_omega, v.b=v))
	
	O <- do.call("rbind", O)
	
	if(verbose) {
# 		message(sprintf("\nCluster overlap relations:\n"))
# 		print(O)
		message(sprintf("Merge cluster threshold: overlap >= %.2f:\n", merge.tol))
	}
	
	adjm <- (O >= merge.tol)
	
	# get connected compoent memberships	
	mem <- clusters(graph.adjacency(adjm, mode="upper"))$membership
	
	V <- lapply(unique(mem), function(cls) {				
		
		V.sel <- V[mem == cls]
		
		v.c <- list(v.in=
			sort(unique(unlist(lapply(V.sel, function(v) v$v.in)))))
		
		update_v(v.c, W, w.sum, alpha)
	})
	
	V
}

# update module info
update_v <- function(v.c, W, w.sum, alpha) {
	
	w.in0 <- W[v.c$v.in, v.c$v.in]
	
	w.in <- sum(w.in0[upper.tri(w.in0)])
	
	# edges connecting nodes within v.c to nodes outside v.c
	w.bound <- sum(w.sum[v.c$v.in]) - sum(w.in0)
	
	v.c$w.in <- w.in
	
	v.c$w.bound <- w.bound
	
	v.c$card <- length(v.c$v.in)
	
	v.c$penalty <- alpha * v.c$card			
		
	v.c$synergy <- get_synergy(v.c)
		
	v.c$density <- get_density(v.c)	
	
	v.c
}


# update boundary nodes for the newly formed module
# update_boundary_entirely <- function(v.c, W, verbose=FALSE, original.id) {
	
	# W.sel <- W[!rownames(W) %in% v.c$v.in, v.c$v.in, drop=FALSE]
	
	# if(nrow(W.sel)>0) {	
		
		# v.c$v.bound <- as.numeric(rownames(W.sel)[apply(W.sel>0, 1, any)])
	# }
	
	# if(nrow(W.sel)==0 || length(v.c$v.bound)==0) {
		
		# v.c["v.bound"] <- list(NULL)
	# }
			
	# if(verbose & !missing(original.id)) {
		# message("update boundary nodes:")
		# cat(as.character(original.id[v.c$v.bound]),"\n")
	# }
			
	# v.c			
# }








