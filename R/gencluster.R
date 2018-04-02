gencluster <- function(V, W, H, alpha, density.tol, verbose) {

        if(!isSymmetric(H)) {
                warning(paste("Gene-gene interaction must be symmetric!",
                "Upper triangle matrix of H is used.", sep="\n"))

                H[lower.tri(H)] <- t(H)[lower.tri(t(H))]

                diag(H) <- 0
        }

        if(class(H)!="dgCMatrix")
                H <- Matrix(as.matrix(H))

        N <- as.numeric(nrow(H))
        M <- as.numeric(ncol(W))

        # H: N x N
        # W: N x M
        # W2: (N+M) x (N+M)
        W2 <- rbind(cbind(H, W),
                cbind(t(W), Matrix(0,nrow=M, ncol=M)))

        # index mRNA & miRNA ID to enable numerical search
        # (more efficient than string search)
        original.id <- rownames(W2)

        names(original.id) <- 1:nrow(W2)

        dimnames(W2) <- list(1:nrow(W2), 1:ncol(W2))

        # precompute weights for downstream speed-up
        w.sum <- as.numeric(apply(W2,1,sum))

        # save neighbour nodes for each node as a list
        # for downstream speed-up
        tmp <- as.matrix(W2)

        n.list <- mclapply(as.numeric(rownames(W2)),
                function(i) as.numeric(which(tmp[i,]>0)))

        V2 <- mclapply(V, function(v.t) {

                if(verbose)
                        message(sprintf("\nForming new module with seed: %s",
                                paste(original.id[v.t$v.in+N],collapse=",")))

                v.t <- init_v(v.t, N, w.sum, n.list, alpha)

                mirID <- v.t$v.in # save original miRNA ID

                v.t.old <- list()

                m <- length(v.t$v.in)

                # grow module v.t
                while(!identical(v.t, v.t.old)) {

                        v.t.old <- v.t

                        v.t <- grow_v(v.t, W2, w.sum, n.list,
                                alpha, verbose, original.id, TRUE, N)
                }

                v.t$miRNA <- original.id[mirID]

                idx <- binary_search(v.t$v.in, N)

                v.t$mRNA <- original.id[v.t$v.in[1:idx]]

                v.t$v.in <- original.id[v.t$v.in]

                v.t$v.bound <- original.id[v.t$v.bound]

                v.t$card.m <- m # number of miRNA

                v.t$card.t <- v.t$card # number of target

                v.t$card <- v.t$card.m + v.t$card.t

                # update final module info
                v.t$density <- get_density2(v.t)

                v.t
        })

        # Step 3: Discard clusters with low density
        V2[sapply(V2, function(v) v$density < density.tol)] <- NULL

        V2
}


# create a new module with seed k
init_v <- function(v.c, N, w.sum, n.list, alpha) {

        # increment miRNA index by N
        v.in <- v.c$v.in + N

        list(v.in=v.in,
                w.in=0, w.bound=sum(w.sum[v.in]),
                synergy=0, density=0,
                # card=0, # only mRNA counts
                # penalty=0,
                penalty=v.c$penalty,
                card=v.c$card, # all count
                v.bound = sort(unique(unlist(n.list[v.in]))))
}


# calculate density of module v.c
get_density2 <- function(v.c) {

        v.c$w.in/(v.c$card.m * v.c$card.t + v.c$card.t * (v.c$card.t-1))
}
