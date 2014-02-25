# perform binary search of x on a sorted numerical array a
binary_search <- function(vec, x) {
		
	# .Internal(findInterval(vec, x, FALSE, FALSE))
	findInterval(x, vec, FALSE, FALSE)
}

binary_insert <- function(vec, x) {
		
	idx <- binary_search(vec, x)
	
	n <- length(vec)
	
	if(idx != 0) {
		
		# cannot insert existing item
		if(vec[idx] == x) return(vec)		
	}
	
	if(idx == 0) {
		
		vec <- c(x, vec)
		
	} else if(idx < n & idx > 0) {
		
		vec <- c(vec[1:idx],x,vec[(idx+1):n])
		
	} else if(idx == n) {
		
		vec <- c(vec, x)
	}
	
	return(vec)
}

binary_delete <- function(vec, x) {
		
	idx <- binary_search(vec, x)
	
	if(idx==0) {
		
		return(vec)
		
	} else if(vec[idx] == x) {
		
		return(vec[-idx])
		
	} else { # x is not in vec
		
		return(vec)
	}
}


# display each cluster by listing miRNA and mRNA IDs
print_modules <- function(V, show.density=FALSE, order=TRUE, original.id) {
	
	names(V) <- paste("M",1:length(V),sep="")
	
	print_helper <- function(m, original.id) {
		
		v <- V[[m]]
		
		v.in <- v$v.in
		
		if(order) v.in <- v.in[order(v.in)]
		
		if(show.density)
			cat(sprintf("%s (density=%s; synergy=%s): \n", 
				m, scientific(v$density), scientific(v$synergy)))
		else
			cat(sprintf("%s (synergy=%s): \n", m, v$synergy))
		
		if(!missing(original.id)) {
			cat(original.id[v.in], "\n")
		} else {
			cat(v.in, "\n")
		}		
	}	
	
	if(!missing(original.id))
		tmp <- lapply(names(V), print_helper, original.id)
	else
		tmp <- lapply(names(V), print_helper)
}


# display each cluster by listing miRNA and mRNA IDs
print_modules2 <- function(V) {
	
	names(V) <- paste("M",1:length(V),sep="")
	
	print_helper <- function(m, original.id) {
		
		v <- V[[m]]
		
		miRNA <- v$miRNA
		mRNA <- v$mRNA
				
		cat(sprintf("%s (density=%s; synergy=%s): \n", 
			m, scientific(v$density), scientific(v$synergy)))
		
		cat(miRNA, "\n")
		cat(mRNA, "\n")		
	}	
	
	tmp <- lapply(names(V), print_helper)
}



# create graph object from data.frame
# with vertex as mRNA/miRNA and edges as edge weights
# NB: for small network only
# For large network, please use our Mirsynergy Cytoscape plugin
create_graph <- function(V, W, H) {
	
	W <- as.matrix(W)
	H <- as.matrix(H[rownames(W),rownames(W)])
	
	N <- as.numeric(nrow(H))
	M <- as.numeric(ncol(W))	
	
	# H: N x N
	# W: N x M
	# W: (N+M) x (N+M)
	W <- rbind(cbind(H, W), 
		cbind(t(W), matrix(0,nrow=M,ncol=M)))
	
	if(is.null(names(V)))
		names(V) <- 1:length(V)
	
	W <- as.matrix(W)
	
	all.id <- rownames(W)
	
	df.v <- data.frame(name=all.id, cluster=NA)	
	
	for(i in names(V)) {
		
		v <- V[[i]]
		
		sel <- df.v$name %in% c(v$miRNA,v$mRNA)
		
		df.v[sel,"cluster"] <- paste(df.v[sel,"cluster"],i,sep=",")
	}
	
	df.v[,"cluster"] <- sub("NA,","",df.v[,"cluster"])
		
	df.e <- melt(W)[,c(2,1,3)]	
	
	colnames(df.e) <- c("from", "to", "weight")
	
	df.e <- df.e[df.e$weight>0,]
	
	g <- graph.data.frame(df.e, directed=FALSE, vertices=df.v)
	
	g
}

# plot graph object
plot_modules <- function(V, W, H, legend.pos="topright", ...) {	
  
	g <- create_graph(V,W,H)
	
	df.v <- get.data.frame(g,"vertices")
	
	df.e <- get.data.frame(g,"edge")
	
	cls <- sort(unique(df.v$cluster))

	if(length(cls) > 12)
		vertex.color <- rainbow(length(cls))
	else if(length(cls) > 2)
		vertex.color <- brewer.pal(length(cls),'Set3')
	else
		vertex.color <- brewer.pal(3,'Set3')[1:length(cls)]
		
	names(vertex.color) <- cls
		
	vertext.shape <- rep("square",nrow(df.e))
	
	vertext.shape[grep("miRNA", df.v$name)] <- "circle"		
	
	par(mar=c(0,0,0,0))
	
	plot(g, vertex.color=vertex.color[df.v$cluster], 
		vertex.shape=vertext.shape,
		# edge.label = round(df.e$weight,2),
		edge.width = df.e$weight,
		edge.arrow.size=0.5, ...)
		
	legend(legend.pos, legend=names(vertex.color), 
		fill=vertex.color[names(vertex.color)])
}


# summarize each module
summary_modules <- function(V) {
		
	df <- do.call("rbind",lapply(V, function(v) 
		data.frame(miRNA=length(v$miRNA), mRNA=length(v$mRNA),
		total=v$card, synergy=v$synergy, density=v$density)))
			
	myrle <- rle(df$miRNA[order(df$miRNA)])
	
	miRNA.cnt <- data.frame(modules=myrle$lengths, miRNA=myrle$values)
	
	myrle <- rle(df$mRNA[order(df$mRNA)])
	
	mRNA.cnt <- data.frame(modules=myrle$lengths, mRNA=myrle$values)	
											
	list(moduleSummaryInfo=df, 
		miRNA.internal=miRNA.cnt, mRNA.internal=mRNA.cnt)
}


# plot module size
plot_module_summary <- function(V) {    
  
	mylist <- summary_modules(V)
			
	df <- mylist$moduleSummaryInfo
		
	df$all <- df$miRNA + df$mRNA
	
	# get around with the warnings in package build
	miRNA <- ""; modules=""; mRNA=""; synergy=0; score=0;	
	
	grid.arrange(	
		ggplot(mylist$miRNA.internal, aes(x=as.factor(miRNA), y=modules)) + 
			xlab("Number of miRNA") + ylab("Frequency of modules") +
			geom_bar(fill="white", color="black", stat="identity") + 
			geom_text(aes(x=as.factor(miRNA), label=modules)) +
			theme_bw(),
			
# 		ggplot(mylist$mRNA.internal, aes(x=as.factor(mRNA), y=modules)) + 
# 			xlab("Number of mRNA") + ylab("Frequency of modules") +
# 			geom_histogram(fill="white", color="black", stat="identity") + 
# 		  geom_text(aes(x=as.factor(mRNA), label=modules)) +
# 			theme_bw(),
    
		ggplot(mylist$mRNA.internal, aes(x=mRNA)) + 
		  xlab("Number of mRNA") + ylab("Frequency of modules") +
		  geom_histogram(fill="white", color="black") + 		  
		  theme_bw(),	
			
		ggplot(df, aes(x=synergy)) + geom_density() + 
			xlab("Synergy score") + ylab("Density") + theme_bw(),
	
		ggplot(df, aes(x=miRNA, y=synergy)) + geom_point() + 
			xlab("Number of miRNA") + ylab("Synergy score") + 
      stat_smooth(method="lm", se=FALSE) + theme_bw()      
	)	
}


tabular_module_helper <- function(i, V, W, H) {	
  
	v <- V[[i]]			
	
	mmi <- W[v$mRNA, v$miRNA, drop=FALSE]
	
	mmi <- melt(mmi)[c(2,1,3)]
	
	colnames(mmi) <- c("A","B","score")
	
	mmi$interaction <- "MMI"
		
	ggi <- as.matrix(H[v$mRNA, v$mRNA, drop=FALSE])
	
	ggi <- as.data.frame(melt(ggi))
	
	colnames(ggi) <- c("A", "B", "score")
	
	ggi$interaction <- "GGI"
		
	edges <- data.frame(rbind(mmi,ggi), module=i)
	
	score <- 0 # get around with package build warning
	
	edges <- subset(edges, abs(score) >0)
		
	nodes <- rbind(data.frame(node=mmi$A, type="miRNA", module=i),
		data.frame(node=ggi$A, type="mRNA", module=i))		
	
	list(edges, nodes)
}


# tabular module miRNA and mRNA
tabular_module <- function(V, W, H, outdir) {
	
	mylist <- lapply(1:length(V), tabular_module_helper, V, W, H)

	edges <- do.call("rbind", lapply(mylist, function(x) x[[1]]))

	nodes <- do.call("rbind", lapply(mylist, function(x) x[[2]]))
			
	if(!missing(outdir)) {
		
		if(!file.exists(outdir)) system(sprintf("mkdir %s", outdir))
		
		write.table(nodes, sprintf("%s/module_nodes.txt",outdir), 
			quote=FALSE, sep='\t', row.names=FALSE)
		
		write.table(edges, sprintf("%s/module_edges.txt",outdir), 
			quote=FALSE, sep='\t', row.names=FALSE)		
	}
	
	list(nodes=nodes, edges=edges)
}

