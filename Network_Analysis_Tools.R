#
# A module of tools to assist in analyzing networks using igraph.
# 
# Justin W. Kenney
# jkenney9a@gmail.com
# 

library(igraph) #For network functions
library(Hmisc) #For generating correlation/p-value matrices
library(boot) #For bootstrapping
library(data.table) #for rbindlist
library(proxy) #For distance measures (e.g, Jaccard distance)
library(statGraph) #For Jensen-Shannon divergence between two graphs

load_data <- function(csv_file){
  #   Input: CSV file name
  #   
  #   Output: Dataframe containing CSV file
  
  df <- read.csv(file=csv_file, header=TRUE, check.names=FALSE)
  return(df)
}

clean_data <- function(df_count, real_zeros=TRUE, missing_thresh=4, fill_missing=TRUE){
  # Input: dataframe of counts to be cleaned, whether zeros are real data points of should be 
  # purged, the maximum threshold of missing values allowed otherwise column is removed
  # 
  #
  # Output: dataframe of counts with columns with more zeros than 
  # zero_threshold removed. Other zeros are replaced with column mean.
  
  #Replace zeros in data with NA to allow easier handling
  if (real_zeros==FALSE){
    df_count[mapply("==", df_count, 0)] <- NA
  }
    
  #Remove columns with more than missing_thresh missing values
  df_out <- df_count[,which(colSums(is.na(df_count)) < missing_thresh)]
    
  if (fill_missing==TRUE){
    col_avgs <- colMeans(df_out, na.rm=TRUE)
    index <- which(is.na(df_out), arr.ind=TRUE)
    df_out[index] <- col_avgs[index[,2]]
  }
  
  return(df_out)
}

normalize_data <- function(df_data, df_norm){
  # Input: dataframe, df, to be normalized and 
  # dataframe to normalize to (df_norm)
  #
  # Output: dataframe with each column divided by the average of the 
  # columns in df_norm. 
  #
  # NOTE: Inputs should already be "cleaned up" the same way using
  # functions above. If columns are not present in both df and df_norm, 
  # they are excluded
  
  #Remove columns from normalization df not in data df and vice versa:
  common.names <- intersect(colnames(df_data), colnames(df_norm))
  df_data <- df_data[,common.names]
  df_norm <- df_norm[,common.names]
  
  #Make sure columns of both dfs are in the same order:
  df_norm <- df_norm[colnames(df_data)]
  
  norm_avgs <- colMeans(df_norm)
  
  df_out <- mapply("/", df_data, norm_avgs)

  return(df_out)
}


corr_matrix <- function(df, p.adjust.method='none', type='pearson'){
  #   Input: Dataframe with headers as titles (brain regions)
  #           Whether or not to apply a p-value adjustment method drawn from the 'p.adjust' function
  #           (e.g, 'fdr', 'bonferroni' etc.). 
  #           whether to use 'pearson' or 'spearman' correlation
  #           NOTE: assume undirected graph!
  #   
  #   Output: List of two dataframes corresponding to 1) all pairwise Pearson 
  #   correlations and 2) all associated un-adjusted p-values


  corr <- rcorr(as.matrix(df), type=type)
  df_corr <- as.data.frame(corr['r'])
  
  #adjust p-values if necessary
  adjusted.P <- corr[['P']]
  adjusted.P[upper.tri(adjusted.P)] <- p.adjust(adjusted.P[upper.tri(adjusted.P)], method=p.adjust.method)
  adjusted.P[lower.tri(adjusted.P)] <- NA
  df_pvalue <- as.data.frame(adjusted.P)
  
  
  names(df_corr) <- gsub("^r.","",colnames(df_corr), fixed=FALSE) #Remove "r." from node names

  
  #Change ... to -; for some reason when R imports "-" it turns it into "..."
  names(df_corr) <- gsub("...","-",colnames(df_corr), fixed=TRUE)
  rownames(df_corr) <- gsub("...", "-", rownames(df_corr), fixed=TRUE)
  names(df_corr) <- gsub("sums.", "", colnames(df_corr), fixed=TRUE)
  rownames(df_corr) <- gsub("sums.", "", rownames(df_corr), fixed=TRUE)
  
  names(df_corr) <- rownames(df_corr)
  
  return(list("corr" = df_corr,"pvalue" = df_pvalue))
}

corr_matrix_threshold <- function(df, negs = FALSE, thresh=0.01, thresh.param='p', p.adjust.method='none',
                                  type='pearson'){
  #   Input: Dataframe with headers as titles (brain regions) of counts etc.
  #           threshold (can be a list of thresholds now) and threshold parameter (p, r, or cost), whether or not to keep negative correlations.
  #           and the p-value adjustment method if using p-value as threshold parameter.
  #           and whether to use pearson or spearman correlation
  #   
  #   Output: Dataframe of correlations thresholded at p < threshold
  #   
  #   NOTE: Removes diagonals
  
  dfs <- corr_matrix(df, p.adjust.method=p.adjust.method, type=type)
  df_corr <- dfs[['corr']]
  df_pvalue <- dfs[['pvalue']]
  
  #remove diagonals (may not be necessary when using igraph...)
  diag(df_corr) <- 0
  diag(df_pvalue) <- NA
  #df_corr[mapply("==", df_corr, 1)] <- 0
  
  #remove any NaNs, infs or NAs (sometimes happens with bootstrapping; not sure why)
  df_pvalue[mapply(is.infinite, df_corr)] <- 1
  df_pvalue[mapply(is.nan, df_corr)] <- 1  
  df_pvalue[mapply(is.na, df_pvalue)] <- 1
  df_corr[mapply(is.infinite, df_corr)] <- 0
  df_corr[mapply(is.nan, df_corr)] <- 0
  df_corr[mapply(is.na, df_corr)] <- 0

  #remove negative correlations
  if(negs == FALSE){
    df_corr[mapply("<", df_corr, 0)] <- 0
  }
  
  if(tolower(thresh.param)=='p'){
    #apply p-value threshold to correlation matrix
    df_corr <- lapply(thresh, function(x) {df_corr[mapply(">=", df_pvalue, x)] <- 0; df_corr})    
  } else if(tolower(thresh.param)=='r'){
    df_corr <- lapply(thresh, function(x) {df_corr[mapply("<=", abs(df_corr), x)] <- 0; df_corr})
  } else if(tolower(thresh.param)=='cost'){
    r.threshold <- quantile(abs(df_corr), probs=1-thresh, na.rm=TRUE)
    df_corr <- lapply(r.threshold, function(x) {df_corr[mapply("<=", abs(df_corr), x)] <- 0; df_corr})
  } else{
    stop("Invalid thresholding parameter")
  }
  
  if(length(thresh) == 1){
    df_corr <- df_corr[[1]]
  }
  
  return(df_corr)
}

threshold_matrix <- function(corr.mat, p.value.mat=NULL, thresh, thresh.param){
  #Input a correlation matrix, a p.value matrix (thresh param=p), a threshold, and a threshold parameter
  # e.g, (r, p, or cost)
  #
  #Output: the correlation matrix thresholded as appropriate
  
  
  
  if(tolower(thresh.param)=='p'){
    #apply p-value threshold to correlation matrix
    corr.mat[mapply(">=", p.value.mat, thresh)] <- 0    
  } else if(tolower(thresh.param)=='r'){
    corr.mat[mapply("<=", abs(corr.mat), thresh)] <- 0
  } else if(tolower(thresh.param)=='cost'){
    r.threshold <- quantile(abs(corr.mat), probs=1-thresh, na.rm=TRUE)
    corr.mat[mapply("<=", abs(corr.mat), r.threshold)] <- 0
  } else{
    stop("Invalid thresholding parameter")
  }
  
  return(corr.mat)
}

CSV_to_igraph <- function(CSV_file, negs = FALSE, thresh=0.01, thresh.param='p', p.adjust.method='none',
                          type='pearson'){
  #Input: CSV_file of counts etc., threshold, whether threshold based on p-value ('p')
  # or correlation(r) value ('r') or cost ('cost'). The p-value adjustment method to use (if any) and 
  # the type of correlation to perform (pearson or spearman)
  #
  #Output: igraph graph object
  
  df <- load_data(CSV_file)
  df_G <- corr_matrix_threshold(df, negs = negs, thresh=thresh, thresh.param=thresh.param, 
                                p.adjust.method=p.adjust.method, type=type)
  
  G <- graph.adjacency(as.matrix(df_G), mode="undirected",
                       weighted=TRUE)
  
  return(G)
}

df_to_igraph <- function(df, negs=FALSE, thresh=0.01, thresh.param='p', p.adjust.method='none',
                         type='pearson'){
  #Input: df of counts per brain region, whether to keep negative
  #correlations and the threshold (r, p, or cost), the p-value adjustment method to use (if any) and
  # the type of correlation (pearson or spearman)
  #
  #Output: igraph graph
  
  df_G <- corr_matrix_threshold(df, negs = negs, thresh=thresh, thresh.param=thresh.param, 
                                p.adjust.method=p.adjust.method, type=type)
  
  #Handle case of single threshold
  if(length(thresh) == 1){
    df_G <- list(df_G)
  }

  G <- lapply(df_G, function(x) {graph.adjacency(as.matrix(x), mode="undirected", 
                       weighted=TRUE, diag=FALSE)})
  #Handle case of single threshold
  if(length(thresh)==1){
    G <- G[[1]]
  }
  
  return(G)
}

get_centrality_measures <-function(G, weighted=FALSE, nodal_efficiency=FALSE, normalized=FALSE, min_max_normalization=FALSE){
  # Input: An igraph graph
  #
  # Output: A dataframe of centrality measures:
  # Degree, betweenness, eigenvector, closeness, 
  # transitivity (~clustering coefficient), nodal efficiency
  # NOTE: nodal efficiency makes use of weights if they exist
  # but has been removed for now b/c of how long the calculation can take!
  if(ecount(G) == 0){
    zeros <- rep(0, vcount(G))
    return(data.frame("degree"=zeros, "betweenness"=zeros, 
                      "eigenvector"=zeros, "closeness"=zeros,
                      "transitivity"=zeros,
                      row.names=V(G)$name))
  }
    
  degree <- igraph::degree(G, normalized=normalized)
  G.pos <- G
  E(G.pos)$weight <- abs(E(G.pos)) #Positive numbers are necessary for betweenness and closeness
  
  if(weighted==FALSE){
    eigenvector <- evcent(G, weights=NA)
    between <- betweenness(G.pos, normalized=normalized, weights=NA)
    close <- closeness(G.pos, weights=NA, normalized=normalized)
  } else{
      eigenvector <- evcent(G)
      between <- betweenness(G.pos, normalized=normalized, weights=NULL)
      close <- closeness(G.pos, normalized=normalized)
  }
  
  trans <- transitivity(G,type="local", isolates='zero')
  
  
  #Need to pull out matrices for efficiency calculations
  #adj_mat <- as.matrix(get.adjacency(G))
  #weight_mat <- as.matrix(get.adjacency(G, attr='weight'))
  #efficiency <- Global_efficiency(G)
  
  output <- data.frame("degree" = degree, "betweenness" = between, 
                       "eigenvector" = eigenvector[[1]], "closeness" = close,
                      "transitivity" = trans)
  
  if(nodal_efficiency==TRUE){
    node_efficiency <- Nodal_efficiency(G, normalized = normalized)
    output$efficiency <- node_efficiency
  }
  
  if(min_max_normalization==TRUE){
    output <- apply(output, MARGIN=2, function(x) {(x-min(x)) / (max(x) - min(x))})
    output <- as.data.frame(output)
    output$degree.betweenness <- output$degree + output$betweenness
    output$transitivity <- trans #Do not normalize transitivity b/c it is already normalized
  }
  
  
  return(output)
}

GC_size <- function(G){
  #Input: igraph graph
  #
  #Output: size of giant component
  
  return(max(clusters(G)[[2]]))
}

get_GC <- function(G){
  clust <- clusters(G)
  GC <- induced.subgraph(G, which(clust$membership == which.max(clust$csize)))
  return(GC)
}

Global_efficiency <- function(G, weighted=TRUE){
  #Input: igraph graph
  #
  #Output: efficiency; 
  #NOTE: default uses weighted matrix and assumes 
  #the attribute is called "weight" in the graph.
  #NOTE: The absolute value of edges is used as negative correlations
  # result in odd behaviour. From a theoretical perspective I believe this is 
  # ok b/c a neg. correlation tells us the same thing as a positive correlation
  # with respect to the transfer of information across a network
  #
  #NOTE (16/03/16): Altered calculation for weighted networks. We assume a 
  # weighted correlation network. The way the shortest path-length algorithm works
  # is that lower numbers = shorter distances. However, with correlations the 
  # opposite is true. To correct for this, will take the inverse of the weights. This 
  # works out to be the equivalent of taking the harmonic mean for the path length part
  # of the calculation.
  

  if(weighted==TRUE & ecount(G) > 0){
    E(G)$weight <- 1 / abs(E(G)$weight)
    eff <- 1/shortest.paths(G)
  }else{
    eff <- 1/shortest.paths(G, weights=NA)
  }
  
  eff[!is.finite(eff)] <- 0
  gl.eff <- mean(eff[upper.tri(eff)])
  
  
  return(gl.eff)
}

get_minimum_connectivity_threshold <- function(mat, densities=seq(1,0.25,-0.01)){
  #Generates a thresholded correlation/connectivity matrix in which all nodes have at 
  #least one connection. 
  #
  #Input: A corrleation matrix, vector of densities to examine
  #
  #Output: a thresholded correlation matrix based on the above inidicated condition. 
  
  #Make sure we have a matrix (in case we get dataframe input)
  mat <- as.matrix(mat)
  diag(mat) <- 0
  
  
  for(cost in densities){
    thresh <- quantile(abs(mat), probs=cost)
    mat.test <- mat
    mat.test[abs(mat) < thresh] <- 0
    
    if(all(colSums(abs(mat.test)) > 0)) break
  }
  
  return(thresh)
  
}

rich.club.coeff <- function(G, k, weighted=FALSE){
  #Calculate the rich club coefficient for graph G given cut-off k
  
  deg <- degree(G)
  deg <- sort(deg, index.return=TRUE)
  G.sub <- induced_subgraph(G, v=deg$ix[deg$x > k])
  phi <- (2*ecount(G.sub))/(vcount(G.sub) * (vcount(G.sub) - 1))
  
  return(list('phi'=phi, 'graph'=G.sub, 'Nk'=vcount(G.sub), 'Ek'=ecount(G.sub)))
}

participation.coeff <- function(G, Comm){
  #Calculate the participation coefficient of all the nodes in a graph (G) given a specific
  # igraph community structure (Comm). This code was taken largely from the brainGraph package
  
  membership <- Comm$membership
  deg <- degree(G)
  nodes <- which(deg > 0)
  part.coeffs <- rep(0, length(deg))
  
  for(n in nodes){
    Kis <- vapply(seq_len(max(membership)), function(x) sum(neighbors(G, n) %in% which(membership == x)),
                  FUN.VALUE=integer(1))
    Ki <- deg[n]
    part.coeffs[n] <- 1 - sum((Kis/Ki)^2)
  }
  
  return(part.coeffs)
}

within.module.deg.z.score <- function(G, Comm){
  #Calculate the within module degree z scores for all the nodes in a graph (G) given a specific
  # igraph community structure (Comm). This code was taken largely from the brainGraph package
  
  membership <- Comm$membership
  deg <- degree(G)
  nodes <- which(deg > 0)
  edges <- E(G)
  z <- Ki <- rep(0, length(deg))
   
  for(n in nodes){
    Ki[n] <- length(edges[n %--% which(membership == membership[n])])
   }
  
  di <- lapply(seq_len(max(membership)), function(x) Ki[membership == x])
  Ksi <- vapply(di, mean, FUN.VALUE=numeric(1))
  sigKsi <- vapply(di, sd, numeric(1))
  z[nodes] <- (Ki[nodes] - Ksi[membership[nodes]])/sigKsi[membership[nodes]]
  z <- ifelse(!is.finite(z), 0, z)
  return(z)
}

Nodal_efficiency <- function(G, weighted=FALSE, normalized=FALSE) {
  #Input: igraph graph, and whether or not to consider weights in the efficiency calculation
  # whether or not to normalize (where max value = 1)
  
  #Output: data frame of nodal efficiency (i.e, average inverse shortest path length from individual node to all other nodes)
  
  if(weighted == TRUE & ecount(G) > 0){
    E(G)$weight <- 1 / abs(E(G)$weight) #Take inverse of weights b/c dealing with correlations
    eff <- 1 / shortest.paths(G)
  }else{
    eff <- 1/shortest.paths(G, weights=NA)
  }
  
  eff[!is.finite(eff)] <- 0
  out <- colSums(eff) / (vcount(G) - 1)
  
  if(normalized == TRUE){
    out <- out / max(out)
  }
  
  return(out)
}


node.distance.comparison <- function(graph.list, method='jaccard', weighted=FALSE){
  #This function computes the distance between nodes of the given network based on their
  # connections
  
  #Input: a list of 2 graphs or adjacency matrices to compare nodes across, and the metric to use
  #
  #Output: a dataframe with the node and distance between nodes in the two dataframes
  if(!is.list(graph.list)){
    stop('Input must be in list form')
  }
  if(length(graph.list) != 2){
    stop('Can only compare nodes between two graphs. More than two graphs have been given.')
  }
  
  if(is.igraph(graph.list[[1]])){
    if(weighted==TRUE){graph.list <- lapply(graph.list, get.adjacency, type='both', attr='weight', sparse=FALSE)}
    else{graph.list <- lapply(graph.list, get.adjacency, type='both', sparse=FALSE)}
  }
  
  #Make sure we have only 0's and 1's in matrix
  if(weighted==FALSE){graph.list <- lapply(graph.list, function(x) {x[x!=0] <- 1; x})}
  
  df.distance <- dist(graph.list[[1]], graph.list[[2]], method = method, pairwise=TRUE)
  df.distance <- data.frame(node = colnames(graph.list[[1]]), distance = as.numeric(df.distance))
  
  return(df.distance)
}


graph_measures_across_threshs <- function(df.counts, thresh.list, thresh.param='cost', negs=TRUE, small.world=TRUE,
                                          small.world.iterations=100, p.adjust.method='none', rand.degree.dist=TRUE,
                                          community=NA){
  
  #thresh.params = cost, r, or p
  #Only need to supply p.value.mat is using thresh.param=p
  
  # diag(corr.mat) <- 0
  # corr.mat <- abs(corr.mat)
  # corr.mat[is.na(corr.mat)] <- 0
  # 
  # 
  # if (thresh.param=='cost'){
  #   l.thresholds <- quantile(corr.mat[upper.tri(corr.mat)], probs=1-thresh.list, na.rm = TRUE)
  # }else if (thresh.param=='r'){
  #   l.thresholds <- thresh.list
  # } 
  # 
  # if (thresh.param=='cost' | thresh.param=='r'){
  #   l.mat <- lapply(l.thresholds, function(x) {corr.mat[corr.mat < x] <- 0; corr.mat})
  # }else if (thresh.param=='p'){
  #   l.mat <- lapply(thresh.list, function(x) {corr.mat[p.value.mat > x] <- 0; corr.mat})
  #   l.r.thresh <- lapply(thresh.list, function(x) {min(corr.mat[p.value.mat > x] <- 0)})
  # }
  
  l.mat <- corr_matrix_threshold(df.counts, negs=negs, thresh=thresh.list, thresh.param=thresh.param, p.adjust.method = p.adjust.method)
  l.mat <- lapply(l.mat, abs)
  
  l.G <- lapply(l.mat, function(x) {G <- graph.adjacency(as.matrix(x), mode='undirected', weighted = TRUE, diag=FALSE); G})
  
  l.GE <- as.numeric(lapply(l.G, Global_efficiency, weighted=FALSE))
  l.trans <- as.numeric(lapply(l.G, transitivity, type='global'))
  
  
  
  #% nodes in connected componenet
  l.GC.sizes <- lapply(l.G, GC_size)
  l.GC.percent <- as.numeric(lapply(l.G, function(x) {(GC_size(x)/vcount(x)) * 100}))
  
  l.density <- unlist(lapply(l.G, graph.density))
  
  #Modularity using Louvian algorithm
  if(is.na(community)){
    l.G.comm <- lapply(l.G, cluster_louvain)
    l.Q <- as.numeric(mapply(function(x,y) {modularity(x,y$membership)}, l.G, l.G.comm, SIMPLIFY=FALSE))
  } else {
    l.Q <- as.numeric(lapply(l.G, function(x) {modularity(x, community$membership)}))
  }
    
  df.out <- data.frame(density=l.density, GE = l.GE, clustering = l.trans, percent.connected = l.GC.percent, 
                       modularity=l.Q, threshold=thresh.list)
  
  if (thresh.param == 'p'){
    l.r.thresh <- lapply(l.mat, function(x) {min(x[x != 0])})
    df.out$r.thresh <- l.r.thresh
  }
  #Small worldness calculations/comparisons to random graphs
  if(small.world==TRUE){
    #l.rand.graphs <- lapply(l.G, Rand_graph_stats, degree_dist=rand.degree.dist, weighted=FALSE, iterations=small.world.iterations)
    l.rand.graphs <- lapply(l.G, function(x) {graph.density(x); out <- Rand_graph_stats(x, degree_dist=rand.degree.dist,
                                                                                        weighted=FALSE, iterations=small.world.iterations); out})
    l.lambda <- mapply(function(x,y) {y$Global.efficiency[1] / x}, l.GE, l.rand.graphs, SIMPLIFY=TRUE)
    l.gamma <- mapply(function(x,y) {x / y$Transitivity[1]}, l.trans, l.rand.graphs, SIMPLIFY=TRUE)
    l.sigma <- mapply(function(x,y) {x/y}, l.gamma, l.lambda, SIMPLIFY=TRUE)
    df.out$lambda <- l.lambda
    df.out$gamma <- l.gamma
    df.out$sigma <- l.sigma
  }
  
  return(df.out)
}


centrality_measures_across_threshs <- function(df.counts, thresh.list, negs=TRUE, normalized=TRUE, thresh.param='cost', p.adjust.method='none',
                                               weighted=FALSE, community=NA){
  
  l.mat <- corr_matrix_threshold(df.counts, negs=negs, thresh=thresh.list, thresh.param=thresh.param, p.adjust.method = p.adjust.method)
  l.mat <- lapply(l.mat, abs)
  
  l.G <- lapply(l.mat, function(x) {G <- graph.adjacency(as.matrix(x), mode='undirected', weighted = TRUE, diag=FALSE); G})
  
  
  
  # diag(corr.mat) <- 0
  # corr.mat <- abs(corr.mat)
  # corr.mat[is.na(corr.mat)] <- 0
  # 
  # if (thresh.param=='cost'){
  #   l.cost.thresholds <- quantile(corr.mat[upper.tri(corr.mat)], probs=1-cost.list, na.rm = TRUE)
  # }else if (thresh.param=='r'){
  #   l.cost.thresholds <- cost.list
  # } 
  # 
  # l.mat <- lapply(l.cost.thresholds, function(x) {corr.mat[corr.mat < x] <- 0; corr.mat})
  # l.G <- lapply(l.mat, function(x) {G <- graph.adjacency(as.matrix(x), mode='undirected', weighted = TRUE, diag=FALSE); G})
  
  if(is.na(community)){
    l.G.comm <- lapply(l.G, cluster_louvain)
    l.part.coef <- mapply(participation.coeff, l.G, l.G.comm, SIMPLIFY=FALSE)
    l.within.mod.z.score <- mapply(within.module.deg.z.score, l.G, l.G.comm, SIMPLIFY=FALSE)
  } else {
    l.part.coef <- lapply(l.G, participation.coeff, community)
    l.within.mod.z.score <- lapply(l.G, within.module.deg.z.score, community)
  }
  
  l.centrality <- lapply(l.G, get_centrality_measures, nodal_efficiency=TRUE, normalized=normalized, weighted=weighted, min_max_normalization=TRUE)
  l.centrality <- mapply(function(x,y) {x$participation.coefficient <- y; x}, l.centrality, l.part.coef, SIMPLIFY=FALSE)
  l.centrality <- mapply(function(x,y) {x$within.module.z.score <- y; x}, l.centrality, l.within.mod.z.score, SIMPLIFY=FALSE)
  l.centrality <- mapply(function(x,y) {x$threshold <- y; x}, l.centrality, thresh.list, SIMPLIFY=FALSE)
  l.centrality <- lapply(l.centrality, function(x) {x$region <- row.names(x) ; x})
  df.centrality <- rbindlist(l.centrality)
  
  
  
  return(df.centrality)
}


Rand_graph_stats <- function(G, iterations = 100, degree_dist=TRUE, weighted=FALSE, rich_club=FALSE, rich_club_k=NULL,
                             rich_club_weighted=FALSE){
  #Input: igraph graph, wehther or not random graphs should have
  #same degree distribution as input graph. If degree_dist=FALSE
  #uses Erdos-Renyi random graph. If rich_club = TRUE, then you must supply a rich_club_k,
  #which can be a single number or a range of numbers. Rich club analysis can only be done if the
  #degree distribution = TRUE
  #
  #Output: If no rich club analysis, then equivalent random graph statistics and stdev in a df:
  #Transitivity (very closely related to clustering coefficient)
  #Global efficiency
  #
  #If rich club analysis, then outputs a list where the first item is the dataframe containing the 
  #info about the transitivity and global efficiency, and the second item is a dataframe of the average 
  #rich club coefficient at each level of k provided
  
  TR <- c() #Transitivity
  GE <- c() #global efficiency (GE)
  path_length <- c() #Path length
  phi <- list() #rich club coefficient
  
  if(degree_dist == TRUE){
    #Get rid of degree zero nodes b/c this causes issues with the
    #generation of random graphs using the "vl" algorithm
    G_1 <- delete.vertices(G, which(igraph::degree(G) == 0))

    degrees <- igraph::degree(G_1)
    
    #Get number of zero degree nodes
    zero_d <- length(igraph::degree(G)) - length(degrees)
    
    l.RG <- lapply(c(1:iterations), function(x) tryCatch(degree.sequence.game(igraph::degree(G_1), method='vl'), error=function(cond) {return(graph.empty(0, directed=FALSE))}))
    TR <- sapply(l.RG, transitivity, type='global')
    GE <- sapply(l.RG, function(x) Global_efficiency(x + zero_d, weighted=FALSE)) #add back in zero deg. nodes
    path_length <- sapply(l.RG, average.path.length)
    phi <- lapply(l.RG, function(x) {lapply(rich_club_k, function(y) rich.club.coeff(x, k=y, weighted=rich_club_weighted)$phi)})
    #process rich club stats
    phi <- rbindlist(phi)
    phi <- colMeans(phi) #get average for each iteration at each k value
    names(phi) <- rich_club_k
    } else{
  
  if(degree_dist == FALSE){
    
    if(rich_club == TRUE){
      stop('cannot run rich club analysis if degree distribution is not degree distribution matched')
    }
    
    n <- vcount(G) #number of nodes
    m <- ecount(G) #number of edges
    
    l.RG <- lapply(c(1:iterations), function(x) erdos.renyi.game(n, m, type='gnm'))
    TR <- sapply(l.RG, transitivity, type='global')
    GE <- sapply(l.RG, function(x) Global_efficiency(x, weighted=FALSE))
    path_length <- sapply(l.RG, average.path.length)
    
    }
    }
  
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  path_length_out <- c(mean(path_length), sd(path_length))
  df.out <- data.frame("Transitivity" = TR_out, 
                       "Global.efficiency" = GE_out,
                       "Path.length" = path_length_out,
                       row.names = c("Average", "stdev"))
  
  if(rich_club==TRUE){
    return(list('small.world'=df.out, 'rich.club'=phi))
  } else {
    return(df.out)
  }
}

Watts_Strogatz_model <- function(G, iterations=1000, trans_match_iter=100){
  #Input: A graph to use as basis for generating Watts-Strogatz small 
  # world model. Number of random Watts-Strogatz graphs to calculate and
  #the number of iterations to use for matching transitivity/clustering
  #coefficient
  #
  #Output: List containing three items:
  # 1) Random graph statistics and stdev in a df for transitivity and global efficiency
  # 2) The average sorted degree distribution for all the WS graphs generated
  # 3) The average sorted clustering distribution for all the WS graphs generated
  # 4) A dataframe containing all degree and transitivity values from every iteration
  #
  #NOTE: wiring probability for graph is set such that the clustering 
  #coefficient of the Watts-Strogatz graphs are approximately the same
  #as the input graph.
  
  n <- vcount(G) #number of nodes
  e <- ecount(G) #number of edges
  G.trans <- transitivity(G, type="global")
  
  nei <- round(e/n) #Calculate number of edges in each node neighborhood
                    #this measure attempts to get as close to the number of 
                    #edges as possible to the original graph, but unlikely to be exact
  
  wiring.ps <- seq(0,1,0.01) #generate re-wiring probabilities to test
  
  #Data frame to hold transitivity/clustering values 
  trans.out <- as.data.frame(matrix(NA, ncol=length(wiring.ps), nrow=trans_match_iter, 
                                    dimnames=list(1:trans_match_iter, wiring.ps)))
  
  #Generate a series of Watts-Strogatz graphs to find out what p-value matches the 
  #transitivity/clustering in the original graph
  for(i in 1:trans_match_iter){
    trans.out[i,] <- sapply(wiring.ps, function(x) 
                            transitivity(watts.strogatz.game(1,n,nei,x), type="global"))
  }
  
  #Get the p-value associated with the graph matching clustering/transitivity of original graph
  trans.out.means <- colMeans(trans.out)
  p.value <- wiring.ps[which(abs(trans.out.means-G.trans) == min(abs(trans.out.means-G.trans)))]
  
  TR <- c()     #Vector to hold transitivity calculations
  GE <- c()     #Vector to hold global efficiency calculations
  deg.dist <- as.data.frame(matrix(NA, ncol=n, nrow=iterations))
  clust.dist <- as.data.frame(matrix(NA, ncol=n, nrow=iterations))
  
  #Generate many graphs
  for(i in 1:iterations){
    W <- watts.strogatz.game(1,n,nei,p.value)
    TR <- append(TR, transitivity(W, type="global"))
    GE <- append(GE, Global_efficiency(W, weighted=FALSE))
    
    deg.dist[i,] <- sort(igraph::degree(W))
    clust.dist[i,] <- sort(transitivity(W, type="local", isolates='zero'))
  }
  
  df_all_out <- data.frame(degree = unlist(deg.dist), transitivity = unlist(clust.dist))
  
  deg.dist <- colMeans(deg.dist)
  clust.dist <- colMeans(clust.dist)
  
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  df_out <- data.frame("Transitivity" = TR_out, "Global efficiency" = GE_out, row.names=c("Average", "stdev"))

  return (list(df_out, deg.dist, clust.dist, df_all_out))
  
}

BA_model <- function(G, iterations=1000){
  #
  # Generate Barabasi-Alberts scale free models based on the input graph G
  #
  # Input: graph G of interest and the number of iterations of BA model
  # to generate
  #
  # Output: List containing three items:
  # 1) Random graph statistics and stdev in a df for transitivity and global efficiency
  # 2) A df containing all the degrees and clustering/transitivity over all the nodes, to be used
  # to generate a histogram.
  #
  # NOTE: BA graphs are generated such that the average edges from the BA graphs that are closest to
  # the input graph is approximately equal to the input graph. The number of iterations used for each
  # BA graph is scaled to generate this result.
  
  n <- vcount(G) #nodes in input graph
  e <- ecount(G) #edges in input graph
  G.density <- e / (((n*(n-1))/2))
  
  nei.out <- 0
  BA.test <- e
  BA.G.diff <- c()
  BA.edges <- c()
  
    
  TR <- c()     #Vector to hold transitivity calculations
  GE <- c()     #Vector to hold global efficiency calculations
  deg.dist <- c()
  clust.dist <- c()
  BA.edges.total <- c()
  
  #Generate graphs for nei.1
  for(i in 1:iterations){
    #Generate random sequence of numbers for feeding into the BA model 
    BA.deg.seq <- runif(n, 0, max(igraph::degree(G)))
    #Ensure number of edges to be added in BA model should be approximately the same 
    #as in the input graph. By generating these randomly each time we add another element
    #of randomization to the generation of the BA model.
    BA.deg.seq <- round(e * (BA.deg.seq / sum(BA.deg.seq)))
    BA.deg.seq <- sort(BA.deg.seq) #sort so low numbers added first and edges are unlikely to be omitted
    
    BA <- barabasi.game(n, directed=FALSE, out.seq = BA.deg.seq)
    TR <- append(TR, transitivity(BA, type="global"))
    GE <- append(GE, Global_efficiency(BA, weighted=FALSE))
    
    deg.dist <- c(deg.dist, igraph::degree(BA))
    clust.dist <- c(clust.dist, transitivity(BA, type="local"))
    BA.edges.total <- c(BA.edges.total, ecount(BA))
  }
  
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  df_out <- data.frame("Transitivity" = TR_out, "Global efficiency" = GE_out, row.names=c("Average", "stdev"))
  df.distributions <- data.frame("degree" = deg.dist, "transitivity" = clust.dist, "edges" = BA.edges.total)
  return (list(df_out, df.distributions))
  
}

get.JS <- function(x,y, bandwidth='Silverman'){
  #Calculate the Jensen-Shannon divergence between two adjacency matrices. See 
  # Takahashi et al, 2012, PLoS One
  # Code more or less taken from the statGraph package but modified for use for 
  # just a pair of matrices
  #
  #Input: adjacency matrices (x, y) and the badnwidth to use for the spectral
  # estimation
  #
  #Output: a number (the Jensen-Shannon divergence)
  #
  #
  adjacencyMatrices <- list(x, y)
  f <- statGraph:::nSpectralDensities(adjacencyMatrices, bandwidth = bandwidth)
  densities <- f$densities
  x <- f$x
  y1 <- densities[,1]
  y2 <- densities[,2]
  result <- statGraph:::JS(list(x = x, y = y1), list(x = x, y = y2))
  return(result)
}

boot_function <- function(df, indices, thresh=0.01){
  df_boot <- df[indices,]
  df_thresh <- corr_matrix_threshold(df_boot, thresh=thresh, neg_Rs=FALSE, thresh.param='p')
  G <- graph.adjacency(as.matrix(df_thresh), mode="undirected",
                       weighted=TRUE)
  
  G <- decompose.graph(G)[[1]]
  GC <- GC_size(G)
  GE <- Global_efficiency(G, weighted=FALSE)
  Trans <- transitivity(G, type="global")
  edges <- length(E(G))
  
  return(c(GC, GE, Trans, edges))
}

bootstrap_measures <- function(df, iterations=500, thresh=0.01, conf_interval=0.95){
  #Input: dataframe of counts (noramlized/cleaned), number of iterations
  #for bootstrapping and p-value for tresholding network, and confidence interval
  #
  #Output: Bootstrap object
  
  return(boot(data=df, statistic=boot_function, R=iterations, thresh=thresh))
  
  
}
