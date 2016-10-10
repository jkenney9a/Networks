#
# A module of tools to assist in analyzing networks using igraph.
# 
# Justin W. Kenney
# jkenney9a@gmail.com
# 

library(igraph) #For network functions
library(Hmisc) #For generating correlation/p-value matrices
library(boot) #For bootstrapping


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


corr_matrix <- function(df){
  #   Input: Dataframe with headers as titles (brain regions)
  #   
  #   Output: List of two dataframes corresponding to 1) all pairwise Pearson 
  #   correlations and 2) all associated un-adjusted p-values


  corr <- rcorr(as.matrix(df), type='pearson')
  df_corr <- as.data.frame(corr['r'])
  df_pvalue <- as.data.frame(corr['P'])
  
  names(df_corr) <- gsub("^r.","",colnames(df_corr), fixed=FALSE) #Remove "r." from node names
  
  #Change ... to -; for some reason when R imports "-" it turns it into "..."
  names(df_corr) <- gsub("...","-",colnames(df_corr), fixed=TRUE)
  rownames(df_corr) <- gsub("...", "-", rownames(df_corr), fixed=TRUE)
  names(df_corr) <- gsub("sums.", "", colnames(df_corr), fixed=TRUE)
  rownames(df_corr) <- gsub("sums.", "", rownames(df_corr), fixed=TRUE)
  
  names(df_corr) <- rownames(df_corr)
  
  return(list("corr" = df_corr,"pvalue" = df_pvalue))
}

corr_matrix_threshold <- function(df, neg_Rs = FALSE, thresh=0.01, thresh.param='p'){
  #   Input: Dataframe with headers as titles (brain regions) of counts etc.
  #           p-value threshold, whether or not to keep negative correlations.
  #   
  #   Output: Dataframe of correlations thresholded at p < alpha
  #   
  #   NOTE: Removes diagonals and negative correlations
  
  dfs <- corr_matrix(df)
  df_corr <- dfs[['corr']]
  df_pvalue <- dfs[['pvalue']]
  
  #remove diagonals (may not be necessary when using igraph...)
  diag(df_corr) <- 0
  diag(df_pvalue) <- NA
  #df_corr[mapply("==", df_corr, 1)] <- 0
  
  #remove any NaNs, infs or NAs (sometimes happens with bootstrapping; not sure why)
  df_pvalue[mapply(is.infinite, df_corr)] <- 1
  df_pvalue[mapply(is.nan, df_corr)] <- 1  
  df_corr[mapply(is.infinite, df_corr)] <- 0
  df_corr[mapply(is.nan, df_corr)] <- 0

  
  if(tolower(thresh.param)=='p'){
    #apply p-value threshold to correlation matrix
    df_corr[mapply(">=", df_pvalue, thresh)] <- 0    
  } else if(tolower(thresh.param)=='r'){
    df_corr[mapply("<=", abs(df_corr), thresh)] <- 0
  } else{
    print("Invalid thresholding paramter")
  }
 
  
  
  
  #remove negative correlations
  if (neg_Rs == FALSE){
    df_corr[mapply("<", df_corr, 0)] <- 0
  }
  
  
  return(df_corr)
}

CSV_to_igraph <- function(CSV_file, negs = FALSE, thresh=0.01, thresh.param='p'){
  #Input: CSV_file of counts etc., threshold, whether threshold based on p-value ('p')
  # or correlation(r) value ('r')
  #
  #Output: igraph graph object
  
  df <- load_data(CSV_file)
  df_G <- corr_matrix_threshold(df, neg_Rs = negs, thresh=thresh, thresh.param=thresh.param)
  
  G <- graph.adjacency(as.matrix(df_G), mode="undirected",
                       weighted=TRUE)
  
  return(G)
}

df_to_igraph <- function(df, negs=FALSE, thresh=0.01, thresh.param='p'){
  #Input: df of counts per brain region, whether to keep negative
  #correlations and the p-value threshold
  #
  #Output: igraph graph
  
  df_G <- corr_matrix_threshold(df, neg_Rs = negs, thresh=thresh, thresh.param=thresh.param)
  G <- graph.adjacency(as.matrix(df_G), mode="undirected", 
                       weighted=TRUE)
  
  return(G)
}

get_centrality_measures <-function(G, weighted=FALSE, normalized=FALSE){
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
  
  trans <- transitivity(G,type="local")
  
  #Need to pull out matrices for efficiency calculations
  #adj_mat <- as.matrix(get.adjacency(G))
  #weight_mat <- as.matrix(get.adjacency(G, attr='weight'))
  #efficiency <- Global_efficiency(G)
  
  output <- data.frame("degree" = degree, "betweenness" = between, 
                       "eigenvector" = eigenvector[[1]], "closeness" = close,
                      "transitivity" = trans)
  
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
    eff[!is.finite(eff)] <- 0
    gl.eff <- mean(eff[upper.tri(eff)])
  }else{
    eff <- 1/shortest.paths(G, weights=NA)
    eff[!is.finite(eff)] <- 0
    gl.eff <- mean(eff[upper.tri(eff)])
  }
  
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

Rand_graph_stats <- function(G, iterations = 100, degree_dist=TRUE, weighted=FALSE){
  #Input: igraph graph, wehther or not random graphs should have
  #same degree distribution as input graph. If degree_dist=FALSE
  #uses Erdos-Renyi random graph
  #
  #Output: Equivalent random graph statistics and stdev in a df:
  #Transitivity (very closely related to clustering coefficient)
  #Global efficiency
  
  TR <- c() #Transitivity
  GE <- c() #global efficiency (GE)
  path_length <- c() #Path length
  
  if(degree_dist == TRUE){
    #Get rid of degree zero nodes b/c this causes issues with the
    #generation of random graphs using the "vl" algorithm
    G_1 <- delete.vertices(G, which(igraph::degree(G) == 0))
    degrees <- igraph::degree(G_1)
    
    #Get number of zero degree nodes
    zero_d <- length(igraph::degree(G)) - length(degrees)
    
    for(i in 1:iterations){
      RG <- degree.sequence.game(igraph::degree(G_1), method='vl')
      TR <- append(TR, transitivity(RG, type="global"))
      #Add back in zero degree nodes to get efficiency of whole graph
      GE <- append(GE, Global_efficiency(RG + zero_d, weighted=FALSE))
      path_length <- append(path_length, average.path.length(RG))
      } 
    } else{
  
  if(degree_dist == FALSE){
    n <- vcount(G) #number of nodes
    m <- ecount(G) #number of edges
    for(i in 1:iterations){
      RG <- erdos.renyi.game(n, m, type="gnm")
      TR <- append(TR, transitivity(RG, type="global"))
      GE <- append(GE, Global_efficiency(RG, weighted=FALSE))
      path_length <- append(path_length, average.path.length(RG))
      }
    }
  }
      
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  path_length_out <- c(mean(path_length), sd(path_length))
  
  return(data.frame("Transitivity" = TR_out, 
                    "Global.efficiency" = GE_out,
                    "Path.length" = path_length_out,
                    row.names = c("Average", "stdev")))
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
    clust.dist[i,] <- sort(transitivity(W, type="local"))
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
