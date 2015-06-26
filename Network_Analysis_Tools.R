#
# A module of tools to assist in analyzing networks using igraph.
# 
# Justin W. Kenney
# jkenney9a@gmail.com
# 

library(igraph) #For network functions
library(brainwaver) #For global efficiency calculation
library(Hmisc) #For generating correlation/p-value matrices
library(boot) #For bootstrapping


load_data <- function(csv_file){
  #   Input: CSV file name
  #   
  #   Output: Dataframe containing CSV file
  
  df <- read.csv(file=csv_file, header=TRUE)
  return(df)
}

clean_data <- function(df_count, missing_thresh=4, fill_missing=TRUE){
  # Input: dataframe of counts to be cleaned, the maximum threshold
  # of missing values allowed otherwise data is removed
  #
  # Output: dataframe of counts with columns with more zeros than 
  # zero_threshold removed. Other zeros are replaced with column mean.
  
  #Replace zeros in data with NA to allow easier handling
  df_count[mapply("==", df_count, 0)] <- NA
  
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
  #   Output: Two dataframes corresponding to 1) all pairwise Pearson 
  #   correlations and 2) all associated un-adjusted p-values


  corr <- rcorr(as.matrix(df), type='pearson')
  df_corr <- as.data.frame(corr['r'])
  df_pvalue <- as.data.frame(corr['P'])
  
  names(df_corr) <- gsub("r.","",colnames(df_corr), fixed=TRUE) #Remove ".r" from node names
  
  #Change ... to -; for some reason when R imports "-" it turns it into "..."
  names(df_corr) <- gsub("...","-",colnames(df_corr), fixed=TRUE)
  rownames(df_corr) <- gsub("...", "-", rownames(df_corr), fixed=TRUE)
  
  return(list("corr" = df_corr,"pvalue" = df_pvalue))
}

corr_matrix_threshold <- function(df, neg_Rs = FALSE, p_threshold=0.01){
  #   Input: Dataframe with headers as titles (brain regions) of counts etc.
  #           p-value threshold, whether or not to keep negative correlations.
  #   
  #   Output: Dataframe of correlations thresholded at p < alpha
  #   
  #   NOTE: Removes diagonals and negative correlations
  
  dfs <- corr_matrix(df)
  df_corr <- dfs[['corr']]
  df_pvalue <- dfs[['pvalue']]
  
  #apply p-value threshold to correlation matrix
  df_corr[mapply(">=", df_pvalue, p_threshold)] <- 0
  
  #remove diagonals (may not be necessary when using igraph...)
  df_corr[mapply("==", df_corr, 1)] <- 0
  
  #remove any NaNs, infs or NAs (sometimes happens with bootstrapping; not sure why)
  df_corr[mapply(is.infinite, df_corr)] <- 0
  
  #remove negative correlations
  if (neg_Rs == FALSE){
    df_corr[mapply("<", df_corr, 0)] <- 0
  }
  
  
  return(df_corr)
}

CSV_to_igraph <- function(CSV_file, negs = FALSE, p_thresh=0.01){
  #Input: CSV_file of counts etc., p-value threshold
  #
  #Output: igraph graph object
  
  df <- load_data(CSV_file)
  df_G <- corr_matrix_threshold(df, neg_Rs = negs, p_threshold=p_thresh)
  
  G <- graph.adjacency(as.matrix(df_G), mode="undirected",
                       weighted=TRUE)
  
  return(G)
}

df_to_igraph <- function(df, negs=FALSE, p_thresh=0.01){
  #Input: df of counts per brain region, whether to keep negative
  #correlations and the p-value threshold
  #
  #Output: igraph graph
  
  df_G <- corr_matrix_threshold(df, neg_Rs = negs, p_threshold=p_thresh)
  G <- graph.adjacency(as.matrix(df_G), mode="undirected", 
                       weighted=TRUE)
  
  return(G)
}

get_centrality_measures <-function(G){
  # Input: An igraph graph
  #
  # Output: A dataframe of centrality measures:
  # Degree, betweenness, eigenvector, closeness, 
  # transitivity (~clustering coefficient), nodal efficiency
  # NOTE: nodal efficiency makes use of weights if they exist
  
  degree <- degree(G)
  between <- betweenness(G, normalized=TRUE, weights=NULL)
  eigenvector <- evcent(G)
  close <- closeness(G)
  trans <- transitivity(G,type="local")
  
  #Need to pull out matrices for efficiency calculations
  adj_mat <- as.matrix(get.adjacency(G))
  weight_mat <- as.matrix(get.adjacency(G, attr='weight'))
  efficiency <- global.efficiency(adj_mat, weight_mat)
  
  output <- data.frame("degree" = degree, "betweenness" = between, 
                       "eigenvector" = eigenvector[[1]], "closeness" = close,
                      "efficiency" = efficiency[[1]], "transitivity" = trans)
  
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
  
  adj_mat <- as.matrix(get.adjacency(G))
  
  if(weighted==TRUE){
    weight_mat <- as.matrix(get.adjacency(G, attr='weight'))
    efficiency <- global.efficiency(adj_mat, weight_mat)
  } else{
    efficiency <- global.efficiency(adj_mat, adj_mat)
  }
    
  return(mean(efficiency[[1]]))
}

Rand_graph_stats <- function(G, iterations = 100, degree_dist=TRUE){
  #Input: igraph graph, wehther or not random graphs should have
  #same degree distribution as input graph. If degree_dist=FALSE
  #uses Erdos-Renyi random graph
  #
  #Output: Equivalent random graph statistics and stdev in a df:
  #Transitivity (very closely related to clustering coefficient)
  #Global efficiency
  
  TR <- c() #Transitivity
  GE <- c() #global efficiency (GE)
  
  if(degree_dist == TRUE){
    #Get rid of degree zero nodes b/c this causes issues with the
    #generation of random graphs using the "vl" algorithm
    G_1 <- delete.vertices(G, which(degree(G) == 0))
    degrees <- degree(G_1)
    
    #Get number of zero degree nodes
    zero_d <- length(degree(G)) - length(degrees)
    
    for(i in 1:iterations){
      RG <- degree.sequence.game(degree(G_1), method='vl')
      TR <- append(TR, transitivity(RG, type="global"))
      #Add back in zero degree nodes to get efficiency of whole graph
      GE <- append(GE, Global_efficiency(RG + zero_d, weighted=FALSE))
      } 
    } else{
  
  if(degree_dist == FALSE){
    n <- vcount(G) #number of nodes
    m <- ecount(G) #number of edges
    for(i in 1:iterations){
      RG <- erdos.renyi.game(n, m, type="gnm")
      TR <- append(TR, transitivity(RG, type="global"))
      GE <- append(GE, Global_efficiency(RG, weighted=FALSE))
      }
    }
  }
      
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  
  return(data.frame("Transitivity" = TR_out, 
                    "Global efficiency" = GE_out,
                    row.names = c("Average", "stdev")))
}

Watts_Strogatz_model <- function(G, iterations=100, trans_match_iter=100){
  #Input: A graph to use as basis for generating Watts-Strogatz small 
  # world model. Number of random Watts-Strogatz graphs to calculate and
  #the number of iterations to use for matching transitivity/clustering
  #coefficient
  #
  #Output: List containing three items:
  # 1) Random graph statistics and stdev in a df for transitivity and global efficiency
  # 2) The average sorted degree distribution for all the WS graphs generated
  # 3) The average sorted clustering distribution for all the WS graphs generated
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
  trans.out <- as.data.frame(matrix(NA, ncol=length(wiring.ps), nrow=iterations, 
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
    
    deg.dist[i,] <- sort(degree(W))
    clust.dist[i,] <- sort(transitivity(W, type="local"))
  }
  
  deg.dist <- colMeans(deg.dist)
  clust.dist <- colMeans(clust.dist)
  
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  df_out <- data.frame("Transitivity" = TR_out, "Global efficiency" = GE_out, row.names=c("Average", "stdev"))

  return (list(df_out, deg.dist, clust.dist))
  
}

BA_model <- function(G, iterations){
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
  
  n <- vcount(G)
  e <- ecount(G)
  G.density <- e / (((n*(n-1))/2))
  nei <- round(e/n) #get approximate number for edges to add per iteration of BA algorithm (m)
  
  nei.out <- 0
  BA.test <- e
  BA.G.diff <- c()
  BA.edges <- c()
  
  # Select nei/m input for BA model that yields closest match to number of edges of input graph
  for(i in 1:(nei+2)){
    BA.temp <- barabasi.game(n, m=i, directed=FALSE)
    diff <- abs(e - ecount(BA.temp))
    BA.edges <- append(BA.edges, ecount(BA.temp))
    BA.G.diff <- append(BA.G.diff, diff) #track differences
  }
  
  names(BA.G.diff) <- seq(1, length(BA.G.diff), 1) #give each diff a name corresponding to the nei
  BA.G.diff <- sort(BA.G.diff, decreasing=FALSE) #Sort the diffs in ascending order
                                                 #first 2 diffs and their name correspond to lowest neis.
  names(BA.edges) <- seq(1, length(BA.G.diff), 1)
  #Sort BA edges in same order as the differences
  BA.edges <- BA.edges[(match(names(BA.G.diff), names(BA.edges)))] 

  nei.1 <- strtoi(names(BA.G.diff)[1])
  nei.2 <- strtoi(names(BA.G.diff)[2])
  
  BA.edge.1 <- BA.edges[1]
  BA.edge.2 <- BA.edges[2]
  
  #Caclulate what fraction of iterations should use nei.1 and 2
  BA.frac.1 <- abs((e - BA.edge.2) / (BA.edge.1 - BA.edge.2))
  BA.frac.2 <- 1 - BA.frac.1
  
  #Choose ceiling of integers so err towards more rather than less iterations
  iter.1 <- ceiling(BA.frac.1 * iterations)
  iter.2 <- ceiling(BA.frac.2 * iterations)
   
  TR <- c()     #Vector to hold transitivity calculations
  GE <- c()     #Vector to hold global efficiency calculations
  deg.dist <- c()
  clust.dist <- c()
  BA.edges.total <- c()
  
  #Generate graphs for nei.1
  for(i in 1:iter.1){
    BA <- barabasi.game(n, m=nei.1, directed=FALSE)
    TR <- append(TR, transitivity(BA, type="global"))
    GE <- append(GE, Global_efficiency(BA, weighted=FALSE))
    
    deg.dist <- c(deg.dist, degree(BA))
    clust.dist <- c(clust.dist, transitivity(BA, type="local"))
    BA.edges.total <- c(BA.edges.total, ecount(BA))
  }
  
  #Generate graphs for nei.2
  for(i in iter.1:(iter.1 + iter.2)){
    BA <- barabasi.game(n, m=nei.2, directed=FALSE)
    TR <- append(TR, transitivity(BA, type="global"))
    GE <- append(GE, Global_efficiency(BA, weighted=FALSE))
    
    deg.dist <- c(deg.dist, degree(BA))
    clust.dist <- c(clust.dist, transitivity(BA, type="local"))
    BA.edges.total <- c(BA.edges.total, ecount(BA))
  }
  
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  df_out <- data.frame("Transitivity" = TR_out, "Global efficiency" = GE_out, row.names=c("Average", "stdev"))
  df.distributions <- data.frame("Degrees" = deg.dist, "Clustering" = clust.dist, "Edges" = BA.edges.total)
  return (list(df_out, df.distributions))
  
}

boot_function <- function(df, indices, p_thresh=0.01){
  df_boot <- df[indices,]
  df_thresh <- corr_matrix_threshold(df_boot, p_threshold=p_thresh, neg_Rs=FALSE)
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
  
  return(boot(data=df, statistic=boot_function, R=iterations, p_thresh=thresh))
  
  
}


  
  
  
  


#   if(plots==TRUE){
#     #Degree distribution plots:
#     h.G.deg <- hist(degree(G), breaks=seq(from=min(degree(G))-0.5, to=max(degree(G))+0.5, by=1), plot=FALSE)
#     h.G.deg$counts <- h.G.deg$counts/sum(h.G.deg$counts) #Recalculate as probability distribution
#     
#     h.WS.deg <- hist(deg.dist, breaks=seq(from=floor(min(deg.dist))-0.5, to=ceiling(max(deg.dist))+0.5, by=1), plot=FALSE)
#     mx.h.WS.deg <- max(h.WS.deg$counts / sum(h.WS.deg$counts))
#                       
#     plot(h.G.deg, ylim=c(0,max(c(max(h.G.deg$counts, mx.h.WS.deg)))))
#     d.deg.dist <- density(deg.dist, na.rm=TRUE, from=min(degree(G), to=max(degree(G))))
#     d.deg.dist$y <- d.deg.dist$y / mx.h.WS.deg
#     lines(d.deg.dist)
#     
#     #Clustering distribution plots:
#     G.trans <- transitivity(G, type="local")
#     h.G.clust <- hist(G.trans, breaks=seq(from=0, to=1, by=0.05), plot=FALSE)
#     h.G.clust$counts <- h.G.clust$counts / sum(h.G.clust$counts)
# 
#     h.WS.clust <- hist(clust.dist, breaks=seq(from=0, to=1, by=0.05), plot=FALSE)
#     mx.h.WS.clust <- max(h.WS.clust$counts / sum(h.WS.clust$counts))
# 
#     plot(h.G.clust, ylim=c(0,max(c(max(h.G.clust$counts, mx.h.WS.clust)))))
#     d.clust.dist <- density(clust.dist, na.rm=TRUE, from=0, to=1)
#     d.clust.dist$y <- d.clust.dist$y / max(h.WS.clust$counts)
#     lines(d.clust.dist)
# 
#   }
