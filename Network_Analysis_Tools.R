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
  names(df_G) <- gsub("r.","",colnames(df_G), fixed=TRUE) #Remove ".r" from node names
  
  #Change ... to -; for some reason when R imports "-" it turns it into "..."
  names(df_G) <- gsub("...","-",colnames(df_G), fixed=TRUE) 
  
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
  #Output: Random graph statistics and stdev in a df:
  #Transitivity (i.e., clustering coefficient)
  #Global efficiency
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
  
  #Generate a series of Watts-Strogatz graphs
  for(i in 1:trans_match_iter){
    trans.out[i,] <- sapply(wiring.ps, function(x) 
                            transitivity(watts.strogatz.game(1,n,nei,x), type="global"))
  }
  
  trans.out.means <- colMeans(trans.out)
  p.value <- wiring.ps[which(abs(trans.out.means-G.trans) == min(abs(trans.out.means-G.trans)))]
  
  TR <- c()
  GE <- c()
  
  for(i in 1:iterations){
    W <- watts.strogatz.game(1,n,nei,p.value)
    TR <- append(TR, transitivity(W, type="global"))
    GE <- append(GE, Global_efficiency(W, weighted=FALSE))
  }
  
  TR_out <- c(mean(TR), sd(TR))
  GE_out <- c(mean(GE), sd(GE))
  
  return (data.frame("Transitivity" = TR_out, "Global efficiency" = GE_out, 
                    row.names=c("Average", "stdev")))
  
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

