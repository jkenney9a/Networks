#
#
#


library(XML)
library(igraph)


load.xml <- function(filename){
  # Load up xml document as a parsed C-type structure
  X <- xmlParse(filename)
  return(X)
}


xml.to.igraph <- function(xml.file){
  #
  #Turn a structured xml Allen Brain Atlas ontology file into an igraph graph
  #
  
  xml <- load.xml(xml.file)
  
  #pull out ids, parent ids, names, and acronyms in order from xml file
  #note that this should return each list in the same order.
  ids <- getNodeSet(xml, "//id")
  parents <- getNodeSet(xml, "//parent-structure-id")
  names <- getNodeSet(xml, "//name")
  acronyms <- getNodeSet(xml, "//acronym")
  
  #Tidy up the xml structures so can be put in a dataframe
  ids <- unlist(xmlApply(ids, xmlValue))

  parents <- unlist(xmlApply(parents, xmlValue))
  
  names <- unlist(xmlApply(names, xmlValue)) #Turn into a vector
  names <- gsub('\"', '', names) #Remove some XML formatting
  
  acronyms <- unlist(xmlApply(acronyms, xmlValue))
  acronyms <- gsub('\"', '', acronyms) #Remove some XML formatting
  
  #Map ids to acronyms
  id.acronyms <- acronyms
  names(id.acronyms) <- ids
  
  #Map ids to full names
  id.names <- names
  names(id.names) <- ids
  
  parents[[1]] <- ids[[1]] #Takes care of null/NA in the root
  
  df.edge.list <- data.frame(parent=acronyms, child=id.acronyms[parents])
  df.attributes <- data.frame(parent=acronyms, full.name=names, node.id=ids)
  
  G <- graph.data.frame(df.edge.list, vertices=df.attributes)
  
  return(G)  
}


get.coarse.ontology <- function(ontology.graph, nodes, coarseness){
  #
  # Input: Ontology graph, a vector of nodes (that should be in the ontology graph)
  # and a vector of the level of coarseness in which to label each node (e.g, 
  # "hippocampus", "cortex" etc.). The coarseness level must correspond to a part
  # of the ontology graph and right now must be the node acronym
  #
  # Output: A dataframe consiting of the nodes and their associated ontology level
  # from the coarseness vector
  
  node.mapping <- data.frame("node"=nodes, "ont.mapping"=NA)
  ont.diameter <- diameter(ontology.graph)
  
  for (C.level in coarseness){
    
    #Get neighborhood of nodes for coarseness level that go into that node
    hood <- graph.neighborhood(ontology.graph, order=ont.diameter, nodes=C.level, mode='in')
    hood.acronym <- V(hood[[1]])$name #Get a list of node names in 'hood
    hood.full.name <- V(hood[[1]])$full.name #Get list of full names in 'hood
    hood <- c(hood.acronym, hood.full.name) #Can use either node acronyms or full names
    
    #Compare nodes in input list to neighborhood and alter ont.mapping accordingly
    node.mapping$ont.mapping[is.element(nodes,hood)] <- C.level
        
  }
  
  return(node.mapping)
}


subset.and.neighbors <- function(ontology.graph, nodes, coarse.only = TRUE){
  # Input: Ontology graph and a vector of nodes to subset 
  # the ontology graph (as a vector). Can be full node names or acronyms.
  #
  # Output: Only the parts of the ontology graph that are connected to the input node
  # names or acronyms.
  
  #Convert node full names/acronyms to node ids
  #This allows function to be more general on what is used for node input
  
  O.G <- ontology.graph #Little easier to work with, even if more gangsta'
  
  #Consider figuring out how to make this more general in the future (JK)
  node.out <- subset(V(O.G), V(O.G)$full.name %in% nodes | V(O.G)$name %in% nodes)
  
  ont.diameter <- diameter(O.G)
  
  #Get node ids for nodes in neighborhood
  hood <- neighborhood(O.G, order=ont.diameter, nodes=node.out, mode='out')
  if(coarse.only == FALSE){
    hood <- append(hood, neighborhood(O.G, order=ont.diameter, nodes=node.out, mode='in'))
  }
  hood <- unique(unlist(hood)) #get set of all node ids to include in output the graph
  
  O.G.out <- induced.subgraph(O.G, vids=hood, impl="copy_and_delete")
        
  return(O.G.out)
}




