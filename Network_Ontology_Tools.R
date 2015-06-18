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
  #ids <- unlist(lapply(ids, as.numeric)) #Turn ids to numbers in a vector

  parents <- unlist(xmlApply(parents, xmlValue))
  #parents <- unlist(lapply(parents, as.numeric)) #Turn parent ids to numbers in vector
  
  names <- unlist(xmlApply(names, xmlValue))
  
  acronyms <- unlist(xmlApply(acronyms, xmlValue))
  
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

