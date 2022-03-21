# Ian Rambo
# Created: April 26, 2021

# Purpose: create an interactive network
# from NetAn undirected network analysis outputs

# Usage: change the paths under INPUT FILES and OUTPUT FILES
#=============================================================================
install.packages("visNetwork")

#for development version
#devtools::install_github("datastorm-open/visNetwork")

library(visNetwork)
library(tidyverse)
library(RColorBrewer)
library(tools)
#-----------------------------------------------------------------------------
create.network <- function(n, e, tmain, tsub,
                           select_by, h = "700px", w = "100%") {
  ### Visualize an interactive network using visNetwork ###
  visNetwork(n, e, height = h, width = w,
             main = tmain,
             submain = tsub) %>%
    visEdges(arrows="to") %>%
    visOptions(selectedBy = select_by,
               highlightNearest = TRUE,
               nodesIdSelection = TRUE,
               manipulation = TRUE) %>%
    visPhysics(stabilization = FALSE) %>%
    visConfigure(enabled = TRUE) }

write.network <- function(network, network_path) {
  ### Write an output visNetwork html ###
  ### In-progress: md5sum of file ###

  #Write the network to output
  visSave(network, file = network_output)
  #Get the md5 checksum
  # network_md5_file <- file.path(network_path, ".md5sum")
  #
  # network_md5sum <- as.vector(tools::md5sum(network_output))
  # write.table(as.data.frame(network_md5sum),
  #             network_md5_file,
  #             quote = FALSE, col.names = FALSE,
  #             row.names = FALSE, fileEncoding = "UTF-8",
  #             sep = "\t")
}


#Create a data frame of node shapes and colors
node_shapes <- c("square", "triangle", "box",
                 "dot", "star", "ellipse", "diamond")

#Qualitative palette colors
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]

qual_col_vec <- unlist(mapply(brewer.pal,
                              qual_col_pals$maxcolors, rownames(qual_col_pals)))

node_col_shape <- data.frame(qual_col_vec,
                             rep(node_shapes, length.out = length(qual_col_vec)))

colnames(node_col_shape) <- c("color", "shape")
#=============================================================================
#NETWORK SETUP
#-----------------------------------------------------------------------------
#Set your working directory for input files
setwd("../Data")

# INPUT FILES

# Tab-delimited output file from extract_netan_groups.py
# Connected components - group numbers and members
# CHANGE this file path to your connected components file
components_file <- "GuaymasC_CRISPRDetect_output_qt3_rl23_spacers___GuaymasC_VIBRANT_extracted_phages_all_high_medium_95pid_network_connected-components-modnames.txt"

# Tab-delimited input file to NetAn containing the nodes information
netan_file <- "GuaymasC_CRISPRDetect_output_qt3_rl23_spacers___GuaymasC_VIBRANT_extracted_phages_all_high_medium_95pid_network.txt"

# Optional taxonomy mapping file
# CHANGE this file path to your taxonomy file
taxonomy_file <- "GuaymasC_Archaea_Bacteria_Taxonomy_Master.txt"
#-----------------------------------------------------------------------------
# OUTPUT FILES

#Set your desired output directory
output_directory = "./visnetwork"

# Network output in HTML format that can be viewed in an internet browser
network_output <- file.path(output_directory,
                          paste("undirected_network_visNetwork_", Sys.Date(), ".html", sep = ''))

#-----------------------------------------------------------------------------
#NetAn input file containing node connections
netan <- read.table(netan_file,
                    col.names = c("from", "to"),
                    sep = "\t")

#Data frame of connected components (i.e. nodes)
cc_groups <- read.table(components_file,
                        col.names = c("group", "Bin"),
                        sep = "\t")

#Add taxonomic labels if an input taxonomy table is specified
add_taxonomy <- FALSE
taxonomy_df <- read.table(taxonomy_file,
                          sep = "\t", header = TRUE) %>%
  dplyr::select(Bin, Taxonomy)

#Optional label for nodes not included in taxonomy_df
taxa_na_label <- "Guaymas21 Virus"

#Merge taxonomy information to the network groups if present
if(exists("taxonomy_df") && is.data.frame(get("taxonomy_df"))) {
  print("INFO: Adding taxonomy labels to network groups")
  add_taxonomy = TRUE
  cc_groups = cc_groups %>% dplyr::left_join(taxonomy_df, by = "Bin")
  #Add a taxonomic identifier for unclassified nodes
  cc_groups$Taxonomy[which(is.na(cc_groups$Taxonomy) & grepl("scaffold_[0-9]+", cc_groups$Bin))] = taxa_na_label

} else {
  print("INFO: Taxonomy data frame not found, ignoring")
}

#-----------------------------------------------------------------------------
#Create a data frame of shapes and colors for taxonomic identifiers
node_taxonomy_uniq <- unique(cc_groups$Taxonomy)

node_taxonomy_col_shape <- data.frame(node_taxonomy_uniq,
                                      node_col_shape$color[1:length(node_taxonomy_uniq)],
                                      node_col_shape$shape[1:length(node_taxonomy_uniq)])
colnames(node_taxonomy_col_shape) <- c("Taxonomy", "color", "shape")
#-----------------------------------------------------------------------------
nodes <- cc_groups
edges <- netan

if(add_taxonomy == TRUE) {
  nodes <- nodes %>% dplyr::left_join(node_taxonomy_col_shape) %>%
    dplyr::distinct() %>%
    tidyr::unite("id", c("Bin", "Taxonomy"),
                 sep = "; ", remove = FALSE)

  edges <- edges %>% dplyr::rename(Bin = from) %>%
    dplyr::left_join(nodes %>% select(id, Bin), by = "Bin") %>%
    dplyr::rename(from = id) %>%
    dplyr::select(-Bin) %>%
    dplyr::rename(Bin = to) %>%
    dplyr::left_join(nodes %>% select(id, Bin), by = "Bin") %>%
    dplyr::rename(to = id) %>%
    dplyr::select(-Bin)

}

#Create the network
network <- NULL

if(add_taxonomy == TRUE) {
  network <- create.network(n = nodes, e = edges, tmain = "Undirected Network",
                            tsub = "NetAn v1.0",
                            select_by = "Taxonomy")
} else {
  network <- create.network(n = nodes, e = edges, tmain = "Undirected Network",
                            tsub = "NetAn v1.0")
}

network

#Export the network as html
write.network(network, network_output)
