library("tidytree")
library("dplyr")
library("castor")
library("phytools")
library("treeio")
library("phangorn")

determine_level <- function(row_ind, spp_df, rank_ids_list){
  if(spp_df[row_ind, "species"] %in% rank_ids_list[["species"]]){
    return("species")
  } else if (spp_df[row_ind, "species_group"] %in% rank_ids_list[["species_group"]]){
    return("species_group")
  } else if (spp_df[row_ind, "genus"] %in% rank_ids_list[["genus"]]){
    return("genus")
  } else if (spp_df[row_ind, "family"] %in% rank_ids_list[["family"]]){
    return("family")
  } else if (spp_df[row_ind, "order"] %in% rank_ids_list[["order"]]){
    return("order")
  } else if (spp_df[row_ind, "class"] %in% rank_ids_list[["class"]]){
    return("class")
  } else if (spp_df[row_ind, "phylum"] %in% rank_ids_list[["phylum"]]){
    return("phylum")
  }else if (spp_df[row_ind, "superkingdom"] %in% rank_ids_list[["superkingdom"]]){
    return("superkingdom")
  }
}

find_node_and_level <- function(row_ind, spp_df, rank_ids_list, phylo_file, phylo_file_with_lineage){
  output <- list()
  level <- determine_level(row_ind, spp_df, rank_ids_list)
  taxid <- spp_df[row_ind, level]
  vector_sharing_the_same_taxon <- phylo_file_with_lineage %>% dplyr::filter(.[[level]] == taxid) %>% pull(node)
  node <- ifelse(length(vector_sharing_the_same_taxon) > 1,
                 getMRCA(phylo_file, vector_sharing_the_same_taxon),
                 tidytree::MRCA(phylo_file,vector_sharing_the_same_taxon))
  output[["Node"]] <- node
  output[["Level"]] <- level
  return(output)
}

change_label_to_node <- function(phylo_file, tip_label_vector){
  phylo_file_tbl <- as_tibble(phylo_file)
  ind <- which(phylo_file_tbl$label %in% tip_label_vector)
  phylo_file_tbl$node[ind]
}

compute_len <- function(node_number, phylo_file, tip_label_vector){
  tryCatch({
    output <- list()
    if(isTip(phylo_file, node_number)){
      output[["matrix"]] <- matrix()
      output[["length"]] <- 0
      output[["subtree"]] <- NA
      return(output)
      } else {
        subtree <- get_subtree_with_tips(phylo_file, tip_label_vector)$subtree
        root_number <- find_root(subtree)
        node_vector <- change_label_to_node(subtree, tip_label_vector)
        distance_matrix <- get_all_pairwise_distances(subtree)
        rownames(distance_matrix) <- as.character(1: dim(distance_matrix)[1])
        colnames(distance_matrix) <- as.character(1: dim(distance_matrix)[2])
        rownames(distance_matrix)[node_vector] <- tip_label_vector
        colnames(distance_matrix)[node_vector] <- tip_label_vector
        distance_vector <- c()
        for(i in node_vector){
          distance_vector <- c(distance_vector, distance_matrix[root_number, i])
          }
        output[["matrix"]] <- distance_matrix
        output[["length"]] <- mean(distance_vector)
        output[["subtree"]] <- subtree
        return(output)
        }
    }, error = function(e){
      output <- list()
      output[["matrix"]] <- matrix()
      output[["length"]] <- NA
      output[["subtree"]] <- NA
      return(output)
    })
}

outliers_two_ref_tips <- function(tip_label_vector, phylo_file, phylo_file_with_lineage){
  tryCatch({
    mapping_level <- list()
    rank_score <- list("species" = 1, "species_group" = 2, "genus" = 3, "family" = 4, 
                       "order" = 5, "class" = 6, "phylum" = 7, "superkingdom" = 8)
    mapping_score <- list()
    mapping_distances <- list()
    temp_df <- phylo_file_with_lineage %>% dplyr::filter(label %in% tip_label_vector)
    for (i in 1: length(tip_label_vector)){
      tip_label <- tip_label_vector[i]
      subtree <- tree_subset(phylo_file, tip_label, levels_back = 2)
      temp_labels <- subtree$tip.label
      temp_lineage_table <- phylo_file_with_lineage %>% 
        dplyr::filter(label %in% temp_labels[!(temp_labels %in% tip_label_vector)])
      temp_ranks <- list()
      for (taxo_rank in c("species", "species_group", "genus", "family", "order", "class", "phylum", "superkingdom")){
        unique_ids <- unique(temp_lineage_table[[taxo_rank]])
        temp_ranks[[taxo_rank]] <- unique_ids[!is.na(unique_ids)]
      }
      mapping_level[[tip_label]] <- determine_level(i, temp_df, temp_ranks)
      mapping_score[[tip_label]] <- rank_score[[mapping_level[[tip_label]]]]
      mapping_distances[[tip_label]] <- max(nodeHeights(subtree))
    }
    tip_at_higher_level <- names(mapping_score)[mapping_score == max(unlist(mapping_score))]
    if (length(tip_at_higher_level) == 1){
      abnormal_tip <- tip_at_higher_level
    } else {
      abnormal_tip <- names(mapping_distances)[mapping_distances == max(unlist(mapping_distances))][1]
    }
    return(abnormal_tip)
    }, error = function(e){
      abnormal_tip <- tip_label_vector[1]
      return(abnormal_tip)
      }
  )
}

outliers_causing_NAs <- function(tip_label_vector, reference_tip, count_limit, phylo_file, reference_node){
  consecutive_success <- 0
  output <- c()
  for (tip_label in tip_label_vector){
    if(consecutive_success >= count_limit){
      break
      } else {
        temp_node <- getMRCA(phylo_file,c(tip_label, reference_tip))
        if (temp_node == reference_node){
          output <- c(output, tip_label)
        } else if (is.na(compute_len(temp_node, phylo_file, c(tip_label, reference_tip))[["length"]]) == TRUE){
          output <- c(output, tip_label)
        } else {
          consecutive_success <- consecutive_success + 1
        }
      }
  }
  return(output)
}

detect_outliers_main <- function(branch_len, tip_label_vector, phylo_file, phylo_file_with_lineage, count_limit, 
                                 reference_node, level, distance_output){
  outliers_names <- c("NAs", "species", "species_group", "genus")
  outliers <- sapply(outliers_names, function(x) NULL)
  if(is.na(branch_len) == TRUE){
    if(length(tip_label_vector) == 2){
      outliers[["NAs"]] <- c(outliers[["NAs"]], outliers_two_ref_tips(tip_label_vector, phylo_file, 
                                                                     phylo_file_with_lineage))
    } else {
      middle_ind <- floor((length(tip_label_vector) + 1) / 2)
      middle_tip_label <- tip_label_vector[middle_ind]
      forward_half_vector <- tip_label_vector[1: (middle_ind - 1)]
      reverse_half_vector <- rev(tip_label_vector[(middle_ind + 1): length(tip_label_vector)])
      forward_NAs <- outliers_causing_NAs(forward_half_vector, middle_tip_label, count_limit, phylo_file, 
                                          reference_node)
      reverse_NAs <- outliers_causing_NAs(reverse_half_vector, middle_tip_label, count_limit, phylo_file, 
                                          reference_node)
      outliers[["NAs"]] <- c(outliers[["NAs"]], forward_NAs, reverse_NAs)
    }
  } else if (level %in% c("species", "species_group", "genus")){
    subtree <- distance_output[["subtree"]]
    if (length(tip_label_vector) == 1){
      invisible()
    } else if (length(tip_label_vector) == 2){
      temp_node_vector_original <- change_label_to_node(phylo_file, tip_label_vector)
      if (max(get_all_distances_to_root(subtree)) / 
          max(get_all_distances_to_root(phylo_file)[temp_node_vector_original]) >= n_outlier_two) {
        temp_name <- outliers_two_ref_tips(tip_label_vector, phylo_file, phylo_file_with_lineage)
        if (level == "species"){
          outliers[["species"]] <- c(outliers[["species"]], temp_name)
        } else if (level == "species_group"){
          outliers[["species_group"]] <- c(outliers[["species_group"]], temp_name)
        } else {
          outliers[["genus"]] <- c(outliers[["genus"]], temp_name)
        }
      }
    } else {
      distance_df <- as.data.frame(distance_output[["matrix"]])[tip_label_vector, tip_label_vector]
      distance_df <- distance_df %>% mutate(meanValue = rowSums(distance_df)/(length(tip_label_vector) - 1))
      mean_dist <- mean(distance_df$meanValue)
      sd_dist <- sd(distance_df$meanValue)
      for(rowname in rownames(distance_df)){
        if((level == "species") & (distance_df[rowname, "meanValue"] > (mean_dist + n_sd_species * sd_dist))){
          outliers[["species"]] <- c(outliers[["species"]], rowname)
        } else if ((level == "species_group") & (distance_df[rowname,"meanValue"] > (mean_dist + 
                                                                                     n_sd_species_group * sd_dist))){
          outliers[["species_group"]] <- c(outliers[["species_group"]], rowname)
        } else if ((level == "genus") & (distance_df[rowname,"meanValue"] > (mean_dist + n_sd_genus * sd_dist))){
          outliers[["genus"]] <- c(outliers[["genus"]], rowname)}
      }
    }
  } else {
    invisible()
  }
  return(outliers)
}

insert_tips <- function(spp_df, phylo_file, phylo_file_with_lineage, rank_ids_list, detect_outliers = TRUE){
  final_output <- list()
  summarization_df <- tibble(taxID = numeric(), name = character(), level = character(), node = numeric(), 
                             branch_length = double())
  ave_proc_time <- 0
  start_time <- Sys.time()
  outliers_names <- c("NAs", "species", "species_group", "genus")
  outliers_summary <- sapply(outliers_names, function(x) NULL)
  n <- nrow(spp_df)
  for (row in 1: n){
    cat("\r>>> The ", format(as.character(row), width = 5, justify = "centre"), "th spp. being processed. Progress: ", 
        format(as.character(round(row * 100 / n, 2)), width = 6, justify = "centre", nsmall = 2), "%; ETR: ", 
        format(as.character(round(((n - row)) * ave_proc_time, 2)), width = 5, justify = "centre", nsmall = 2), 
        " min. <<<", sep = "")
    taxid <- spp_df[row, "taxID"]
    node_level_result <- find_node_and_level(row, spp_df, rank_ids_list, phylo_file, phylo_file_with_lineage)
    level <- node_level_result[["Level"]]
    node <- node_level_result[["Node"]]
    temp_vector <- phylo_file_with_lineage %>% dplyr::filter(.[[level]] == spp_df[row, level]) %>% pull(label)
    distance_output <- compute_len(node, phylo_file, temp_vector)
    len <- distance_output[["length"]]
    name <- spp_df[row, "Scientific_Name"]
    record <- tibble(taxID = taxid, name = name, level = level, node = node, branch_length = len)
    summarization_df <- rbind(summarization_df, record)
    if(detect_outliers == TRUE){
      outliers <- detect_outliers_main(len, temp_vector, phylo_file, phylo_file_with_lineage, 5, node, level, 
                                       distance_output)
      outliers_summary[["NAs"]] <- c(outliers_summary[["NAs"]], outliers[["NAs"]])
      outliers_summary[["species"]] <- c(outliers_summary[["species"]], outliers[["species"]])
      outliers_summary[["species_group"]] <- c(outliers_summary[["species_group"]], outliers[["species_group"]])
      outliers_summary[["genus"]] <- c(outliers_summary[["genus"]], outliers[["genus"]])
    }
    current_time <- Sys.time()
    elapsed_time <- current_time - start_time
    ave_proc_time <- as.numeric(elapsed_time, units = "mins") / row
  }
  cat("\n")
  final_output[["Table"]] <- summarization_df
  final_output[["Outliers"]] <- outliers_summary
  return(final_output)
}

check_spp_source <- function(taxid, external_taxids, phylo_tree_taxids){
  Both <- intersect(external_taxids, phylo_tree_taxids)
  external_only <- setdiff(external_taxids, Both)
  phylo_tree_only <- setdiff(phylo_tree_taxids, Both)
  ifelse(taxid %in% Both, "Both", ifelse(taxid %in% external_only, "Added Source Only", "Phylogeny Only"))
}

summarize_phylo <- function(phylo_file, phylo_file_lineage ){
  phylo_with_lineage <- left_join(as_tibble(phylo_file), phylo_file_lineage, by = c("label" = "Tip_Label"))
  tips <- length(phylo_file$tip.label)
  phylo_species <- unique(phylo_with_lineage$species)[!is.na(unique(phylo_with_lineage$species))]
  phylo_species_group <- unique(phylo_with_lineage$species_group)[!is.na(unique(phylo_with_lineage$species_group))]
  phylo_genus <- unique(phylo_with_lineage$genus)[!is.na(unique(phylo_with_lineage$genus))]
  phylo_family <- unique(phylo_with_lineage$family)[!is.na(unique(phylo_with_lineage$family))]
  phylo_order <- unique(phylo_with_lineage$order)[!is.na(unique(phylo_with_lineage$order))]
  phylo_class <- unique(phylo_with_lineage$class)[!is.na(unique(phylo_with_lineage$class))]
  phylo_phylum <- unique(phylo_with_lineage$phylum)[!is.na(unique(phylo_with_lineage$phylum))]
  phylo_superkingdom <- unique(phylo_with_lineage$superkingdom)[!is.na(unique(phylo_with_lineage$superkingdom))]
  summary_df <- data.frame(tips, length(phylo_species), length(phylo_species_group), length(phylo_genus), 
                           length(phylo_family), length(phylo_order), length(phylo_class), length(phylo_phylum), 
                           length(phylo_superkingdom))
  names(summary_df) <- c("Tips", "Species", "SpeciesGroup", "Genus", "Family", "Order", "Class", "Phylum", 
                         "Superkingdom")
  return(summary_df)
}

summarize_insertion <- function(summary_df){
  mapped_level <- c("species", "species_group", "genus", "family", "order", "class", "phylum")
  stats_df <- data.frame(mapped_level)
  count_n <- function(mapped_level, df){
    filtered_df <- df %>% dplyr::filter(level == mapped_level)
    return(nrow(filtered_df))
  }
  mean_branch_len <- function(mapped_level, df){
    filtered_df <- df %>% dplyr::filter(mapped_level == level)
    return(mean(filtered_df$branch_length))
  }
  stats_df <- stats_df %>% mutate(Count_spp = sapply(stats_df$mapped_level, count_n, df = summary_df))
  stats_df <- stats_df %>% mutate(Avg_branch_len = sapply(stats_df$mapped_level, mean_branch_len, df = summary_df))
  vec_br <- summary_df$branch_length
  stats_df <- stats_df %>% add_row(mapped_level = "overall", Count_spp = nrow(summary_df), 
                                   Avg_branch_len = mean(vec_br))
  stats_df
}

args <- commandArgs(trailingOnly = TRUE)
input_phylo_path <- args[1]
output_dir_path <- args[2]
n_sd_species <- as.numeric(args[3])
n_sd_species_group <- as.numeric(args[4])
n_sd_genus <- as.numeric(args[5])
n_outlier_two <- as.numeric(args[6])

if (!dir.exists(output_dir_path)){
  dir.create(output_dir_path)
}

input_tip_lineage_path <- file.path(output_dir_path, "tip_lineage.tsv")
input_phylo_spp_path <- file.path(output_dir_path, "phylo_spp_lineage.tsv")
input_external_spp_path <- file.path(output_dir_path, "added_spp_lineage.tsv") 
output_summary_path <- file.path(output_dir_path, "adding_tip_summary.tsv")
output_name_phylo_path <- file.path(output_dir_path, "expanded_names.tree")
output_id_phylo_path <- file.path(output_dir_path, "expanded_taxIDs.tree")
output_simp_name_phylo_path <- file.path(output_dir_path, "simplified_names.tree")
output_simp_id_phylo_path <- file.path(output_dir_path, "simplified_taxIDs.tree")

phylo_tree <- read.tree(input_phylo_path)
phylo_tree_tbl <- as_tibble(phylo_tree)
phylo_tree_lineage <- read.delim(input_tip_lineage_path)
phylo_tree_spp <- read.delim(input_phylo_spp_path)
external_spp <- read.delim(input_external_spp_path)
phylo_tree_with_lineage <- left_join(phylo_tree_tbl, phylo_tree_lineage, by = c("label" = "Tip_Label"))

cat("Basic statistics about the original phylogeny:\n")
print(summarize_phylo(phylo_tree, phylo_tree_lineage), row.names = FALSE)
cat("\n")

phylo_tree_rank_ids <- list()
for (taxo_rank in c("species", "species_group", "genus", "family", "order", "class", "phylum", "superkingdom")){
  unique_ids <- unique(phylo_tree_with_lineage[[taxo_rank]])
  phylo_tree_rank_ids[[taxo_rank]] <- unique_ids[!is.na(unique_ids)]
}

cat("Checking if any species from both sources can only be mapped at the superkingdom level where the MRCA node is ", 
    "the root node and calculating the branch length is computationally intensive and biologically meaningless:\n", 
    sep = "")

external_superkingdom <- c()
phylo_tree_superkingdom <- c()
for (row in 1: nrow(external_spp)){
  if (determine_level(row, external_spp, phylo_tree_rank_ids) == "superkingdom"){
    external_superkingdom <- append(external_superkingdom, external_spp[row,"taxID"])
  }
}
for (row in 1: nrow(phylo_tree_spp)){
  if (determine_level(row, phylo_tree_spp, phylo_tree_rank_ids) == "superkingdom"){
    phylo_tree_superkingdom <- append(phylo_tree_superkingdom, phylo_tree_spp[row,"taxID"])
  }
}

if((length(external_superkingdom) == 0) & (length(external_superkingdom) == 0)){
  cat("All species can be mapped at taxonomic ranks lower than superkingdom.\n")
} else {
  for (item in external_superkingdom){
    ind <- which(external_spp$taxID == item)
    name <- external_spp[ind, "Scientific_Name"]
    cat("TaxID ", item, " (", name, ") is removed from the user-provided source and will not be added to the ", 
        "phylogeny as it can only be mapped at the superkingdom level.\n", sep = "")}
  for (item in phylo_tree_superkingdom){
    ind <- which(phylo_tree_spp$taxID == item)
    name <- phylo_tree_spp[ind, "Scientific_Name"]
    cat("TaxID ", item, " (", name, ") is removed from the phylogenetic source and will not be added to the ", 
        "phylogeny as it can only be mapped at the superkingdom level.\n", sep = "")}
}
cat("\n")

external_spp <- external_spp %>% dplyr::filter(!taxID %in% external_superkingdom)
phylo_tree_spp <- phylo_tree_spp %>% dplyr::filter(!taxID %in% phylo_tree_superkingdom)

outliers <- list()
summary_table_initial <- tibble(taxID = numeric(), name = character(), level = character(), node = numeric(),
                                branch_length = double())

cat("*** Step 1/5 ***\n")
cat("Initial mapping of species from the user-provided source to identify potential outlier tips.\n")
add_external_initial <- insert_tips(external_spp, phylo_tree, phylo_tree_with_lineage, phylo_tree_rank_ids, TRUE)
summary_table_initial <- rbind(summary_table_initial, add_external_initial$Table)
cat("Identification of potential outlier tips completed!\n\n")

for (rank_name in names(add_external_initial$Outliers)){
  outliers[[rank_name]] <- c(outliers[[rank_name]], add_external_initial$Outliers[[rank_name]])
}

cat("*** Step 2/5 ***\n")
cat("Initial mapping of species from the original phylogeny to identify potential outlier tips.\n")
phylo_tree_spp_left <- phylo_tree_spp %>% 
  dplyr::filter(taxID %in% setdiff(phylo_tree_spp$taxID, external_spp$taxID))
if (nrow(phylo_tree_spp_left) == 0){
  add_phylo_tree_initial <- list("Table" = tibble(taxID = numeric(), name = character(), level = character(), 
                                                  node = numeric(), branch_length = double()),
                                 "Outliers" = list())
} else {
  add_phylo_tree_initial <- insert_tips(phylo_tree_spp_left, phylo_tree, phylo_tree_with_lineage, 
                                      phylo_tree_rank_ids, TRUE)
}
summary_table_initial <- rbind(summary_table_initial, add_phylo_tree_initial$Table)


cat("Identification of potential outlier tips completed!\n\n")

cat("Basic statistics about insertion of species taxIDs before removal of potential outlier tips.\n", sep = "")
print(summarize_insertion(summary_table_initial), row.names = FALSE)
cat("\n")

for (rank_name in names(add_phylo_tree_initial$Outliers)){
  outliers[[rank_name]] <- c(outliers[[rank_name]], add_external_initial$Outliers[[rank_name]])
}

to_delete_tips <- unique(unlist(outliers))

phylo_tree_updated <- ape::drop.tip(phylo_tree, to_delete_tips)
phylo_tree_tbl_updated <- as_tibble(phylo_tree_updated)
phylo_tree_with_lineage_updated <- left_join(phylo_tree_tbl_updated, phylo_tree_lineage, by = c("label" = "Tip_Label"))

phylo_tree_rank_ids_updated <- list()
for (taxo_rank in c("species", "species_group", "genus", "family", "order", "class", "phylum", "superkingdom")){
  unique_ids <- unique(phylo_tree_with_lineage_updated[[taxo_rank]])
  phylo_tree_rank_ids_updated[[taxo_rank]] <- unique_ids[!is.na(unique_ids)]
}

cat("### Basic statistics about the phylogeny after removal of outlier tips ###\n")
print(summarize_phylo(phylo_tree_updated, phylo_tree_lineage), row.names = FALSE)
cat("\n")

cat("Checking if any species from both sources can only be mapped at the superkingdom level after deletion of ", 
    "potential outlier tips from the original phylogeny:\n", sep = "")

external_superkingdom_updated <- c()
phylo_tree_superkingdom_updated <- c()
for (row in 1: nrow(external_spp)){
  if (determine_level(row, external_spp, phylo_tree_rank_ids_updated) == "superkingdom"){
    external_superkingdom_updated <- append(external_superkingdom_updated, external_spp[row,"taxID"])
  }
}
for (row in 1: nrow(phylo_tree_spp)){
  if (determine_level(row, phylo_tree_spp, phylo_tree_rank_ids_updated) == "superkingdom"){
    phylo_tree_superkingdom_updated <- append(phylo_tree_superkingdom_updated, phylo_tree_spp[row,"taxID"])
  }
}

if((length(external_superkingdom_updated) == 0) & (length(external_superkingdom_updated) == 0)){
  cat("All species can be mapped at taxonomic ranks lower than superkingdom.\n")
} else {
  for (item in external_superkingdom_updated){
    ind <- which(external_spp$taxID == item)
    name <- external_spp[ind, "Scientific_Name"]
    cat("TaxID ", item, " (", name, ") is removed from the user-provided source and will not be added to the ", 
        "phylogeny as it can only be mapped at the superkingdom level.\n", sep = "")}
  for (item in phylo_tree_superkingdom_updated){
    ind <- which(phylo_tree_spp$taxID == item)
    name <- phylo_tree_spp[ind, "Scientific_Name"]
    cat("TaxID ", item, " (", name, ") is removed from the phylogenetic source and will not be added to the ", 
        "phylogeny as it can only be mapped at the superkingdom level.\n", sep = "")}
}
cat("\n")

external_spp_updated <- external_spp %>% dplyr::filter(!taxID %in% external_superkingdom_updated)
phylo_tree_spp_updated <- phylo_tree_spp %>% dplyr::filter(!taxID %in% phylo_tree_superkingdom_updated)

summary_table_final <- tibble(taxID = numeric(), name = character(), level = character(), node = numeric(),
                              branch_length = double())

cat("*** Step 3/5 ***\n")
cat("Final mapping of species from the user-provided source after removal of outlier tips.\n")
add_external_final <- insert_tips(external_spp_updated, phylo_tree_updated, phylo_tree_with_lineage_updated, 
                                  phylo_tree_rank_ids_updated, FALSE)
cat("Done!\n\n")
summary_table_final <- rbind(summary_table_final, add_external_final$Table)

cat("*** Step 4/5 ***\n")
cat("Final mapping of species from the original phylogeny after removal of outlier tips.\n")
phylo_tree_spp_updated_left <- phylo_tree_spp_updated %>% 
  dplyr::filter(taxID %in% setdiff(phylo_tree_spp_updated$taxID, external_spp_updated$taxID))
if (nrow(phylo_tree_spp_updated_left) == 0){
  add_phylo_tree_final <- list("Table" = tibble(taxID = numeric(), name = character(), level = character(), 
                                                node = numeric(), branch_length = double()),
                               "Outliers" = list())
} else{
  add_phylo_tree_final <- insert_tips(phylo_tree_spp_updated_left, phylo_tree_updated, phylo_tree_with_lineage_updated, 
                                      phylo_tree_rank_ids_updated, FALSE)
}
summary_table_final <- rbind(summary_table_final, add_phylo_tree_final$Table)

cat("Done!\n\n")

cat("Basic statistics about insertion of species taxIDs after removal of potential outlier tips.\n", sep = "")
print(summarize_insertion(summary_table_final), row.names = FALSE)
cat("\n")

cat("*** Step 5/5 ***\n")
cat("Generating output files.\n")
output_names_phylo <- add.tips(phylo_tree_updated, tips = summary_table_final$name, where = summary_table_final$node,
                               edge.length = summary_table_final$branch_length)
output_names_phylo <- ape::drop.tip(output_names_phylo, phylo_tree_updated$tip.label)

output_taxid_phylo <- output_names_phylo
for (ind in 1: length(output_taxid_phylo$tip.label)){
  label <- output_taxid_phylo$tip.label[ind]
  summary_table_ind <- which(summary_table_final$name == label)
  taxid <- summary_table_final$taxID[summary_table_ind]
  output_taxid_phylo$tip.label[ind] <- taxid
}

output_summary_table <- summary_table_final %>% select(-node) %>% 
  mutate(source = sapply(summary_table_final$taxID, check_spp_source, external_taxids = external_spp_updated$taxID, 
                         phylo_tree_taxids = phylo_tree_with_lineage_updated$species))

output_simp_id_phylo <- ape::keep.tip(output_taxid_phylo, as.character(external_spp_updated$taxID))
output_simp_name_phylo <- output_simp_id_phylo
for (ind in 1: length(output_simp_name_phylo$tip.label)){
  label <- output_simp_name_phylo$tip.label[ind]
  summary_table_ind <- which(summary_table_final$taxID == label)
  sciname <- summary_table_final$name[summary_table_ind]
  output_simp_name_phylo$tip.label[ind] <- sciname
}

write.table(output_summary_table, file = output_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.tree(output_names_phylo, file = output_name_phylo_path)
write.tree(output_taxid_phylo, file = output_id_phylo_path)
write.tree(output_simp_name_phylo, file = output_simp_name_phylo_path)
write.tree(output_simp_id_phylo, file = output_simp_id_phylo_path)
cat("Done!\n")
