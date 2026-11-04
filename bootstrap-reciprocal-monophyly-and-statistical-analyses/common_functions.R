root_tree <- function(input_tree, drop_outgroup = NULL) {
  
  # remove non-rodent sequences
  
  remove_leptomys <- input_tree$tip.label[grep("EC5_Leptomys_paulus_M23838.", input_tree$tip.label)] # actually a bird
  remove_planigale <- input_tree$tip.label[grep("EC6_Pseudomys_delicatulus_M56791_WA_WarrambooCreek.", input_tree$tip.label)]
  
  remove <- c(remove_leptomys, remove_planigale)
  trimmed_tree <- drop.tip(input_tree, input_tree$tip.label[match(remove, input_tree$tip.label)])
  
  if (!is.null(drop_outgroup)) {
    
    # exclude additional tips
    
    remove <- c()
    
    for (t in 1:length(drop_outgroup)) {
      
      # find tips to remove
      
      remove_tip <- trimmed_tree$tip.label[grep(drop_outgroup[t], trimmed_tree$tip.label)]
      
      remove <- c(remove, remove_tip)
      
    }
    
    # drop tips
    trimmed_tree <- drop.tip(trimmed_tree, trimmed_tree$tip.label[match(remove, trimmed_tree$tip.label)])
    
  }
  
  # root on outgroup
  
  root_gerbilliscus <- trimmed_tree$tip.label[grep("Gerbilliscus", trimmed_tree$tip.label)]
  root_deomys <- trimmed_tree$tip.label[grep("Deomys", trimmed_tree$tip.label)]
  
  root_labels <- c(root_gerbilliscus, root_deomys)
  
  rooted_tree <- RootTree(trimmed_tree, root_labels)
  
  return(rooted_tree)
  
}

rename_tips <- function(input_tree, label_old, label_new) {
  
  label_list <- label_new[match(input_tree$tip.label, label_old)]
  input_tree$tip.label <- label_list
  
  return(input_tree)
  
}

rename_tips_sample_ID <- function(input_tree, gene_name = NULL, add_names_old = NULL, add_names_new = NULL) {
  
  # create lookup
  rename_lookup <- tibble(old = input_tree$tip.label)
  
  if (!is.null(gene_name)) {
  
    # create delimiter for removing gene information
    gene_delimiter <- paste0("_", gene_name)
  
    # use gene delimiter to cleave off gene information
    rename_lookup <- rename_lookup %>%
    mutate(new2 = str_split_fixed(old, gene_delimiter, 2)[,1])
  
    # by reversing string of current tip names can cleave off .method_bm.contigNumber from sample_ID without cleaving sample_IDs containing periods
    rename_lookup <- rename_lookup %>%
      mutate(new3 = strReverse(new2)) %>%
      mutate(new4 = str_split_fixed(new3, "[.]", 3)[,3]) %>%
      mutate(new = strReverse(new4)) %>%
      select(old, new)
    
  }
  
  if (is.null(gene_name)) {
   
    rename_lookup <- rename_lookup %>%
      mutate(new3 = strReverse(old)) %>%
      mutate(new4 = str_split_fixed(new3, "[.]", 3)[,3]) %>%
      mutate(new5 = strReverse(new4)) %>%
      select(old, new5)
    
    # remove _R_ from start of tip names
    rename_lookup <- rename_lookup %>%
      mutate(new2 = str_split_fixed(new5, "_R_", 2)[,2]) %>%
      mutate(new = case_when(!new2 == "" ~ new2, 
                             .default = new5)) %>%
      select(old, new)
     
  }
  
  if (!is.null(add_names_old)) {
    
    for (t in 1:length(add_names_old)) {
      
      # create new name for non_sample_ID sequences
      rename_lookup <- rename_lookup %>% mutate(new = case_when(grepl(add_names_old[t], rename_lookup$old) ~ add_names_new[t],
                                                                 .default = new))
      
    }
  }
  
  # restore mito tip labels to original sample names
  renamed_tree <- rename_tips(input_tree, rename_lookup$old, rename_lookup$new)
  
  return(renamed_tree)
  
}

#https://stackoverflow.com/questions/13612967/how-to-reverse-a-string-in-r
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")


average_bootstrap <- function(loci_list, trees_folder_path, trees_prefix, trees_suffix) {
  
  # calculate average bootstrap value
  
  # create empty tibble to store values
  bootstrap_values_tibble <- tibble(loci = c(""), bootstrap = c(""))
  
  # loop through loci trees
  for (n in 1:length(loci_list)) {
    
    gene <- loci_list[n]
    
    # read in loci tree file
    gene_tree <- read.tree(paste0(trees_folder_path, trees_prefix, gene, trees_suffix))
    
    # create temporary tibble to store values for loop iteration
    tmp_bootstrap_values_tibble <- tibble(loci = c(gene), shalrt_bootstrap = c(gene_tree$node.label))
    
    # exclude empty node labels
    tmp_bootstrap_values_tibble <- tmp_bootstrap_values_tibble %>% 
      filter(!shalrt_bootstrap == "")
    
    # extract bootstrap score from SH-aLRT score
    tmp_bootstrap_values_tibble <- tmp_bootstrap_values_tibble %>% 
      mutate(bootstrap = str_split_i(tmp_bootstrap_values_tibble$shalrt_bootstrap, pattern = "/", i = 2)) %>% 
      select(!shalrt_bootstrap)
    
    # append loop iteration values to all iteration tibble
    bootstrap_values_tibble <- rbind(bootstrap_values_tibble, tmp_bootstrap_values_tibble)
    
  }
  
  # drop empty first row
  bootstrap_values_tibble <- bootstrap_values_tibble %>% 
    filter(!loci == "")
  
  # convert bootstrap score to numeric
  bootstrap_values_tibble$bootstrap <- as.numeric(bootstrap_values_tibble$bootstrap)
  
  # calculate average bootstrap score per loci
  average_bootstrap <- bootstrap_values_tibble %>% 
    group_by(loci) %>% 
    summarise(average_bootstrap_score = mean(bootstrap),
              variance = var(bootstrap),
              standard_deviation = sd(bootstrap),
              standard_error = sd(bootstrap)/sqrt(length(bootstrap))) %>% 
    ungroup()
  
  #return(average_bootstrap)
  
  return(list(average_bootstrap_table = average_bootstrap, full_bootstrap_values = bootstrap_values_tibble))
  
}

# wrapper around getParent to get all node labels to root
getParent_all <- function(tree, node_label) {
  
  iteration <- 0
  
  # create empty vector to store parent nodes
  parent_nodes <- c()
  
  # arbitrarily set parent node to start loop
  parent_node <- 1
  
  while (!is.null(parent_node)) {
    
    # iteration counter
    iteration <- iteration + 1
    
    if (iteration == 1) {
      
      # if first iteration get parent of clade node label
      parent_node <- getParent(tree, node_label)
      
      # write node to vector
      if (!is.null(parent_node)) {
        parent_nodes <- c(parent_nodes, parent_node) 
      }
      
    } else {
      
      # if subsequent iteration get parent of last parent
      parent_node <- getParent(tree, parent_node)
      
      if (!is.null(parent_node)) {
        parent_nodes <- c(parent_nodes, parent_node) 
      }
      
    }
    
  }
  
  return(parent_nodes)
  
}

primer_vs_binding_sites <- function(fasta_path, sample_clade_lookup, primer_sequence, primer_name, reverse_complement = FALSE) {
  
  aligned_sequences <- read_csv(fasta_path, col_names = "fasta")
  
  # separate headers and sequence data
  aligned_sequences_headers <- aligned_sequences %>% 
    filter(str_detect(fasta, ">"))
  aligned_sequences_sequence <- aligned_sequences %>% 
    filter(!str_detect(fasta, ">"))
  
  # replace original sequence list
  aligned_sequences <- tibble(sample = aligned_sequences_headers[[1]], sequence = aligned_sequences_sequence[[1]])
  
  # remove > markers, _R_, jalview sequence information, and method information
  aligned_sequences <- aligned_sequences %>% 
    mutate(sample = case_when(!is.na(str_split_i(sample, ">_R_", 2)) ~ paste0(">", str_split_i(sample, ">_R_", 2)),
                              .default = sample)) %>% # replace >_R_ with >
    mutate(sample = str_split_i(sample, ">", 2)) %>% # remove >
    mutate(sample = str_split_i(sample, "/", 1)) %>% # remove / and anything after
    mutate(sample = strReverse(sample)) %>% # reverse string to enable removal of method info without splitting sample names
    mutate(sample = str_split_fixed(sample, "\\.", 3)[,3]) %>%
    mutate(sample = strReverse(sample)) %>% 
    rename(sample_ID = "sample")
  
  # add clade info
  aligned_sequences <- inner_join(aligned_sequences, sample_clade_lookup, join_by(sample_ID))
  
  # add primer sequence as row
  primer_sequence <- tibble(sample_ID = NA, sequence = primer_sequence, clade = NA)
  
  aligned_sequences <- rbind(primer_sequence, aligned_sequences)
  
  # get length of primer sequence
  primer_length <- nchar(primer_sequence$sequence[1])
  
  # split nucleotides into separate variables
  aligned_sequences <- cbind(aligned_sequences,str_split_fixed(aligned_sequences$sequence, "", primer_length))
  
  # separate primer sequence
  primer_sequence <- aligned_sequences %>% filter(is.na(sample_ID))
  aligned_sequences <- aligned_sequences %>% filter(!is.na(sample_ID))
  
  # pivot position names into position variable
  aligned_sequences <- pivot_longer(aligned_sequences,
                                    cols = -c(sample_ID, sequence, clade),
                                    names_to = "position",
                                    values_to = "nucleotide")
  
  primer_sequence <- primer_sequence %>% 
    pivot_longer(cols = -c(sample_ID, sequence, clade),
                 names_to = "position",
                 values_to = "original_primer_nucleotide") %>% 
    select(!c(sample_ID, sequence, clade)) %>% 
    # convert IUPAC nucleotides to string of eligible nucleotides
    mutate(primer_nucleotide = case_when(original_primer_nucleotide == "r" ~ "ag",
                                         original_primer_nucleotide == "y" ~ "ct",
                                         original_primer_nucleotide == "s" ~ "gc",
                                         original_primer_nucleotide == "w" ~ "at",
                                         original_primer_nucleotide == "k" ~ "gt",
                                         original_primer_nucleotide == "m" ~ "ac",
                                         original_primer_nucleotide == "b" ~ "cgt",
                                         original_primer_nucleotide == "d" ~ "agt",
                                         original_primer_nucleotide == "h" ~ "act",
                                         original_primer_nucleotide == "v" ~ "acg",
                                         original_primer_nucleotide == "n" ~ "atcg",
                                         .default = original_primer_nucleotide))
  
  # join primer sequence by position
  aligned_sequences <- left_join(aligned_sequences, primer_sequence, join_by(position))
  
  # change nucleotide to . if matches primer sequence
  aligned_sequences <- aligned_sequences %>% 
    mutate(nucleotide = case_when(str_detect(primer_nucleotide, nucleotide) ~ ".",
                                  .default = nucleotide))
  
  # make position a factor to control order
  aligned_sequences$position <- factor(aligned_sequences$position, levels = sort(as.numeric(unique(aligned_sequences$position))))
  
  # add primer sequence back as row in dataframe
  aligned_sequences <- rbind(aligned_sequences, primer_sequence %>% mutate(sample_ID = primer_name,
                                                                           clade = primer_name,
                                                                           sequence = NA,
                                                                           nucleotide = original_primer_nucleotide))
  
  # make clade and sample id factors - order by clade alphabetically except with primr on top
  aligned_sequences$clade <- factor(aligned_sequences$clade, levels = c(primer_name, unique((aligned_sequences %>% filter(clade != primer_name) %>% arrange(clade))$clade)))
  
  # sample id factor order follows clades
  aligned_sequences <- aligned_sequences %>% 
    arrange(clade, sample_ID) %>% 
    mutate(sample_ID = factor(sample_ID, levels = unique(sample_ID)))
  
  # create a degenerate primer variable
  aligned_sequences <- aligned_sequences %>% 
    mutate(degenerate_base = case_when(clade == primer_name & nchar(primer_nucleotide) == 2 ~ "2-fold degenerate",
                                       clade == primer_name & nchar(primer_nucleotide) == 3 ~ "3-fold degenerate",
                                       clade == primer_name & nchar(primer_nucleotide) == 4 ~ "4-fold degenerate",
                                       .default = "fixed"))
  
  # create a purine vs pyrimidine match variable (conserved purine and conserved pyrimidine are more tolerable mismatches)
  aligned_sequences <- aligned_sequences %>%
    mutate(purine_pyrimidine = case_when(clade != primer_name & str_detect(nucleotide, "a|g") ~ "purine",
                                         clade != primer_name & str_detect(nucleotide, "c|t") ~ "pyrimidine",
                                         .default = NA),
           primer_purine_pyrimidine = case_when(clade != primer_name & str_detect(original_primer_nucleotide, "a|g|r") ~ "purine",
                                                clade != primer_name & str_detect(original_primer_nucleotide, "c|t|y") ~ "pyrimidine",
                                                .default = NA),
           conserved_purine_pyrimidine = case_when(clade != primer_name & nucleotide != "." & purine_pyrimidine == primer_purine_pyrimidine ~ "conserved",
                                                   clade != primer_name & nucleotide != "." & purine_pyrimidine != primer_purine_pyrimidine ~ "not conserved",
                                                   .default = NA)
           )
  
  # identify the 3' primer end
  
  if (reverse_complement == FALSE) {
    # 3' end is on right
    primer_3_prime <- primer_sequence %>% 
      arrange(desc(as.numeric(position))) %>% 
      slice(1:5)
  } else {
    # 3' end is on LEFT - sequence written 5' to 3' but primer is reverse complement
    primer_3_prime <- primer_sequence %>% 
      arrange(as.numeric(position)) %>% 
      slice(1:5)
  }
  
  # plot data
  plot <- ggplot(aligned_sequences, y = sample_ID, x = position) +
    # indicate 3' end - rectangle defines by primer seq position and 5x 3' sites
    annotate("rect",
             ymin = length(unique(aligned_sequences$sample_ID)) - 1.25,
             ymax = length(unique(aligned_sequences$sample_ID)) + 0.75,
             xmin = min(as.numeric(primer_3_prime$position)) - 0.5,
             xmax = max(as.numeric(primer_3_prime$position)) + 0.5,
             fill = "darkgray") +
    # plot nucleotides
    geom_text(aes( y = sample_ID, x = position, label = nucleotide)) +
    # colour primer nucleotides by degeneracy
    geom_text(data = aligned_sequences %>% filter(clade == primer_name), aes(y = sample_ID, x = position, label = nucleotide, colour = degenerate_base)) +
    scale_y_discrete(labels = function(x) aligned_sequences$clade[match(x, aligned_sequences$sample_ID)],
                     limits=rev) +
    scale_colour_manual(
      values = c("2-fold degenerate" = "#FFC100", "3-fold degenerate" = "#FF7400", "4-fold degenerate" = "#FF0000", "fixed" = "black"),
      name = "Degenerate Primer Bases"
    ) +
    ggnewscale::new_scale_colour() + # facilitates separate scales
    # plot tolerable mismatches in a light grey
    geom_text(data = aligned_sequences %>% filter(clade != primer_name, !is.na(conserved_purine_pyrimidine)), aes(y = sample_ID, x = position, label = nucleotide, colour = conserved_purine_pyrimidine)) +
    scale_colour_manual(
      values = c("conserved" = "lightgrey", "not conserved" = "black"),
      name = "Conserved Purine or Pyrimidine\n(tolerable mismatch)"
    ) +
    theme_minimal() +
    ylab("clade") +
    ggtitle(paste("Unique clade haplotype variation against", primer_name))
  
  #return(primer_3_prime)
  return(plot)
  
}

# add division to list so can order by division then by clade to see if different divisions have different patterns
primer_vs_binding_sites_2 <- function(fasta_path, sample_clade_lookup, primer_sequence, primer_name, reverse_complement = FALSE) {
  
  aligned_sequences <- read_csv(fasta_path, col_names = "fasta")
  
  # separate headers and sequence data
  aligned_sequences_headers <- aligned_sequences %>% 
    filter(str_detect(fasta, ">"))
  aligned_sequences_sequence <- aligned_sequences %>% 
    filter(!str_detect(fasta, ">"))
  
  # replace original sequence list
  aligned_sequences <- tibble(sample = aligned_sequences_headers[[1]], sequence = aligned_sequences_sequence[[1]])
  
  # remove > markers, _R_, jalview sequence information, and method information
  aligned_sequences <- aligned_sequences %>% 
    mutate(sample = case_when(!is.na(str_split_i(sample, ">_R_", 2)) ~ paste0(">", str_split_i(sample, ">_R_", 2)),
                              .default = sample)) %>% # replace >_R_ with >
    mutate(sample = str_split_i(sample, ">", 2)) %>% # remove >
    mutate(sample = str_split_i(sample, "/", 1)) %>% # remove / and anything after
    mutate(sample = strReverse(sample)) %>% # reverse string to enable removal of method info without splitting sample names
    mutate(sample = str_split_fixed(sample, "\\.", 3)[,3]) %>%
    mutate(sample = strReverse(sample)) %>% 
    rename(sample_ID = "sample")
  
  ### EDITING TO ADD DIVISION ###
  
  # add clade info
  aligned_sequences <- left_join(aligned_sequences, sample_clade_lookup, join_by(sample_ID))
  
  # add primer sequence as row
  primer_sequence <- tibble(sample_ID = NA, sequence = primer_sequence, clade = NA, final_division = NA)
  
  aligned_sequences <- rbind(primer_sequence, aligned_sequences)
  
  # get length of primer sequence
  primer_length <- nchar(primer_sequence$sequence[1])
  
  # split nucleotides into separate variables
  aligned_sequences <- cbind(aligned_sequences,str_split_fixed(aligned_sequences$sequence, "", primer_length))
  
  # separate primer sequence
  primer_sequence <- aligned_sequences %>% filter(is.na(sample_ID))
  aligned_sequences <- aligned_sequences %>% filter(!is.na(sample_ID))
  
  # pivot position names into position variable
  aligned_sequences <- pivot_longer(aligned_sequences,
                                    cols = -c(sample_ID, sequence, clade, final_division),
                                    names_to = "position",
                                    values_to = "nucleotide")
  
  primer_sequence <- primer_sequence %>% 
    pivot_longer(cols = -c(sample_ID, sequence, clade, final_division),
                 names_to = "position",
                 values_to = "original_primer_nucleotide") %>% 
    select(!c(sample_ID, sequence, clade, final_division)) %>% 
    # convert IUPAC nucleotides to string of eligible nucleotides
    mutate(primer_nucleotide = case_when(original_primer_nucleotide == "r" ~ "ag",
                                         original_primer_nucleotide == "y" ~ "ct",
                                         original_primer_nucleotide == "s" ~ "gc",
                                         original_primer_nucleotide == "w" ~ "at",
                                         original_primer_nucleotide == "k" ~ "gt",
                                         original_primer_nucleotide == "m" ~ "ac",
                                         original_primer_nucleotide == "b" ~ "cgt",
                                         original_primer_nucleotide == "d" ~ "agt",
                                         original_primer_nucleotide == "h" ~ "act",
                                         original_primer_nucleotide == "v" ~ "acg",
                                         original_primer_nucleotide == "n" ~ "atcg",
                                         .default = original_primer_nucleotide))
  
  # join primer sequence by position
  aligned_sequences <- left_join(aligned_sequences, primer_sequence, join_by(position))
  
  # change nucleotide to . if matches primer sequence
  aligned_sequences <- aligned_sequences %>% 
    mutate(nucleotide = case_when(str_detect(primer_nucleotide, nucleotide) ~ ".",
                                  .default = nucleotide))
  
  # make position a factor to control order
  aligned_sequences$position <- factor(aligned_sequences$position, levels = sort(as.numeric(unique(aligned_sequences$position))))
  
  # add primer sequence back as row in dataframe
  aligned_sequences <- rbind(aligned_sequences, primer_sequence %>% mutate(sample_ID = primer_name,
                                                                           clade = primer_name,
                                                                           sequence = NA,
                                                                           final_division = NA,
                                                                           nucleotide = original_primer_nucleotide))
  
  # make clade and sample id factors - order by DIVISION, then by clade alphabetically except with primer on top
  aligned_sequences$clade <- factor(aligned_sequences$clade, levels = c(primer_name, unique((aligned_sequences %>% filter(clade != primer_name) %>% arrange(final_division, clade))$clade)))
  
  # sample id factor order follows DIVISION then clades
  aligned_sequences <- aligned_sequences %>% 
    arrange(sample_ID, clade, final_division) %>% 
    mutate(sample_ID = factor(sample_ID, levels = unique(sample_ID)))
  
  # create a degenerate primer variable
  aligned_sequences <- aligned_sequences %>% 
    mutate(degenerate_base = case_when(clade == primer_name & nchar(primer_nucleotide) == 2 ~ "2-fold degenerate",
                                       clade == primer_name & nchar(primer_nucleotide) == 3 ~ "3-fold degenerate",
                                       clade == primer_name & nchar(primer_nucleotide) == 4 ~ "4-fold degenerate",
                                       .default = "fixed"))
  
  # create a purine vs pyrimidine match variable (conserved purine and conserved pyrimidine are more tolerable mismatches)
  aligned_sequences <- aligned_sequences %>%
    mutate(purine_pyrimidine = case_when(clade != primer_name & str_detect(nucleotide, "a|g") ~ "purine",
                                         clade != primer_name & str_detect(nucleotide, "c|t") ~ "pyrimidine",
                                         .default = NA),
           primer_purine_pyrimidine = case_when(clade != primer_name & str_detect(original_primer_nucleotide, "a|g|r") ~ "purine",
                                                clade != primer_name & str_detect(original_primer_nucleotide, "c|t|y") ~ "pyrimidine",
                                                .default = NA),
           conserved_purine_pyrimidine = case_when(clade != primer_name & nucleotide != "." & purine_pyrimidine == primer_purine_pyrimidine ~ "conserved",
                                                   clade != primer_name & nucleotide != "." & purine_pyrimidine != primer_purine_pyrimidine ~ "not conserved",
                                                   .default = NA)
    )
  
  # identify the 3' primer end
  
  if (reverse_complement == FALSE) {
    # 3' end is on right
    primer_3_prime <- primer_sequence %>% 
      arrange(desc(as.numeric(position))) %>% 
      slice(1:5)
  } else {
    # 3' end is on LEFT - sequence written 5' to 3' but primer is reverse complement
    primer_3_prime <- primer_sequence %>% 
      arrange(as.numeric(position)) %>% 
      slice(1:5)
  }
  
  # plot data
  plot <- ggplot(aligned_sequences, y = sample_ID, x = position) +
    # indicate 3' end - rectangle defines by primer seq position and 5x 3' sites
    annotate("rect",
             ymin = length(unique(aligned_sequences$sample_ID)) - 1.25,
             ymax = length(unique(aligned_sequences$sample_ID)) + 0.75,
             xmin = min(as.numeric(primer_3_prime$position)) - 0.5,
             xmax = max(as.numeric(primer_3_prime$position)) + 0.5,
             fill = "darkgray") +
    # plot nucleotides
    geom_text(aes( y = sample_ID, x = position, label = nucleotide)) +
    # colour primer nucleotides by degeneracy
    geom_text(data = aligned_sequences %>% filter(clade == primer_name), aes(y = sample_ID, x = position, label = nucleotide, colour = degenerate_base)) +
    scale_y_discrete(labels = function(x) aligned_sequences$clade[match(x, aligned_sequences$sample_ID)],
                     limits=rev) +
    scale_colour_manual(
      values = c("2-fold degenerate" = "#FFC100", "3-fold degenerate" = "#FF7400", "4-fold degenerate" = "#FF0000", "fixed" = "black"),
      name = "Degenerate Primer Bases"
    ) +
    ggnewscale::new_scale_colour() + # facilitates separate scales
    # plot tolerable mismatches in a light grey
    geom_text(data = aligned_sequences %>% filter(clade != primer_name, !is.na(conserved_purine_pyrimidine)), aes(y = sample_ID, x = position, label = nucleotide, colour = conserved_purine_pyrimidine)) +
    scale_colour_manual(
      values = c("conserved" = "lightgrey", "not conserved" = "black"),
      name = "Conserved Purine or Pyrimidine\n(tolerable mismatch)"
    ) +
    theme_minimal() +
    ylab("clade") +
    ggtitle(paste("Unique clade haplotype variation against", primer_name))
  
  #return(primer_3_prime)
  return(plot)
  
}