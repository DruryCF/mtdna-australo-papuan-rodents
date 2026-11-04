#!/usr/bin/Rscript

library(tidyverse)
library(sjmisc)

OUTPREFIX <- Sys.getenv("OUTPREFIX_VALUE_EXPORT")
MT_REFERENCE <- Sys.getenv("MT_REFERENCE_VALUE_EXPORT")

samples <- read_lines(paste0("./", OUTPREFIX, "tmp/sampleIDs.txt"))
methods <- read_lines(paste0("./", OUTPREFIX, "tmp/benchmark_methods.txt"))
expected_genes <- read_csv(paste0("./", OUTPREFIX, "tmp/expected_genes.csv"), show_col_types = F)

# if exclude genes read depth file exists, read lines
path_exclude_genes_read_depths <- paste0("./", OUTPREFIX, "tmp/exclude_genes_read_depths.txt")

if (file.exists(path_exclude_genes_read_depths)) {
  exclude_read_depth_genes <- read_lines(path_exclude_genes_read_depths)
}

# calculate reference bp

# replace path to present working directory with a period
mt_reference_path <- paste0(".", str_remove(MT_REFERENCE, pattern = getwd()))

# read lines of reference fasta
mt_reference_fasta <- read_lines(mt_reference_path)

# remove header and confirm uppercase nucleotides
mt_reference_no_header <- toupper(mt_reference_fasta[! str_detect(mt_reference_fasta, ">")])

# remove all white space (if any)
mt_reference_no_header <- na_if(as.vector(mt_reference_no_header), c(""))
mt_reference_no_header <- na_if(as.vector(mt_reference_no_header), c(" "))
mt_reference_no_header <- mt_reference_no_header[!is.na(mt_reference_no_header)]
mt_reference_no_header <- str_squish(mt_reference_no_header)

# get num_bp
ref_bp <- sum(str_count(mt_reference_no_header))

for (s in 1:length(samples)) {
  
  for (m in 1:length(methods)) {
    
    #store path to input fasta file
    path_input_fasta <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".input.fasta")
    
    #check if input fasta file exists, if not write NA for all relevant stats
    if (file.exists(path_input_fasta)) {
      
      # load fasta
      input_fasta <- read_lines(path_input_fasta)
      
      # get stats
      
      # subset non-headers and convert/confirm all nucleotides to uppercase
      non_header_strings <- toupper(input_fasta[! str_detect(input_fasta, ">")])
      
      # count number of bp
      num_bp <- sum(str_count(non_header_strings))
      
      # count number of Ns
      
      num_Ns <- sum(str_count(non_header_strings, "N"))

      # calculate identified nucleotides
      num_bp_minus_Ns <- num_bp - num_Ns
      
    } else {
      
      # write NA for relevant stats
      num_Ns <- "NA"
      num_bp <- "NA"
      num_bp_minus_Ns <- "NA"

    }
    
    #store path to contig summary file
    path_contig_summary <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".contig_summary.csv")
    
    #check if contig summary file exists, if not write NA for all relevant stats
    if (file.exists(path_contig_summary)) {
      
      # load csv
      contig_summary <- read_csv(path_contig_summary, show_col_types = F)
      
      # get stats
      num_contigs <- nrow(contig_summary)
      
      # get circularity
      circularity_all_contigs <- contig_summary$Circularity
      circularity <- NULL
      
      # check if contains yes
      if (str_contains(circularity_all_contigs, "yes")) {
        # if single row mark as yes else as inconclusive (covers multi yes and mixed yes no)
        if (nrow(contig_summary) == 1) {
          circularity <- "yes"
        } else {
          circularity <- "inconclusive"
        }
      }
      
      # check if contains no, if not and NULL mark as NA
      if (str_contains(circularity_all_contigs, "no")) {
        # confirm circularity not already set (prevents overwriting of inconclusive mixed yes no) and mark as no
        if (is.null(circularity)) {
          circularity <- "no"
        }
      } else if (is.null(circularity)) {
        circularity <- "NA" 
      }
      
    } else {
      
      # write NA for relevant stats
      num_contigs <- "NA"
      circularity <- "NA" 
      
    }
    
    #store path to gene summary file
    path_gene_summary <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".gene_summary.csv")
    
    #check if gene summary file exists, if not write NA for all relevant stats
    if (file.exists(path_gene_summary)) {
      
      # load csv
      gene_summary <- read_csv(path_gene_summary, show_col_types = F)
      
      # get stats
      
      # check if any gene is detected more than expected number of times
      
      potential_contamination <- "not_detected"
      
      for (g in 1:nrow(expected_genes)) {
        
        gene_check <- expected_genes[g,]
        
        # subset gene_summary by gene and return number of rows
        gene_num_detected <- nrow(gene_summary[str_detect(gene_summary$Gene_name, gene_check$gene),])
        
        # check if gene partial matches to any other expected gene
        if(length(expected_genes$gene[str_detect(expected_genes$gene, gene_check$gene)]) > 1){
          
          # create match detected tibble
          match_detected <- gene_summary[str_detect(gene_summary$Gene_name, gene_check$gene),]
          
          # subset non_exact matches from expected genes
          detected_matches <- expected_genes$gene[str_detect(expected_genes$gene, gene_check$gene)]
          incorrect_matches <- detected_matches[! detected_matches == gene_check$gene]
          
          # find partial matches of non_interest genes in detected rows
          for (d in length(incorrect_matches)) {
            non_match_indices <- which(str_detect(match_detected$Gene_name, incorrect_matches[d]))
            remaining_indices <- setdiff(1:(nrow(match_detected)), non_match_indices)
            
            # subset detected matches by remaining row indices
            match_detected <- match_detected[remaining_indices,]
          }
          # calculate new number of genes detected
          gene_num_detected <- nrow(match_detected)
        }
        
        #print(paste(gene_check$gene, "detected:", gene_num_detected, "expected:", gene_check$num_copies))
        
        # if number of rows greater than expected then mark as potential contamination and break loop
        if (gene_num_detected > gene_check$num_copies) {
          potential_contamination <- "detected"
          break
        }
      }
      
      # count number of missing stop codons -> distinguish by direction because Ns at start also represented
      reverse_genes <- gene_summary %>% 
        filter(Direction == "-")
      forward_genes <- gene_summary %>% 
        filter(Direction == "+")
      missing_stop_codons <- sum(sum(str_count(forward_genes$End, ">")), sum(str_count(reverse_genes$Start, "<")))
      
    } else {
      
      # write NA for relevant stats
      potential_contamination <- "NA"
      missing_stop_codons <- "NA"
      
    }
    
    #store path to detected_gene_type summary file
    path_detected_gene_type_summary <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".detected_gene_type_summary.csv")
    
    #check if detected_gene_type summary file exists, if not write NA for all relevant stats
    if (file.exists(path_detected_gene_type_summary)) {
      
      # load csv
      detected_gene_type_summary <- read_csv(path_detected_gene_type_summary, show_col_types = F)
      
      # get stats
      genes_detected <- sum(detected_gene_type_summary$total_detected)
      CDS_genes <- detected_gene_type_summary$total_detected[detected_gene_type_summary$gene_type == "CDS"]
      tRNA_genes <- detected_gene_type_summary$total_detected[detected_gene_type_summary$gene_type == "tRNA"]
      rRNA_genes <- detected_gene_type_summary$total_detected[detected_gene_type_summary$gene_type == "rRNA"]
      
    } else {
      
      # write NA for relevant stats
      genes_detected <- "NA"
      CDS_genes <- "NA"
      tRNA_genes <- "NA"
      rRNA_genes <- "NA"
      
    }
    
    # calculate dependant stats
    if (! num_bp == "NA") {
      #calculate bp ber ref bp
      bp_per_refbp <- num_bp/ref_bp

      bp_minus_Ns_per_refbp <- num_bp_minus_Ns/ref_bp
      
      if (! num_Ns == "NA") {
        #calculate Ns per bp
        Ns_per_bp <- num_Ns/num_bp
      } else {
        Ns_per_bp <- "NA"
        bp_minus_Ns_per_refbp <- "NA"
      }
      
    } else {
      bp_per_refbp <- "NA"
      bp_minus_Ns_per_refbp <- "NA"
      Ns_per_bp <- "NA"
    }
    
    #store path to read_depth file
    path_read_depths <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".read_depth.txt")
    
    #check if read_depth file exists, if not write NA for all relevant stats
    if (file.exists(path_read_depths)) {
      
      # load txt file
      read_depths <- read_lines(path_read_depths)
      
      # convert to tibble
      read_depths_str_split <- str_split_fixed(read_depths, pattern = " ", n = 4)
      read_depths_tibble <- tibble(contig = read_depths_str_split[,1], bp_1 = read_depths_str_split[,2], bp_2 = read_depths_str_split[,3], read_depth = read_depths_str_split[,4])
      
      # get stats
      avg_read_depth <- mean(as.numeric(read_depths_tibble$read_depth))
      med_read_depth <- median(as.numeric(read_depths_tibble$read_depth))
      min_read_depth <- min(as.numeric(read_depths_tibble$read_depth))
      max_read_depth <- max(as.numeric(read_depths_tibble$read_depth))
      
      if (file.exists(path_exclude_genes_read_depths)) {
        
        # check if gene_summary file exists -> required for gene exclusion
        if (file.exists(path_gene_summary)) {
          
          # store gene summary information of genes to be excluded from read depths
          exclude_read_depths <- gene_summary %>% 
            filter(Gene_name %in% exclude_read_depth_genes)
          
          # if genes to be excluded exist
          if (nrow(exclude_read_depths) > 0) {
            
            # remove indicators of missing stops (>, <) from exclude list
            exclude_read_depths$Start <- str_replace_all(exclude_read_depths$Start, ">", "")
            exclude_read_depths$Start <- str_replace_all(exclude_read_depths$Start, "<", "")
            exclude_read_depths$End <- str_replace_all(exclude_read_depths$End, ">", "")
            exclude_read_depths$End <- str_replace_all(exclude_read_depths$End, "<", "")
            
            # create tmp read depths tibble
            tmp_read_depths_tibble <- read_depths_tibble
            
            # sort contig names by length -> how ordered in read depths file
            ordered_contigs  <- contig_summary %>% arrange(`Length(bp)`)
            # store unique contig names
            contig_names <- ordered_contigs$`#Seq_id`
            
            # store all unique read depth contig names
            read_depth_contig_names <- unique(read_depths_tibble$contig)
            
            # loop through entries in genes to be excluded
            for (n in 1:nrow(exclude_read_depths)) {
              
              # store start and end bp
              start_bp <- exclude_read_depths[n,]$Start
              end_bp <- exclude_read_depths[n,]$End
              
              # if more than one contig, map contig name to name in read depths file
              if (length(contig_names) > 1) {
                # find contig name index in unique names list
                contig_name <- exclude_read_depths[n,]$`#Seq_id`
                
                contig_name_index <- which(contig_names == contig_name)
                
                # relate to name in read depths file
                read_depth_contig_name <- read_depth_contig_names[contig_name_index]
                
                # remove read_depth rows matching contig name and bp range
                
                # store rows not matching contig name
                read_depths_tibble_not_contig_match <- tmp_read_depths_tibble %>% 
                  filter(! contig == read_depth_contig_name)
                
                # store rows matching contig name, not matching excluded bp
                read_depths_tibble_contig_match <- tmp_read_depths_tibble %>%
                  filter(contig == read_depth_contig_name) %>% 
                  filter(! bp_1 %in% c(start_bp:end_bp))
                
                tmp_read_depths_tibble <- rbind(read_depths_tibble_not_contig_match, read_depths_tibble_contig_match)
              }
              
              # if only one contig, just look for bp
              if (length(contig_names) == 1) {
                
                tmp_read_depths_tibble <- tmp_read_depths_tibble %>%
                  filter(! bp_1 %in% c(start_bp:end_bp))
              }
            }       
            
            avg_read_depth_excluded_genes <- mean(as.numeric(tmp_read_depths_tibble$read_depth))
            med_read_depth_excluded_genes <- median(as.numeric(tmp_read_depths_tibble$read_depth))
            min_read_depth_excluded_genes <- min(as.numeric(tmp_read_depths_tibble$read_depth))
            max_read_depth_excluded_genes <- max(as.numeric(tmp_read_depths_tibble$read_depth))
            
          } else {
            
            # if genes to be excluded do not exist, calculate from full dataset
            avg_read_depth_excluded_genes <- mean(as.numeric(read_depths_tibble$read_depth))
            med_read_depth_excluded_genes <- median(as.numeric(read_depths_tibble$read_depth))
            min_read_depth_excluded_genes <- min(as.numeric(read_depths_tibble$read_depth))
            max_read_depth_excluded_genes <- max(as.numeric(read_depths_tibble$read_depth))
            
          }
          
        } else {
          
          # if gene summary does not exist
          avg_read_depth_excluded_genes <- "NA"
          med_read_depth_excluded_genes <- "NA"
          min_read_depth_excluded_genes <- "NA"
          max_read_depth_excluded_genes <- "NA"
          
        } 
        
      } else {
        
        # if exclude genes read depths does not exist
        avg_read_depth_excluded_genes <- "NA"
        med_read_depth_excluded_genes <- "NA"
        min_read_depth_excluded_genes <- "NA"
        max_read_depth_excluded_genes <- "NA"
        
      }
      
    } else {
      
      # if read depths file does not exist
      # write NA for relevant stats
      avg_read_depth <- "NA"
      med_read_depth <- "NA"
      min_read_depth <- "NA"
      max_read_depth <- "NA"
      
      if (file.exists(path_exclude_genes_read_depths)) {
        avg_read_depth_excluded_genes <- "NA"
        med_read_depth_excluded_genes <- "NA"
        min_read_depth_excluded_genes <- "NA"
        max_read_depth_excluded_genes <- "NA"
      }
    }

    #store path to translated fasta file
    path_translated_fasta <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".translated.fasta")

    #check if input fasta file exists, if not write NA for all relevant stats
    if (file.exists(path_translated_fasta)) {

	# load fasta
	translated_fasta <- read_lines(path_translated_fasta)

	# get stats

	# subset non-headers
	non_header_strings <- translated_fasta[! str_detect(translated_fasta, ">")]

	# count number of asterisks
	num_internal_stops <- sum(str_count(non_header_strings, "\\*"))

    } else {

	# write NA for relevant stats
	    num_internal_stops <- "NA"

    }
      
    # compile stats as comma delimited strings
    whole_genome_stats <- paste(samples[s], methods[m], num_Ns, num_bp, num_bp_minus_Ns, Ns_per_bp, bp_per_refbp, bp_minus_Ns_per_refbp, genes_detected, CDS_genes, rRNA_genes, tRNA_genes, circularity, num_contigs, potential_contamination, missing_stop_codons, num_internal_stops, avg_read_depth, med_read_depth, min_read_depth, max_read_depth, sep = ",")
    expected_commas <- 20
    
    if (file.exists(path_exclude_genes_read_depths)) {
      whole_genome_stats <- paste(whole_genome_stats, avg_read_depth_excluded_genes, med_read_depth_excluded_genes, min_read_depth_excluded_genes, max_read_depth_excluded_genes, sep = ",")
      expected_commas <- 24
    }
    
    # store NA string for expected number of commas
    na_string <- NULL
    
    for (c in 2:expected_commas) {
      na_string <- paste0(na_string, ",NA")
    }
    
    # check strings have expected number of commas
    # if yes, append strings to appropriate stats .csv file
    if (str_count(whole_genome_stats, ",") == expected_commas) {
      write_lines(whole_genome_stats, paste0("./benchmark_stats/", OUTPREFIX, "whole_genome_criteria.csv"), append = T)
    }
    # if longer, write NA for all entries and report error
    if (str_count(whole_genome_stats, ",") > expected_commas) {
      
      whole_genome_stats <- paste0(samples[s], ",", methods[m], na_string)
      write_lines(whole_genome_stats, paste0("./benchmark_stats/", OUTPREFIX, "whole_genome_criteria.csv"), append = T)
      
      print(paste("ERROR: more stats compiled than expected for", samples[s], ",", methods[m], sep = " "))
    }
    # if shorter, write NA for all entries and report error
    
    if (str_count(whole_genome_stats, ",") < expected_commas) {
      whole_genome_stats <- paste0(samples[s], ",", methods[m], na_string)
      write_lines(whole_genome_stats, paste0("./benchmark_stats/", OUTPREFIX, "whole_genome_criteria.csv"), append = T)
      
      print(paste("ERROR: fewer stats compiled than expected for", samples[s], ",", methods[m], sep = " "))
    }
    
  }
  
}
