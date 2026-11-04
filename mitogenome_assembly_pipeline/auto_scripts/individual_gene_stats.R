#!/usr/bin/Rscript

library(tidyverse)
library(sjmisc)

#print(paste("WORKING DIR:", getwd()))

OUTPREFIX <- Sys.getenv("OUTPREFIX_VALUE_EXPORT")
MT_REFERENCE_ANNOTATED <- Sys.getenv("MT_REFERENCE_ANNOTATED_VALUE_EXPORT")

samples <- read_lines(paste0("./", OUTPREFIX, "tmp/sampleIDs.txt"))
methods <- read_lines(paste0("./", OUTPREFIX, "tmp/benchmark_methods.txt"))
expected_genes <- read_csv(paste0("./", OUTPREFIX, "tmp/expected_genes.csv"), show_col_types = F)

# retrieve annotated reference and extract gene num bp

# replace path to present working directory with a period
mt_reference_annotated_path <- paste0(".", str_remove(MT_REFERENCE_ANNOTATED, pattern = getwd()))

# read lines of reference fasta
mt_reference_annotated_fasta <- read_lines(mt_reference_annotated_path)

# store header line indices
annotated_header_lines <- which(str_detect(mt_reference_annotated_fasta, ">"))

# create tibble of reference gene lengths
reference_gene_bp <- tibble(gene="",num_bp="")

# loop through genes

for (g in 1:nrow(expected_genes)) {
  
  # get current gene
  current_gene <- expected_genes[g,]
  
  # find location of header line that matches gene name -> case insensitive
  gene_header <- which(str_detect(tolower(mt_reference_annotated_fasta), tolower(current_gene$gene)))
  
  # check if gene partial matches to any other expected gene
  if(length(expected_genes$gene[str_detect(expected_genes$gene, current_gene$gene)]) > 1){
    
    # subset non_exact matches from expected genes
    detected_matches <- expected_genes$gene[str_detect(expected_genes$gene, current_gene$gene)]
    incorrect_matches <- detected_matches[! detected_matches == current_gene$gene]
    
    # find partial matches of non_interest genes in detected rows
    for (d in length(incorrect_matches)) {
      non_match_indices <- which(str_detect(tolower(mt_reference_annotated_fasta), tolower(incorrect_matches[d])))
      remaining_indices <- setdiff(gene_header, non_match_indices)
      
      # subset detected matches by remaining indices
      gene_header <- gene_header[gene_header %in% remaining_indices]
    }
  }
  
  # relate index back to all header lines and retrieve index of next header
  next_gene_header <- annotated_header_lines[(which(annotated_header_lines %in% gene_header) + 1)]
  
  # subset strings from matching header to next gene header -> need to loop if multiple headers
  for (h in 1:length(gene_header)) {
    
    if (! is.na(next_gene_header[h])) {
      gene_sequence <- mt_reference_annotated_fasta[(gene_header[h] + 1):(next_gene_header[h] - 1)]
    } else {
      gene_sequence <- mt_reference_annotated_fasta[(gene_header[h] + 1):length(mt_reference_annotated_fasta)]
    }
    
    # count number of nucleotides
    ref_gene_bp <- sum(str_count(gene_sequence))
    
    # create tibble
    gene_ref_bp <- tibble(gene=current_gene$gene,num_bp=ref_gene_bp)
    
    # append tibble to full tibble
    
    reference_gene_bp <- rbind(reference_gene_bp, gene_ref_bp)
    
  }
}

# remove initial empty row and change num_bp to numeric class
reference_gene_bp <- reference_gene_bp[2:nrow(reference_gene_bp),]
reference_gene_bp <- reference_gene_bp %>% mutate(across(2, as.numeric))

# retrieve stats

#print("error after point 1")

for (s in 1:length(samples)) {
  
  for (m in 1:length(methods)) {
    
    # check file paths and write NA for relevant stats when files do not exist
    
    # store paths to files
    
    path_annotated_fasta <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".annotated.fasta")
    path_gene_summary <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".gene_summary.csv")
    path_contig_summary <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".contig_summary.csv")
    path_read_depths <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".read_depth.txt")
    path_translated_fasta <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".translated.fasta")
    
    # check existence and write NA where relevant
    
    if (file.exists(path_gene_summary)) {
      # load csv
      gene_summary <- read_csv(path_gene_summary, show_col_types = F)

      contig_summary <- read_csv(path_contig_summary, show_col_types = F)
      
      # set existence to true
      gene_summary_exists <- TRUE
    } else {
      # write NA for relevant stats
      stop_codon_detected <- "NA"
      
      # set existence to false
      gene_summary_exists <- FALSE
    }
    
    if (file.exists(path_read_depths)) {
      # load txt file
      read_depths <- read_lines(path_read_depths)
      
      # convert to tibble
      read_depths_str_split <- str_split_fixed(read_depths, pattern = " ", n = 4)
      read_depths_tibble <- tibble(contig = read_depths_str_split[,1], bp_1 = read_depths_str_split[,2], bp_2 = read_depths_str_split[,3], read_depth = read_depths_str_split[,4])
      
      # set existence to true
      read_depth_exists <- TRUE
    } else {
      # write NA for relevant stats
      avg_read_depth <- "NA"
      med_read_depth <- "NA"
      min_read_depth <- "NA"
      max_read_depth <- "NA"
      
      # set existence to false
      read_depth_exists <- FALSE
    }

    if (file.exists(path_translated_fasta)) {
      # load txt file
      translated_fasta <- read_lines(path_translated_fasta)

      # set existence to true
      translated_fasta_exists <- TRUE
    } else {
      # write NA for relevant stats
      num_internal_stops <- "NA"

      # set existence to false
      translated_fasta_exists <- FALSE
    }
    
    if (file.exists(path_annotated_fasta)) {
      # load fasta
      annotated_fasta <- read_lines(path_annotated_fasta)
      
      # get stats
      
      # store header line indices of annotated fasta
      annotated_header_lines <- which(str_detect(annotated_fasta, ">"))
      
      # if translated fasta exists, store header indices
      if (file.exists(path_translated_fasta)) {
	      translated_header_lines <- which(str_detect(translated_fasta, ">"))
      }

      for (g in 1:nrow(expected_genes)) {
        
        #get current gene
        current_gene <- expected_genes[g,]
        
        # find location of header line that matches gene name -> case insensitive
        gene_header <- which(str_detect(tolower(annotated_fasta), tolower(current_gene$gene)))

#print("error after point 2")

	# if translated fasta exists, repeat header retrieval -> appending semicolon to expected gene prevents partial matches like ND4 to ND4L
	#if (file.exists(path_translated_fasta)) {
	#	translated_gene_header <- which(str_detect(tolower(translated_fasta), paste0(tolower(current_gene$gene), ";")))
	#}
        
        # if no header found, continue to next gene
        if (is.na(gene_header[1])) {
          next
        }

	# if no header found in translated fasta - will occur with non-CDS genes - write NA for num_internal stops
	#if (file.exists(path_translated_fasta)) {
	#	if (is.na(translated_gene_header[1])) {
	#		num_internal_stops <- "NA"
	#	}
	#}
        
        # check if gene partial matches to any other expected gene
        if(length(expected_genes$gene[str_detect(expected_genes$gene, current_gene$gene)]) > 1){
          
          # subset non_exact matches from expected genes
          detected_matches <- expected_genes$gene[str_detect(expected_genes$gene, current_gene$gene)]
          incorrect_matches <- detected_matches[! detected_matches == current_gene$gene]
          
          # find partial matches of non_interest genes in detected rows
          for (d in length(incorrect_matches)) {
            non_match_indices <- which(str_detect(tolower(annotated_fasta), tolower(incorrect_matches[d])))
            remaining_indices <- setdiff(gene_header, non_match_indices)
            
            # subset detected matches by remaining indices
            gene_header <- gene_header[gene_header %in% remaining_indices]
          }
        }
        
        # if no header found, continue to next gene
        if (is.na(gene_header[1])) {
          next
        }
        
        # count headers for number of gene copies
        num_gene_copies <- length(gene_header)
        
        # relate index back to all header lines and retrieve index of next header
        next_gene_header <- annotated_header_lines[(which(annotated_header_lines %in% gene_header) + 1)]
        
	# relate index back to all header lines and retrieve index of next header - for translated fasta
	#if (file.exists(path_translated_fasta)) {
	#	if (!is.na(translated_gene_header[1])) {
	#		next_translated_gene_header <- translated_header_lines[(which(translated_header_lines %in% translated_gene_header) + 1)]
	#	}
	#}

        # subset strings from matching header to next gene header -> need to loop if multiple headers
        for (h in 1:length(gene_header)) {
          if (! is.na(next_gene_header[h])) {
            gene_sequence <- annotated_fasta[(gene_header[h] + 1):(next_gene_header[h] - 1)]
          } else {
            gene_sequence <- annotated_fasta[(gene_header[h] + 1):length(annotated_fasta)]
          }

#	print("error after point 3")

	  # if translated fasta exists, match exactly to gene header and relate index back
	  if (file.exists(path_translated_fasta)) {
		  # confirm headers exist in translation file -> will not exist if no cds genes (file may or may not exist, unsure)
		  if (!is.na(translated_header_lines[1])) {
			  # get header line indice for exact CONTIG AND POSITION match to current gene header - note length is different between fasta and amino acid seq.
			  #print("error after point 4")
			  # raw header text
			  gene_header_text <- annotated_fasta[gene_header[h]]
			  # contig
			  gene_header_contig <- str_split_i(gene_header_text, ";", i=1)
			  # position
			  gene_header_position <- str_split_i(gene_header_text, ";", i=4)
			  # search for header in transladed fasta that matches both contig and position
			  translated_gene_header <- which(str_detect(translated_fasta, fixed(gene_header_contig)) & str_detect(translated_fasta, fixed(gene_header_position)))

			  # check if translated gene header is empty -> can happen when headers in gene fasta vs translated fasta have discrepencies e.g. length recorded wrong in translated file -> CORRECTED ERROR
			  if (length(translated_gene_header) == 0) {
				 num_internal_stops <- "NA"
			  } else { 

#print("error after point 5")

			  	# confirm header matches to translated header -> will not if non cds gene
			  	if (!is.na(translated_gene_header)) { 

				  	# get line indice for next header in fasta file by relating back to header line indices
				  	next_translated_gene_header <- translated_header_lines[(which(translated_header_lines %in% translated_gene_header) + 1)]

			  	  	# subset strings from matching header to next gene header
		  		  	if (! is.na(next_translated_gene_header)) {
			  			translated_gene_sequence <- translated_fasta[(translated_gene_header + 1):(next_translated_gene_header - 1)]
			  	  	} else {
		  		  		translated_gene_sequence <- translated_fasta[(translated_gene_header + 1):length(translated_fasta)]
	  		  	  	}

		 		  	# count number of asterisks in translated gene sequence
		  		  	num_internal_stops <- sum(str_count(translated_gene_sequence, "\\*"))

		  	   	} else {
				   	# set num_internal_stops to NA if no headers in file
				   	num_internal_stops <- "NA"
			   	}
			  }
		  } else {
			  # set num_internal_stops to NA if header does not match to header in translation file
			  num_internal_stops <- "NA"
		  }

	  }
	  
          
          # count number of nucleotides
          num_bp <- sum(str_count(gene_sequence))
          
          # count number of Ns
          num_Ns <- sum(str_count(toupper(gene_sequence), "N"))

	  # calculate number of nucleotides minus Ns
	  num_bp_minus_Ns <- num_bp - num_Ns
         
#print("error after point 6")

          # calculate dependant stats
          if (! num_bp == "NA") {
            
            # retrieve ref_bp -> will be mean of duplicate genes i.e. trnL and trnS
            current_ref_gene <- reference_gene_bp[reference_gene_bp$gene %in% current_gene$gene,]
            ref_bp <- mean(current_ref_gene$num_bp)
            
            # calculate bp ber ref bp
            bp_per_refbp <- num_bp/ref_bp

	    # calculate bp minus Ns per ref bp
            bp_minus_Ns_per_refbp <- num_bp_minus_Ns/ref_bp
            
            if (! num_Ns == "NA") {
              # calculate Ns per bp
              Ns_per_bp <- num_Ns/num_bp
            } else {
              Ns_per_bp <- "NA"
            }
          } else {
            bp_per_refbp <- "NA"
            Ns_per_bp <- "NA"
	    bp_minus_Ns_per_refbp <- "NA"
          }
          
          if (gene_summary_exists == TRUE) {
            # subset gene summary by gene
#print("error after point 7")
            gene_summary_subset <- gene_summary[str_detect(gene_summary$Gene_name, current_gene$gene),]
            
            # select observation matching sequence header
            specific_gene_summary <- gene_summary_subset[h,]

#print("error after point 8")
#print(paste(samples[s], methods[m], current_gene$gene))
#print(specific_gene_summary)

            # look for missing stop codon -> num_missing_stops should be 1 or 0
            
            if (specific_gene_summary$Direction == "-") {
              
              # if direction is negative then stop codon is on start
              num_missing_stops <- sum(str_count(specific_gene_summary$Start, "<"))
            }
            if (specific_gene_summary$Direction == "+") {
              
              # if direction is negative then stop codon is on start
              num_missing_stops <- sum(str_count(specific_gene_summary$End, ">"))
            }
            
            if (num_missing_stops == 0) {
              stop_codon_detected <- "yes"
            } else if (num_missing_stops == 1) {
              stop_codon_detected <- "no"
            } else {
              stop_codon_detected <- "NA"
              print(paste("ERROR: number of missing stop codons for", samples[s], methods[m], current_gene$gene, "is not equal to 1 or 0, equals:", num_missing_stops))
            }

	    # sort contig names by length -> how ordered on read depths file and how numbered in tree files
            ordered_contigs  <- contig_summary %>% arrange(`Length(bp)`)
            
	    # subset to contig name, and length, add contig number (from zero)
            contig_names_and_numbers <- ordered_contigs %>% select(`#Seq_id`,`Length(bp)`) %>% mutate(contig_number = row_number() - 1)

	    # get specific contig information for current gene assembly
	    specific_contig_stats <- contig_names_and_numbers %>% filter(`#Seq_id` == specific_gene_summary$`#Seq_id`)

	    if (nrow(specific_contig_stats != 0)) {
		    specific_contig_name=specific_contig_stats$`#Seq_id`[1]
		    specific_contig_number=specific_contig_stats$contig_number[1]
		    specific_contig_length=specific_contig_stats$`Length(bp)`[1]
	    } else {
		    specific_contig_name="NA"
                    specific_contig_number="NA"
                    specific_contig_length="NA"

		    print(paste("ERROR: contig statistics for", samples[s], methods[m], current_gene$gene, "are missing"))
            }
            
            if (read_depth_exists == TRUE) {
              
              # remove indicators of missing stops (>, <) from specific gene summary
              specific_gene_summary$Start <- str_replace_all(specific_gene_summary$Start, ">", "")
              specific_gene_summary$Start <- str_replace_all(specific_gene_summary$Start, "<", "")
              specific_gene_summary$End <- str_replace_all(specific_gene_summary$End, ">", "")
              specific_gene_summary$End <- str_replace_all(specific_gene_summary$End, "<", "")
              
              # sort contig names by length -> how ordered on read depths file
              #ordered_contigs  <- contig_summary %>% arrange(`Length(bp)`)
              # store unique contig names
              contig_names <- ordered_contigs$`#Seq_id`
              
              # store all unique read depth contig names
              read_depth_contig_names <- unique(read_depths_tibble$contig)
              
              # store start and end bp
              start_bp <- specific_gene_summary$Start
              end_bp <- specific_gene_summary$End
              
              # find contig name index in unique names list
              contig_name <- specific_gene_summary$`#Seq_id`
                
              contig_name_index <- which(contig_names == contig_name)
                
              # relate to name in read depths file
              read_depth_contig_name <- read_depth_contig_names[contig_name_index]
                
              # select read_depth rows matching contig name and bp range; store rows matching contig name and bp range
              read_depths_tibble_gene_match <- read_depths_tibble %>%
                filter(contig == read_depth_contig_name) %>% 
                filter(bp_1 %in% c(start_bp:end_bp))

if (nrow(read_depths_tibble_gene_match) > 0) {
              # calculate read_depth stats
              avg_read_depth <- mean(as.numeric(read_depths_tibble_gene_match$read_depth))
              med_read_depth <- median(as.numeric(read_depths_tibble_gene_match$read_depth))
              min_read_depth <- min(as.numeric(read_depths_tibble_gene_match$read_depth))
              max_read_depth <- max(as.numeric(read_depths_tibble_gene_match$read_depth))
} else {
      avg_read_depth <- "NA"
      med_read_depth <- "NA"
      min_read_depth <- "NA"
      max_read_depth <- "NA"
}
              
            }
            
          }
        
          # compile stats as comma delimited strings
          individual_gene_stats <- paste(samples[s], methods[m], current_gene$gene, num_Ns, num_bp, num_bp_minus_Ns, Ns_per_bp, bp_per_refbp, bp_minus_Ns_per_refbp, stop_codon_detected, num_gene_copies, num_internal_stops, avg_read_depth, med_read_depth, min_read_depth, max_read_depth, specific_contig_name, specific_contig_number, specific_contig_length, sep = ",")
        
          # check strings have expected number of commas
          # if yes, append strings to appropriate stats .csv file
          if (str_count(individual_gene_stats, ",") == 18) {
            write_lines(individual_gene_stats, paste0("./benchmark_stats/", OUTPREFIX, "gene_criteria.csv"), append = T)
          }
          # if longer, write NA for all entries and report error
          if (str_count(individual_gene_stats, ",") > 18) {
            individual_gene_stats <- paste(samples[s], methods[m], current_gene$gene, "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA", sep = ",")
            write_lines(individual_gene_stats, paste0("./benchmark_stats/", OUTPREFIX, "gene_criteria.csv"), append = T)
          
            print("ERROR: more stats compiled than expected")
          }
          # if shorter, write NA for all entries and report error
        
          if (str_count(individual_gene_stats, ",") < 18) {
            individual_gene_stats <- paste(samples[s], methods[m], current_gene$gene, "NA,NA,NA,NA,NA,NA,NA,NA,NA,NA", sep = ",")
            write_lines(individual_gene_stats, paste0("./benchmark_stats/", OUTPREFIX, "gene_criteria.csv"), append = T)
          
            print("ERROR: fewer stats compiled than expected")
          }
        
        }
        
      } 
    }
  }
}
