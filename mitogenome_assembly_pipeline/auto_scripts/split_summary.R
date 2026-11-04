library(tidyverse)
library(sjmisc)
library(janitor)

split_summary <- function(summary_file_path) {
  
  #read lines of text file
  summary <- read_lines(summary_file_path)
  
  #transform to tibble
  summary_t <- tibble(summary)
  
  if (nrow(summary_t) != 0) {
    #get table headers
    t_2_h <- NULL
    t_3_h <- NULL
    t_4_h <- NULL #stored as not all summaries have third table
    
    for (i in 1:nrow(summary_t)) {
      
      contains_pattern_2 <- str_contains(summary_t$summary[[i]], c("#Seq_id", "Start", "End", "Length(bp)", "Direction", "Type", "Gene_name", "Gene_prodcut", "Total_freq_occurred"), logic = "&")
      contains_pattern_3 <- str_contains(summary_t$summary[[i]], c("Protein coding genes totally found:"))
      contains_pattern_4 <- str_contains(summary_t$summary[[i]], c("#Gene", "total_missing_number"), logic = "&")
      
      if (contains_pattern_2 == TRUE) {
        t_2_h <- i
      }
      if (contains_pattern_3 == TRUE) {
        t_3_h <- i
      }
      if (contains_pattern_4 == TRUE) {
        t_4_h <- i
      }
      
    }
    
    #split tibble by header
    summary_A <- summary_t[row.names(summary_t) %in% c(1:(t_2_h-1)), ]
    summary_B <- summary_t[row.names(summary_t) %in% c(t_2_h:(t_3_h-1)), ]
    if (!is.null(t_4_h)) {
      summary_C <- summary_t[row.names(summary_t) %in% c(t_4_h:nrow(summary_t)), ]
    }
    
    #store table variable names
    if (!is.null(t_4_h)) {
      table_values <- c("summary_A", "summary_B", "summary_C")
    }
    if (is.null(t_4_h)) {
      table_values <- c("summary_A", "summary_B")
    }
    
    #create empty vector to store strings
    string_s <- c()
    
    #loop through tables
    for (i in 1:length(table_values)) {
      table_v <- table_values[i]
      
      #loop through lines in table
      for (n in 1:length(get(table_v)[[1]])) {
        if (n == 1) {
          #split header row
          string <- str_split_fixed(get(table_v)[[1]][[n]], pattern = " ", n = Inf)
        }
        if (n != 1) {
          #split observation rows; double space prevents splitting of multi-word terms
          string <- str_split_fixed(get(table_v)[[1]][[n]], pattern = "  ", n = Inf)
        }
        
        #replace non-values with NA
        string <- na_if(as.vector(string), c(""))
        string <- na_if(as.vector(string), c(" "))
        string <- na_if(as.vector(string), c("-----"))
        
        #remove NAs
        string <- string[!is.na(string)]
        
        #remove whitespace
        string <- str_squish(string)
        
        # check if string has expected length
        if (table_v == "summary_B") {
          
          #store gene names that contain spaces
          space_name <- c("NADH dehydrogenase subunit 1",
                          "NADH dehydrogenase subunit 2",
                          "NADH dehydrogenase subunit 3",
                          "NADH dehydrogenase subunit 4",
                          "NADH dehydrogenase subunit 4L",
                          "NADH dehydrogenase subunit 5",
                          "NADH dehydrogenase subunit 6",
                          "cytochrome c oxidase subunit I",
                          "cytochrome c oxidase subunit II",
                          "cytochrome c oxidase subunit III",
                          "ATP synthase F0 subunit 6",
                          "ATP synthase F0 subunit 8",
                          "cytochrome b",
                          "12S ribosomal RNA",
                          "16S ribosomal RNA")
          
          if (! length(string) %in% c(9,0)) {
            # find space name in string if exists
            contains_space_name <- FALSE
            check_space_name <- FALSE
            
            for (v in 1:length(string)) {
              check_space_name <- string[v] %in% space_name
              
              if (check_space_name == "TRUE") {
                contains_space_name <- TRUE
                space_name_location <- v
              }
              
            }
            
            #if string contains space name then split by single space individually
            if (contains_space_name == TRUE) {
              
              #create empty string storage
              string_ss <- c()
              
              for (v in 1:length(string)) {
                #if value doesn't contain space name then split
                if (v != space_name_location) {
                  string_v <- str_split_fixed(string[v], pattern = " ", n = Inf)
                  #bind split string value to stored string values
                  string_ss <- c(string_ss, string_v)
                }
                #if value does contain space name then add to string storage
                if (v == space_name_location) {
                  string_ss <- c(string_ss, string[space_name_location])
                }
              }
              #overwrite string with string_ss
              string <- string_ss
            }
            
            #if string contains no space name then split by single space
            if (contains_space_name == FALSE) {
              string <- str_split_fixed(get(table_v)[[1]][[n]], pattern = " ", n = Inf)
              string <- na_if(as.vector(string), c(""))
              string <- na_if(as.vector(string), c(" "))
              string <- na_if(as.vector(string), c("-----"))
              string <- string[!is.na(string)]
              string <- str_squish(string)
            }
            
          }
        }
        
        #bind split string vector to stored strings
        string_s <- as_tibble(rbind(string_s, string), .name_repair = "minimal")
        
      }
      
      #write stored tibble of strings to appropriate table variable
      if (i == 1) {
        summary_A_t <- string_s
        #replace column names with values in first row
        colnames(summary_A_t) <- as.character(summary_A_t[1,])
        #remove first row
        summary_A_t <- summary_A_t[-1,]
        #empty string storage
        string_s <- c()
      }
      if (i == 2) {
        summary_B_t <- string_s
        colnames(summary_B_t) <- as.character(summary_B_t[1,])
        #additionally remove second row as it is a divider between header and observations
        summary_B_t <- summary_B_t[-c(1,2),]
        string_s <- c()
      }
      if (i == 3) {
        summary_C_t <- string_s
        colnames(summary_C_t) <- as.character(summary_C_t[1,])
        summary_C_t <- summary_C_t[-c(1,2),]
        string_s <- c()
        
        #remove rows containing notes if present
        
        #create empty vector to store indices of to be removed rows
        remove <- c()
        
        for (r in 1:nrow(summary_C_t)) {
          
          contains_pattern <- str_contains(summary_C_t$total_missing_number[[r]], c("you got", "module of", "from the", "directory!"), logic = "or")
          
          if (contains_pattern == TRUE) {
            remove <- c(remove, r)
          }
          
        }
        
        if (! is.null(remove)) {
          #remove stored rows
          summary_C_t <- summary_C_t[-c(remove),]
        }
        
      }
      
    }
    
    if (!is.null(t_4_h)) {
      summary_total <- list(summary_A_t, summary_B_t, summary_C_t)
    }
    if (is.null(t_4_h)) {
      summary_total <- list(summary_A_t, summary_B_t)
    }
    
    #Add missing gene and detected gene summary tables
    
    #Create missing gene summary tibble
    
    if (length(summary_total) == 3){
      #change total_missing variable to numeric class
      summary_total[[3]] <- summary_total[[3]] %>% mutate(across(2, as.numeric))
      
      summary_total[[4]] <- summary_total[[3]] %>% 
        mutate(gene_type = case_when(str_detect(`#Gene`, "rRNA") ~ "rRNA",
                                     str_detect(`#Gene`, "tRNA") ~ "tRNA",
                                     .default = "CDS"
        ),
        .before = 2
        ) %>% 
        group_by(gene_type) %>% 
        summarise(sum(total_missing_number)) %>% 
        rename(total_missing_number = "sum(total_missing_number)")
      
      gene_types <- c("rRNA", "tRNA", "CDS")
      #identify gene types not in gene_type variable of summary_total[[4]]
      complete_genes <- gene_types[! gene_types %in% summary_total[[4]][[1]]]
      append_cg <- tibble(gene_type = complete_genes, total_missing_number = 0)
      
      summary_total[[4]] <- rbind(summary_total[[4]], append_cg)
    }
    
    if (length(summary_total) == 2){
      #add empty missing gene tibble
      summary_total[[3]] <- tibble(`#Gene` = NULL, total_missing_number = NULL)
      #add empty missing gene summary tibble
      summary_total[[4]] <- tibble(gene_type = c("rRNA", "tRNA", "CDS"), total_missing_number = c(0,0,0))
    }
    
    #Create detected gene summary tibble
    summary_total[[5]] <- summary_total[[4]] %>% 
      mutate(total_detected = case_when(gene_type == "rRNA" ~ 2 - total_missing_number,
                                        gene_type == "tRNA" ~ 22 - total_missing_number,
                                        gene_type == "CDS" ~ 13 - total_missing_number,
      ),
      .before = 2
      ) %>% 
      select(!total_missing_number)
    
  }
  
  #if summary file is empty, return empty variable; this should only occur when no CDS genes detected and using unmodified mitoz version
  if (!(nrow(summary_t) > 0)) {
    summary_total <- NULL
  }
  
  return(summary_total)
  
}




