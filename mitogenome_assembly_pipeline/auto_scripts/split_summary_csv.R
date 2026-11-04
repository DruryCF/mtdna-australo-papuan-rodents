#!/usr/bin/Rscript

source("auto_scripts/split_summary.R")

OUTPREFIX <- Sys.getenv("OUTPREFIX_VALUE_EXPORT")

samples <- read_lines(paste0("./", OUTPREFIX, "tmp/sampleIDs.txt"))
methods <- read_lines(paste0("./", OUTPREFIX, "tmp/benchmark_methods.txt"))

for (s in 1:length(samples)) {
  
  for (m in 1:length(methods)) {
    
    #store path to summary file
    path_summary_file <- paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".summary.txt")
    
    #check if summary file exists then split summary
    if (file.exists(path_summary_file)) {

      print(paste("splitting", path_summary_file))
      summary_tables <- split_summary(path_summary_file)

      #check that summary_tables is not equal to NULL, indicates summary file exists but is empty
      if (! is.null(summary_tables)) {
        write.csv(summary_tables[[1]], paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".contig_summary.csv"), row.names = F)
        write.csv(summary_tables[[2]], paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".gene_summary.csv"), row.names = F)
        write.csv(summary_tables[[3]], paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".missing_gene_summary.csv"), row.names = F)
        write.csv(summary_tables[[4]], paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".missing_gene_type_summary.csv"), row.names = F)
        write.csv(summary_tables[[5]], paste0("./", OUTPREFIX, "tmp/", samples[s], ".", methods[m], ".detected_gene_type_summary.csv"), row.names = F)
      } else {

        # if value is NULL then print and log error

        # get method name without _bm
	method <- str_remove(string = methods[m], pattern = "_bm")

	# write log message
	error_log <- paste0(samples[s], ",", method, "_mitoz_annotation:empty_summary.txt\n")

	# store paths to error logs
	failed_samples_csv <- paste0("./logs/", OUTPREFIX, "failed_samples_log.csv")
	failed_samples_txt <- paste0("./logs/", OUTPREFIX, "failed_samples.txt")

	# set just_created to FALSE
	just_created <- FALSE

	# check if error log .txt file exists, if not create it (.csv file created on script start)
	if (! file.exists(failed_samples_txt)) {

	  write("", failed_samples_txt)

	  # set just_created to TRUE
	  just_created <- TRUE
	
	}

	# append error information to log files
	cat(error_log, file = failed_samples_csv, append = T)
	cat(paste0(samples[s], "\n"), file = failed_samples_txt, append = T)

	# check if error log .txt file has just been created, if yes remove first blank line and set just_created to false
	
	if (just_created == TRUE) {

	  # read lines in file
	  empty_lines_failed_txt <- read_lines(failed_samples_txt)
	  # subset non_empty lines
	  no_empty_lines_failed_txt <- empty_lines_failed_txt[! empty_lines_failed_txt %in% ""]
	  # rewrite file with non_empty lines
	  write(no_empty_lines_failed_txt, failed_samples_txt)

	  # set just_created to FALSE
	  just_created <- FALSE

	}

	# print error message
	
	print(paste("ERROR:", samples[s], "assembled via", method, "has an empty summary.txt file"))

      }

    }
    
  }
  
}