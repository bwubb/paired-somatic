#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
acceptable_callers <- c("sequenza","cnvkit","ascat","purecn")

process_file <- function(file_path) {
  caller <- sub(".*/([^/]+)/[^/]+\\.hrd\\.txt$","\\1",file_path)

  # Check if the caller is acceptable
  if(!caller %in% acceptable_callers) {
    warning(paste("Skipping file:",file_path,"- Caller",caller,"not recognized."))
    return(NULL)
  }

  data <- read.csv(file_path,stringsAsFactors=FALSE)
  data <- data.frame(data[1],Caller=caller,data[-1])
  return(data)
}

data_list <- lapply(args,process_file)
data_list <- Filter(Negate(is.null),data_list)

combined_data <- do.call(rbind,data_list)
current_date <- format(Sys.Date(),"%Y%m%d")
output_file_name <- paste0("combined_HRDex_output_",current_date,".csv")

write.csv(combined_data,output_file_name,row.names=FALSE,quote=FALSE)
cat("Combined data written to:",output_file_name,"\n")
