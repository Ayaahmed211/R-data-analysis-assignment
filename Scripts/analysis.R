
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1
  if (logFC > 1 && padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 && padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}


files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")
dir.create("Results", showWarnings = FALSE)

deg_results <- list()

for (f in files_to_process) {
  dat <- read.csv(file.path("Raw_Data", f), header = TRUE, stringsAsFactors = FALSE)
  dat$padj[is.na(dat$padj)] <- 1
  dat$status <- mapply(classify_gene, dat$logFC, dat$padj)
  outname <- file.path("Results", paste0("Processed_", f))
  write.csv(dat, outname, row.names = FALSE)
  print(table(dat$status))
  deg_results[[f]] <- dat
}


all_dat <- do.call(rbind, deg_results)
summary_df <- as.data.frame(table(all_dat$status))
write.csv(summary_df, file = file.path("Results", "summary_counts.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), con = file.path("Results", "sessionInfo.txt"))
