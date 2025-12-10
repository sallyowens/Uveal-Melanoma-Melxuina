#!/usr/bin/env Rscript

# FACETS Analysis for WGS Data

library(facets)

# Set paths using HOME environment variable
home_dir <- Sys.getenv("HOME")
project_dir <- file.path(home_dir, "Project")

pileup_dir <- file.path(project_dir, "FACETS/pileups")
results_dir <- file.path(project_dir, "FACETS/results")
plots_dir <- file.path(project_dir, "FACETS/plots")
reports_dir <- file.path(project_dir, "FACETS/reports")

# Create output directories
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== FACETS WGS Analysis ===\n")
cat("Pileup directory:", pileup_dir, "\n")
cat("Results directory:", results_dir, "\n\n")

# Process function
process_sample <- function(pileup_file) {
    patient_id <- gsub("\\.sorted\\.pileup$", "", basename(pileup_file))
    cat("\n==========================================\n")
    cat("Processing:", patient_id, "\n")
    cat("==========================================\n")
    
    tryCatch({
        # Read data
        cat("Reading SNP matrix...\n")
        rcmat <- readSnpMatrix(pileup_file)
        
        if (nrow(rcmat) == 0) {
            cat("ERROR: No data in pileup file\n")
            return(NULL)
        }
        
        cat("SNPs loaded:", nrow(rcmat), "\n")
        
        # Preprocess (WGS parameters)
        cat("Preprocessing sample...\n")
        preproc <- preProcSample(
            rcmat,
            gbuild = "hg38",
            cval = 50,           # More stringent for WGS
            ndepth = 25,         # Higher minimum depth
            het.thresh = 0.25
        )
        
        cat("Segments after preprocessing:", nrow(preproc$out), "\n")
        
        # Process
        cat("Processing with FACETS algorithm...\n")
        out <- procSample(preproc, min.nhet = 15)
        
        # Fit model
        cat("Fitting copy number model...\n")
        fit <- emcncf(out, min.nhet = 15, maxiter = 30)
        
        if (is.null(fit$cncf)) {
            cat("ERROR: FACETS fitting failed\n")
            return(NULL)
        }
        
        # Save results
        cat("Saving results...\n")
        saveRDS(fit, file.path(results_dir, paste0("fit_", patient_id, ".rds")))
        saveRDS(out, file.path(results_dir, paste0("out_", patient_id, ".rds")))
        
        # Save segments
        cncf_df <- fit$cncf
        cncf_df$sample <- patient_id
        cncf_df$purity <- fit$purity
        cncf_df$ploidy <- fit$ploidy
        
        write.csv(cncf_df, 
                  file.path(results_dir, paste0("segments_", patient_id, ".csv")),
                  row.names = FALSE)
        
        # Plot
        cat("Generating plots...\n")
        png(file.path(plots_dir, paste0("facets_", patient_id, ".png")),
            width = 1200, height = 900, res = 150)
        plotSample(out, emfit = fit, sname = patient_id)
        dev.off()
        
        # Spider plot
        if (!is.null(fit$dipLogR)) {
            png(file.path(plots_dir, paste0("spider_", patient_id, ".png")),
                width = 800, height = 600, res = 150)
            logRlogORspider(out$out, out$dipLogR)
            title(paste("Spider Plot -", patient_id))
            dev.off()
        }
        
        # Summary
        cat("\n=== Summary for", patient_id, "===\n")
        cat("Purity:", round(fit$purity, 3), "\n")
        cat("Ploidy:", round(fit$ploidy, 3), "\n")
        cat("Diploid LogR:", round(fit$dipLogR, 3), "\n")
        cat("Segments:", nrow(cncf_df), "\n")
        cat("===================================\n")
        
        # Return summary
        data.frame(
            sample = patient_id,
            purity = fit$purity,
            ploidy = fit$ploidy,
            dipLogR = fit$dipLogR,
            n_segments = nrow(cncf_df),
            stringsAsFactors = FALSE
        )
        
    }, error = function(e) {
        cat("ERROR in", patient_id, ":", e$message, "\n")
        return(NULL)
    })
}

# Find pileup files
setwd(pileup_dir)
pileup_files <- list.files(pattern = ".*\\.sorted\\.pileup$", full.names = TRUE)

if (length(pileup_files) == 0) {
    cat("ERROR: No pileup files found in:", pileup_dir, "\n")
    quit(status = 1)
}

cat("Found", length(pileup_files), "pileup files:\n")
for (file in pileup_files) {
    cat(" -", basename(file), "\n")
}
cat("\n")

# Process all samples
all_results <- lapply(pileup_files, process_sample)
all_results <- do.call(rbind, Filter(Negate(is.null), all_results))

# Save combined summary
if (!is.null(all_results) && nrow(all_results) > 0) {
    write.csv(all_results,
              file.path(reports_dir, "facets_summary_all_samples.csv"),
              row.names = FALSE)
    
    cat("\n==========================================\n")
    cat("FINAL SUMMARY - ALL SAMPLES\n")
    cat("==========================================\n")
    print(all_results)
    
    cat("\nSummary Statistics:\n")
    cat("Samples processed:", nrow(all_results), "\n")
    cat("Purity range:", round(min(all_results$purity), 3), "-", 
        round(max(all_results$purity), 3), "\n")
    cat("Ploidy range:", round(min(all_results$ploidy), 3), "-", 
        round(max(all_results$ploidy), 3), "\n")
    cat("Mean purity:", round(mean(all_results$purity), 3), "\n")
    cat("Mean ploidy:", round(mean(all_results$ploidy), 3), "\n")
} else {
    cat("\nWARNING: No samples were successfully processed\n")
}

cat("\n=== FACETS Analysis Complete! ===\n")
cat("Results saved in:\n")
cat(" - Results:", results_dir, "\n")
cat(" - Plots:", plots_dir, "\n")
cat(" - Reports:", reports_dir, "\n")
