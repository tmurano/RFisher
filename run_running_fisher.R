run_running_fisher <- function(b1_file, b2_file, P1, P2, output_file = "summary_overlap.txt") { 
  # Load data
  b1 <- read.csv(b1_file, stringsAsFactors = FALSE)
  b2 <- read.csv(b2_file, stringsAsFactors = FALSE)
  
  # Set label from file name
  bioset1_label <- tools::file_path_sans_ext(b1_file)
  bioset2_label <- tools::file_path_sans_ext(b2_file)
  
  # Format
  read_and_format <- function(df) {
    df <- df[!duplicated(df$Gene), c("Gene", "Fold.Change")]
    colnames(df) <- c("Gene", "FC")
    df$Gene <- tolower(trimws(df$Gene))
    df
  }
  b1 <- read_and_format(b1)
  b2 <- read_and_format(b2)
  
  # Rank
  rank_fc <- function(df, dir) {
    df <- df[df$FC != 0, ]
    df <- if (dir == "up") df[df$FC > 0, ] else df[df$FC < 0, ]
    df[order(if (dir == "up") -df$FC else df$FC), ]
  }
  
  # Fisher test
  run_fisher_bidirectional <- function(b1_dir, b2_dir, P1, P2) {
    g1 <- rank_fc(b1, b1_dir)
    g2 <- rank_fc(b2, b2_dir)
    if (nrow(g1) == 0 || nrow(g2) == 0) return(list(M = 0, P = 1, Overlap = 0))
    K1 <- nrow(g2); K2 <- nrow(g1); ref1 <- g2$Gene; ref2 <- g1$Gene
    best_p1 <- 1; best_m1 <- 0
    for (N1 in 1:K2) {
      M1 <- length(intersect(g1$Gene[1:N1], ref1))
      pval1 <- phyper(M1 - 1, K1, P1 - K1, N1, lower.tail = FALSE)
      if (pval1 < best_p1) { best_p1 <- pval1; best_m1 <- M1 }
    }
    best_p2 <- 1; best_m2 <- 0
    for (N2 in 1:K1) {
      M2 <- length(intersect(g2$Gene[1:N2], ref2))
      pval2 <- phyper(M2 - 1, K2, P2 - K2, N2, lower.tail = FALSE)
      if (pval2 < best_p2) { best_p2 <- pval2; best_m2 <- M2 }
    }
    mean_logP <- (-log10(best_p1) + -log10(best_p2)) / 2
    list(M = round((best_m1 + best_m2) / 2),
         P = 10^(-mean_logP),
         Overlap = length(intersect(g1$Gene, g2$Gene)))
  }
  
  # Run comparisons
  directions <- list(
    upup     = c("up", "up", +1),
    updown   = c("up", "down", -1),
    downup   = c("down", "up", -1),
    downdown = c("down", "down", +1)
  )
  
  res <- lapply(names(directions), function(name) {
    d <- directions[[name]]
    r <- run_fisher_bidirectional(d[1], d[2], P1, P2)
    data.frame(Pair = name, OverlapGenes = r$Overlap, BestM = r$M,
               Pvalue = r$P, Score = -log10(r$P) * as.numeric(d[3]))
  })
  res_df <- do.call(rbind, res)
  
  S_total <- sum(res_df$Score)
  P_total <- exp(-abs(S_total))
  res_df <- rbind(res_df, data.frame(Pair = "Total", OverlapGenes = NA,
                                     BestM = NA, Pvalue = P_total, Score = S_total))
  
  # Output
  write.table(res_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(res_df)
}
