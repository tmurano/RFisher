rm(list = ls())

wd <- selectDirectory("Select working directory", path = getwd())
setwd(wd)

# ==== Load data ====

get_bioset <- function(label) {
  file <- selectFile(paste0("📂 Select ", label), path = getwd())
  if (is.null(file)) stop("❌ Cancelled.")
  p <- if (showQuestion(label, "Mouse?", "Mouse", "Others")) 20195 else
    if (showQuestion(label, "Human?", "Human", "Others")) 21196 else {
      val <- as.integer(showPrompt(label, "Enter gene count:"))
      if (is.na(val) || val <= 0) stop("❌ Invalid P.") else val
    }
  list(file = file, P = p)
}

b1 <- get_bioset("Bioset1")
b2 <- get_bioset("Bioset2")

b1_file <- basename(b1$file)
b2_file <- basename(b2$file)


# ==== Function  ====

source("run_running_fisher.R")
res_df <- run_running_fisher(b1_file, b2_file, P1 = b1$P, P2 = b2$P)  

source("generate_plot.R")
generate_plot(b1_file, b2_file, res_df)


# ==== end of program  ====


