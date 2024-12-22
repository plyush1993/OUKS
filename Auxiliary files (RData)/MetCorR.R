##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Installation                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(!"pacman" %in% rownames(installed.packages())){
  install.packages("pacman")
}
cat("Checking required packages (auto-installing if missing)\n")
pacman::p_load("crayon", "dplyr", "mgcv", "pbapply")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                            Correction Function                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MetCorR <- function(method = 2, int_data, order, class, batch, qc_label) {
  cat("\n")
  cat(crayon::red(
    c("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
  ),
  sep = "\n")
  
  cat(crayon::blue(
    c(
      "  .--.   ,--.         ,--.   ,-----.              ,------.    ",
      "  |   '.'   | .---. .-'  '-.'  .--.; .---. ,--.--.|  .--. '   ",
      "  |  |'.'|  || .-. :'-.  .-'|  |    | .-. ||  .--'|  '--'.'   ",
      "  |  |   |  ||  ---.  |  |  |  '--.;' '-' '|  |   |  |.  .    ",
      "  '--'   '--' `----'  '--'   '-----' '---' '--'   '--' '--'   ")
  ),
  sep = "\n")
  
  cat(crayon::red(
    c("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
  ),
  sep = "\n")
  cat("\n")

  if (method == 1) {
    cat(crayon::cyan("Method " %+% green$underline$bold ("1") %+% " has been selected\n"))
    cat(crayon::cyan("Used formula: y ~ s(order)\n"))
    }
  
  if (method == 2) {
    cat(crayon::cyan("Method " %+% green$underline$bold ("2") %+% " has been selected\n"))
    cat(crayon::cyan("Used formula: y ~ s(order, batch)\n"))
  }
  
  if (method != 1 & method != 2) {
    cat(crayon::cyan("Method " %+% green$underline$bold ("2") %+% " has been selected by default\n"))
    cat(crayon::cyan("Used formula: y ~ s(order, batch)\n"))
  }
  
  cat("\n")
  print("Fitting:")
  pboptions(style = 1, char = "-")
  if (method == 1) {
    order <- as.numeric(order)
    int_data_log <- int_data %>% mutate_all( ~log2(.+1.1))
    qc_id <-  grep(qc_label, class)
    int_data_log_ro <- as.data.frame(cbind(order, int_data_log))
    int_data_log_ro <- int_data_log_ro %>% mutate_all(as.numeric)
    int_data_log_ro_qc <- int_data_log_ro[qc_id,]
    fit_gam <- pblapply(2:ncol(int_data_log_ro_qc), function(t) gam(int_data_log_ro_qc[,t] ~ s(order), 
                                                                      data = int_data_log_ro_qc, method = "REML")) 
    
  print("Correction:")
    predict_gam <- pblapply(1:length(fit_gam), function(t) predict(fit_gam[[t]], data.frame(order)))
  } else { (method == 2)
    order <- as.numeric(order)
    batch <- as.numeric(batch)
    int_data_log <- int_data %>% mutate_all( ~log2(.+1.1))
    qc_id <-  grep(qc_label, class)
    int_data_log_ro_b <- as.data.frame(cbind(order, batch, int_data_log))
    int_data_log_ro_b <- int_data_log_ro_b %>% mutate_all(as.numeric)
    int_data_log_ro_b_qc <- int_data_log_ro_b[qc_id,]
    fit_gam <- pblapply(3:ncol(int_data_log_ro_b_qc), function(t) gam(int_data_log_ro_b_qc[,t] ~ s(order, batch), 
                                                                      data = int_data_log_ro_b_qc, method = "REML")) 
    
  print("Correction:")
    predict_gam <- pblapply(1:length(fit_gam), function(t) predict(fit_gam[[t]], data.frame(order, batch))) }
    
  cor_gam <- lapply(1:ncol(int_data_log), function(t) int_data_log[,t]-predict_gam[[t]]+mean(int_data_log[qc_id,t]))
  res_gam <- as.data.frame(t(do.call(rbind, cor_gam)))
  res_gam <- res_gam %>% mutate_all( ~2^(.))
  res_gam <- res_gam %>% mutate_all(as.numeric)
  rownames(res_gam) <- rownames(int_data_log)
  colnames(res_gam) <- colnames(int_data_log)
  return(res_gam)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~