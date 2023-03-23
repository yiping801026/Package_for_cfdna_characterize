#' tss_top_genes
#'
#' NO2.STEP
#' @param orig_file_path original file path for all frag.csv files
#' @param N_health_samples number of health samples
#' @param N_cancer_samples number of cancer samples
#'
#' @return df_ratio
#' @export
#' @importFrom utils read.csv
#' @examples
frag_ratio <- function(orig_file_path,N_health_samples,N_cancer_samples){
  ### read in the csv files
  file_names_path <- dir(orig_file_path, pattern = "*.csv", recursive = F, full.names = T)

  ### DataFrame to store all ratio
  df_ratio_all <- read.csv(file_names_path[1],sep="")

  for (i in 2:length(file_names_path)) {
    df_ratio_all <- cbind(df_ratio_all, read.csv(file_names_path[i],sep=""))}

  ### change the file_names
  file_names <- dir(orig_file_path, pattern = "*.csv", recursive = F, full.names = F)

  file_names_less <- c(1:length(file_names))
  for (i in 1:length(file_names)){
    m <- strsplit(file_names[i],'[.]')[[1]][1]
    file_names_less[i] <- strsplit(m,'_')[[1]][1]}

  ### remove 0s of short && long
  df_ratio_all <- df_ratio_all[which(df_ratio_all$short_sum != 0),]

  ### get the ratios
  df_ratio <- df_ratio_all[,seq(from=3, to=ncol(df_ratio_all), by=3)]
  colnames(df_ratio) <- file_names_less

  ### without pos
  df_ratio <- data.frame(df_ratio)
  rownames(df_ratio) <- 1:nrow(df_ratio)

  ### scale for individuals
  df_ratio[nrow(df_ratio)+1,] <-  apply(df_ratio,2,mean)

  ### subtract the mean of individual
  nsample = nrow(df_ratio)-1
  for(i in 1:ncol(df_ratio)){
    df_ratio[1:nsample,i] <-  (df_ratio[1:nsample,i] - df_ratio[nsample+1,i])
  }

  ### remove the mean
  df_ratio <- df_ratio[-nrow(df_ratio),]

  return(df_ratio)
}
