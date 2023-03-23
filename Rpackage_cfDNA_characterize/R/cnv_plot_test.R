#' cnv_plot_test
#'
#' NO4.STEP
#' @param cnv_test_matrix matrix to plot cnv
#' @param he_raw_out to store the control median
#' @param he_zero_row to store the rows zero exists for control&&test groups (should be a list)
#'
#' @return he_plot
#' @export
#' @importFrom stats median
#' @examples

cnv_plot_test <- function(cnv_test_matrix,he_raw_out,he_zero_row){
  #cnv_test_matrix = cnv_test_matrix
  #he_raw_out = he_raw_out
  #he_zero_row = he_zero_row$V2
  ### pre-treat the matrix (remove the position name & add colnames)
  cnv_test_matrix <- cnv_test_matrix[,-2]
  colnames(cnv_test_matrix)[1] <- 'chr'
  raw_save_patient <- cnv_test_matrix

  ### for patients
  #pa_raw <- data.frame(raw_save_patient[,2:8])
  pa_raw <- data.frame(raw_save_patient[,2:ncol(raw_save_patient)])
  #colnames(pa_raw)
  ### remove rows including '0' from healthy
  pa_zero_pos <- which(pa_raw==0)
  pa_zero_row <- as.integer(unique(sort(pa_zero_pos%%nrow(pa_raw))) )
  #setdiff(pa_zero_row,he_zero_row)
  #setdiff(he_zero_row,pa_zero_row)
  pa_raw <- pa_raw[-he_zero_row,]
  # raw_save_try[,2:36] <- raw_save_try[,2:36] + 0.01

  ### normalize: calculate sum of patient && percentage to replace every cell
  t_pa_raw <- data.frame(t(pa_raw))
  t_pa_raw$sum <- apply(t_pa_raw,1,sum)
  pasum_col <- ncol(t_pa_raw)
  t_pa_raw_out <- t_pa_raw[,-pasum_col]/t_pa_raw[,pasum_col]
  pa_raw_out <- data.frame(t(t_pa_raw_out))

  #pa_raw_out$median <- apply(pa_raw_out,1,median)####should be saved later

  ###  delete by median
  per_pa_100CNV_wide <- data.frame(pa_raw_out/he_raw_out$median)

  per_pa_100CNV_wide$pos <- 1:nrow(per_pa_100CNV_wide)
  ################################################################
  ### log2-------healthy did log2 or not?
  ### without log2 range :  7.137759e-04 4.562386e+03
  ### with log2 range : -10.45224  12.15557   ???????
  per_pa_100CNV_wide[,-ncol(per_pa_100CNV_wide)] <- log2(per_pa_100CNV_wide[,-ncol(per_pa_100CNV_wide)])
  ################################################################
  ### wide to long
  per_pa_100CNV_wide$pos <- factor(per_pa_100CNV_wide$pos,levels = 1:nrow(per_pa_100CNV_wide))
  #summary(per_pa_100CNV_wide)
  N_samples_pa = ncol(per_pa_100CNV_wide)-1
  #per_he_100CNV_long <- tidyr::gather(per_he_100CNV_wide,sample_ID,CNV_ratio,1:N_samples,factor_key=TRUE)
  per_pa_100CNV_long <- tidyr::gather(per_pa_100CNV_wide,sample_ID,CNV_ratio,1:N_samples_pa,factor_key=TRUE)

  #per_pa_100CNV_long <- gather(per_pa_100CNV_wide,sample_ID,CNV_ratio,1:19,factor_key=TRUE)
  #range(per_pa_100CNV_long$CNV_ratio)
  ### change position back to integer
  per_pa_100CNV_long$pos <- as.integer(per_pa_100CNV_long$pos )
  ### only use pos && CNV_ratio
  pa_plot <- per_pa_100CNV_long[,c(1,3)]

  ### plot 2

  p_he <-
    ggplot2::ggplot(pa_plot, ggplot2::aes(pa_plot$pos, pa_plot$CNV_ratio) ) +
    #xlim(0,27633)+
    ggplot2::ylim(-0.5,0.5)+ggplot2::ylab('log2 CNV ratio')+ ggplot2::xlab('bin')+
    #scale_fill_continuous(type = "viridis") +
    scico::scale_fill_scico(palette = "lajolla")+
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), geom = "polygon"
                             ,colour="white")

  p5 <- p_he+ggplot2::ggtitle("cancer samples") + ggplot2::theme_bw() + ggplot2::theme(panel.grid=ggplot2::element_blank())

  pdf(file='cnv3.pdf',bg = "white")
  print(p5)
  dev.off()


}



