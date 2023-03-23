#' cnv_plot_control
#'
#' NO4.STEP
#' @param X_healthy_100KCNV matrix to plot cnv
#'
#' @return he_plot
#' @export
#' @importFrom stats median
#' @examples


cnv_plot_control <- function(X_healthy_100KCNV){
  #X_healthy_100KCNV = cnv_control_matrix
  ### pre-treat the matrix (remove the position name & add colnames)
  X_healthy_100KCNV <- X_healthy_100KCNV[,-2]
  colnames(X_healthy_100KCNV)[1] <- 'chr'
  colnames(X_healthy_100KCNV)[2:ncol(X_healthy_100KCNV)] <- 1:(ncol(X_healthy_100KCNV)-1)
  raw_save <- X_healthy_100KCNV
  he_raw <- data.frame(raw_save[,2:ncol(X_healthy_100KCNV)])

  ### remove rows including '0'
  he_zero_pos <- which(he_raw==0)
  he_zero_row <- as.integer(unique(sort(he_zero_pos%%nrow(he_raw))) )
  write.table(he_zero_row,'he_zero_row.txt',col.names = F)
  he_raw <- he_raw[-he_zero_row,]

  ### normalize: calculate sum of patient && percentage to replace every cell
  t_he_raw <- data.frame(t(he_raw))
  t_he_raw$sum <- apply(t_he_raw,1,sum)
  sum_col <- ncol(t_he_raw)
  t_he_raw_out <- t_he_raw[,-sum_col]/t_he_raw[,sum_col]
  he_raw_out <- data.frame(t(t_he_raw_out))

  he_raw_out$median <- apply(he_raw_out,1,median)####should be saved later

  ###  delete by median
  per_he_100CNV_wide <- data.frame(he_raw_out[,-ncol(he_raw_out)]/he_raw_out$median)

  per_he_100CNV_wide$pos <- 1:nrow(per_he_100CNV_wide)
  ################################################################
  ### log2-------healthy did log2 or not?
  ### without log2 range :  7.137759e-04 4.562386e+03
  ### with log2 range : -10.45224  12.15557   ???????
  per_he_100CNV_wide[,-ncol(per_he_100CNV_wide)] <- log2(per_he_100CNV_wide[,-ncol(per_he_100CNV_wide)])
  ################################################################
  ### wide to long
  per_he_100CNV_wide$pos <- factor(per_he_100CNV_wide$pos,levels = 1:nrow(per_he_100CNV_wide))
  #summary(per_he_100CNV_wide)
  N_samples = ncol(per_he_100CNV_wide)-1
  per_he_100CNV_long <- tidyr::gather(per_he_100CNV_wide,sample_ID,CNV_ratio,1:N_samples,factor_key=TRUE)
  range(per_he_100CNV_long$CNV_ratio)
  ### change position back to integer
  per_he_100CNV_long$pos <- as.integer(per_he_100CNV_long$pos )
  ### only use pos && CNV_ratio
  he_plot <- per_he_100CNV_long[,c(1,3)]

  ### plot 1
  p3 <- ggplot2::ggplot(data = he_plot , ggplot2::aes(he_plot$pos, he_plot$CNV_ratio)) +
    ggplot2::geom_hex(bins = 150) +
    ggplot2::scale_fill_continuous(type = "viridis") +
    ggplot2::theme_bw()

  pdf(file='cnv1.pdf',bg = "white")
  print(p3)
  dev.off()

  ### plot 2

  p_he <-
    ggplot2::ggplot(he_plot, ggplot2::aes(he_plot$pos, he_plot$CNV_ratio) ) +
    #xlim(0,27633)+
    ggplot2::ylim(-0.5,0.5)+ggplot2::ylab('log2 CNV ratio')+ ggplot2::xlab('bin')+
    #scale_fill_continuous(type = "viridis") +
    scico::scale_fill_scico(palette = "lajolla")+
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), geom = "polygon"
                    ,colour="white")

  p4 <- p_he+ggplot2::ggtitle("Healthy samples") + ggplot2::theme_bw() + ggplot2::theme(panel.grid=ggplot2::element_blank())

  pdf(file='cnv2.pdf',bg = "white")
  print(p4)
  dev.off()

  return(he_raw_out)
}

