#' frag_plot
#'
#' NO3.STEP
#' @param frag_matrix the all matrix for the plot including control && test samples
#' @param N_health_samples number of samples of control type
#' @param N_cancer_samples number of samples of health type
#' @param type type of samples selected 'control' or 'test' (for the color of plot)
#'
#' @return df_ratio
#' @export
#' @examples


frag_plot <- function(frag_matrix,N_health_samples,N_cancer_samples,type){

  if(type == 'control'){
    start_pos = 1
    end_pos = N_health_samples
    line_color = 'navy'
    plot_title = "Fragment ratio pattern of control samples"}
  else {
    start_pos = 1 + N_health_samples
    end_pos = N_health_samples + N_cancer_samples
    line_color = 'brown3'
    plot_title = "Fragment ratio pattern of test samples"}

  ### plot the samples
  df_ratio_plot2 <- data.frame(frag_matrix[,start_pos:end_pos])

  df_ratio_plot2$pos <- 1:nrow(df_ratio_plot2)
  df_ratio_plot2 <- data.frame(reshape2::melt(df_ratio_plot2,id="pos"))
  colnames(df_ratio_plot2)[2:3] <- c('sample','ratio')

  df_ratio_plot2$pos <- as.integer(df_ratio_plot2$pos)
  df_ratio_plot2$sample <- as.character(df_ratio_plot2$sample)

  p2 <- ggplot2::ggplot(data = df_ratio_plot2,ggplot2::aes(x=df_ratio_plot2$pos,y=df_ratio_plot2$ratio,group = sample,color=sample))+
    ggplot2::ggtitle(plot_title)+
    ggplot2::xlab('Bin positions')+
    ggplot2::ylab('Normalized fragment ratio')+
    ggplot2::ylim(-0.025,0.025)+
    ggplot2::geom_line(color=line_color)+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major=ggplot2::element_line(colour=NA),
                   panel.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                   plot.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = c(.075,.915),
                   legend.box.background = ggplot2::element_rect(color="black"))

  pdf(file=paste0(plot_title,'.pdf'),width=20,height=4,bg = "white")
  print(p2)
  dev.off()
}
