#' tss_top_genes
#'
#' NO1.STEP
#' @param test patient matrix
#' @param control healthy matrix
#' @param ntop_gene the gene numbers used in the heatmap
#'
#'
#' @return table_cal_all
#' @export
#'
#' @importFrom stats predict
#' @importFrom caret preProcess
#' @importFrom stats t.test
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @examples
tss_top_genes <- function(test,control,ntop_gene){


  all_tss_raw <- dplyr::inner_join(test,control,by=c('gene','transcript'))
  ### change -1 into 0
  all_tss_raw[all_tss_raw == -1] = 0
  ### remove all 0 regions
  all_tss_raw$sum <- apply(all_tss_raw[,3:ncol(all_tss_raw)],1,sum)
  all_tss_rezero <- data.frame(all_tss_raw[which(all_tss_raw$sum != 0),])

  ### remove sum
  all_tss_rezero <- all_tss_rezero[,-which(colnames(all_tss_rezero)=='sum')]

  ### get the number of rows && columns
  N_gene <- nrow(all_tss_rezero)
  N_column <- ncol(all_tss_rezero)
  N_health_samples <- ncol(control)-2
  N_test_samples <- ncol(test)-2
  new_start = N_test_samples+1
  new_end = N_test_samples+N_health_samples
  ### change the matrix from vertical into horizontal
  table_cal_all_r <-  data.frame(t(all_tss_rezero[,3:N_column]))

  ### scale && center
  preProcValues <- caret::preProcess(table_cal_all_r, method = c("center", "scale"))
  table_cal_all_nor <- predict(preProcValues, table_cal_all_r)
  table_cal_all <- data.frame(t(table_cal_all_nor))

  ### calculate t-test of 2 group(cancer vs healthy) for every row (gene)
  tout <- rep(1,N_gene)

  for(i in 1:N_gene){tout[i]<-t.test(table_cal_all[i,1:N_test_samples], table_cal_all[i,new_start:new_end])$p.value}
  tout_withID<-data.frame(cbind(as.numeric(tout),all_tss_rezero$gene))
  tout_withID[,1] <- as.double(tout_withID[,1])
  tout_sort<-data.frame(tout_withID[order(as.numeric(tout_withID[,1])),])

  ### genes used for further heatmap plot
  top_selected_genes <- tout_sort[c(1:as.numeric(ntop_gene)),]
  rowtop_selected_genes <- rownames(top_selected_genes)

  ### if just use the raw_matrix to plot heatmap
  raw_plot_table <- all_tss_rezero[rowtop_selected_genes,3:N_column]

  ### plot start
  raw_plot_data <- as.matrix(raw_plot_table)

  #### anno
  annotation_col <- data.frame(
    group = c(rep('cancer',N_test_samples),rep('healthy',N_health_samples))
  )
  row.names(annotation_col) <- colnames(raw_plot_data)

  annotation_colors = list(
    group = c(cancer='dark red',healthy='navy')
  )

  ### plot to store the cluster inform of this plot
  #result=pheatmap::pheatmap(raw_plot_data,
  #                legend_breaks = 0:1,
  #                cellwidth = 5, cellheight = 9,
  #                scale = 'none',
  #                show_colnames = F,show_rownames = F,
  #                color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  #               annotation_col=annotation_col
  #)

  ### use the cluster inform
  #order_row = result$tree_row$order
  #data.2 = data.frame(raw_plot_data[order_row,])

  ### plot
  p1 <- pheatmap::pheatmap(raw_plot_data,
                            cluster_rows = F,
                            cluster_cols = F,
                            cellwidth = 5, cellheight = 9,
                            scale = 'none',
                            show_colnames = F,show_rownames = F,
                            color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                            annotation_col=annotation_col,
                            annotation_colors=annotation_colors,
                            border_color=NA )
  pdf(file=paste0('heatmap.pdf'),bg = "white")
  print(p1)
  dev.off()


  #return(raw_plot_data)

}














