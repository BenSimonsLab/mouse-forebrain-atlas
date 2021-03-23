plot_UMAP_expr <- function(geneID){
  
  tbl_expr = get_adata_raw_expr(geneID)
  
  metadata = metadata %>%
              left_join(tbl_expr, by="cellID") %>%
              mutate(expr_Gene = get(paste0("expr_", geneID)))
              
  # expr_Gene
  plt_UMAP_expr = ggplot(slice(metadata, sample(1:n())), aes(x=UMAP_dim1, y=UMAP_dim2, colour=expr_Gene)) +
                              # stat_summary_hex(bins=100) +
                              # geom_point(size=0.2, alpha=0.8, shape=16) +
                              geom_point_rast(size=0.1, stroke=0, alpha=0.8, shape=16, raster.dpi=600) +
                              theme_publication +
                              scale_colour_viridis(option="magma") +
                              labs(x="UMAP 1", y="UMAP 2") +
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank()) +
                              labs(colour=paste0(geneID, " [log1p]")) +
                              guides(colour = guide_colorbar(barwidth = 0.3,
                                                           barheight = 2)) #+
                              # theme(legend.position="top")


  return(plt_UMAP_expr)
}



