library(reshape2)
library(readr)
library(dplyr)

library(ggthemes)
library(ggrepel)

get_adata_raw_expr <- function(vec_geneIDs){
  # extracts raw count data from anndata object
  # input: vector of gene names
  # output: tibble containing cellIDs and gene expression
  assign("vec_geneIDs", vec_geneIDs, envir=.GlobalEnv)
  
  # raise error when geneID are missing/not detected
  py_run_string("array_geneIDs = np.atleast_1d(r.vec_geneIDs)")
  py_run_string("geneIDs_missing = [x for x in array_geneIDs if x not in adata.raw.var_names]")
  if (length(py$geneIDs_missing) > 0){
    stop(paste0('geneIDs missing: ', py$geneIDs_missing, ' '))
  }
  
  py_run_string("filter_genes = [x in array_geneIDs for x in adata.raw.var_names]")
  py_run_string("cmat = adata.raw.X[:, np.nonzero(filter_genes)[0]].todense()")
  py_run_string("cellIDs_cmat = adata.raw.obs_names.values")
  py_run_string("geneIDs_cmat = adata.raw.var_names[filter_genes].values")
  
  cmat_expr = py$cmat
  colnames(cmat_expr) = paste0("expr_", py$geneIDs_cmat)
  tbl_expr = as_tibble(cmat_expr) %>%
                add_column(cellID = py$cellIDs_cmat, .before=1)
  
  return(tbl_expr)
}

plot_UMAP_leiden <- function(leiden_key){
  
  metadata = metadata %>%
              mutate(leiden = get(leiden_key))
              
  cl_cent = metadata %>%
                  select(UMAP_dim1, UMAP_dim2, leiden) %>%
                  group_by(leiden) %>%
                  summarize_all(mean)
  
  # randomise order for plot
  plt_UMAP_leiden = ggplot(slice(metadata, sample(1:n())), aes(x=UMAP_dim1, y=UMAP_dim2, colour=leiden)) +
                              geom_point(size=0.5, alpha=0.8, shape=16) +
                              geom_label_repel(aes(x=UMAP_dim1, y=UMAP_dim2, label=leiden, fill=leiden),
                                      label.size=NA, alpha=0.6, label.padding = 0.15, point.padding = NA,
                                      data=cl_cent, size=2.5, show.legend=FALSE, inherit.aes=F) +
                              theme_publication +
                              labs(x="UMAP 1", y="UMAP 2") +
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank()) +
                              # guides(color = guide_legend(override.aes = list(size=3))) +
                              # labs(colour="Leiden")
                              labs(title=paste0(leiden_key)) +
                              theme(legend.position="none")

  return(plt_UMAP_leiden)
}


get_marker_genes <- function(fileID, leiden_res, filter_log2fc=-0.25, filter_frac = 0.25, filter_qval=0.05, n_DE = 30){
  file_paths = list.files(paste0("output/", fileID, "/DE-diffxpy"), patter=paste0(leiden_res, ".+_rest.csv.gz"), full.names=T)
  
  tbl_results = tibble(ranking = 1:n_DE)
  
  for (path in file_paths){
    tbl_DE_cl = read_csv(path)
    tbl_DE_cl = tbl_DE_cl %>%
          arrange(qval) %>%
          mutate(frac_max = pmax(frac_1, frac_2)) %>%
          filter(log2fc < filter_log2fc) %>%
          filter(frac_max >= filter_frac) %>%
          filter(qval < filter_qval) %>%
          arrange(qval, log2fc)
    geneIDs_cl = tbl_DE_cl %>% slice(1:n_DE) %>% pull(geneID)
    if (length(geneIDs_cl) < n_DE){
      geneIDs_cl = c(geneIDs_cl, rep(NA, n_DE - length(geneIDs_cl)))
    }
    clusterID = gsub("_vs_rest.csv.gz", "", gsub(".*_cl_", "", path))
    tbl_results = tbl_results %>% mutate(!!clusterID := geneIDs_cl)
  }
  tbl_results = tbl_results %>% select(-ranking)
  
  
  write_csv(tbl_results, paste0("output/", fileID, "/DE-diffxpy/DE_wilcoxon_", fileID, "_DE_", leiden_res, "_top100.csv.gz"))
  return()
}


plot_DE_heatmap <- function(fileID, leiden_res){
  
  de_wilcoxon = read_csv(paste0("output/", fileID, "/DE-diffxpy/DE_wilcoxon_", fileID, "_DE_", leiden_res, "_top100.csv.gz")) %>%
                                  slice(1:20)
  
  DE_genes = de_wilcoxon %>%
                tidyr::pivot_longer(everything(), names_to="clusterID", values_to="geneID") %>%
                arrange(desc(row_number())) %>%
                arrange(clusterID) %>%
                pull(geneID) %>%
                unique()  
  
  tbl_expr = get_adata_raw_expr(DE_genes)
  
                              
  tbl_DE_expr = metadata %>%
                  select(cellID, !!(leiden_res_annot), !!(leiden_res)) %>%
                  mutate(clusterID = get(leiden_res)) %>%
                  left_join(tbl_expr, by="cellID") %>%
                  tidyr::pivot_longer(starts_with("expr_"), names_to="geneID", values_to="expr") %>%
                  mutate(geneID = gsub("expr_", "", geneID)) %>%
                  group_by(geneID) %>%
                  mutate(expr_scaled = (expr - min(expr))/(max(expr) - min(expr))) %>%
                  mutate(expr_zscore = scale(expr)) %>%
                  mutate(expr_scaled_zscore = scale(expr_scaled)) %>%
                  mutate(expr_zscore_cut = ifelse(expr_zscore > 2.5, 2.5, ifelse(expr_zscore < -2.5, -2.5, expr_zscore))) %>%
                  ungroup() %>%
                  arrange(clusterID) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID))) %>%
                  mutate(geneID = factor(geneID, levels = DE_genes))
                  

  plt_DE_heatmap = ggplot(tbl_DE_expr, aes(x=cellID, y=geneID, fill=expr_zscore_cut)) +
                              geom_tile() +
                              scale_fill_distiller(palette = 'RdYlBu') +
                              # scale_fill_gradient2(low = "magenta", mid = "black", high = "yellow") +
                              theme_publication +
                              labs(x=NULL, y=NULL) +
                              scale_x_discrete(expand = c(0,0)) +
                              scale_y_discrete(expand = c(0,0)) +
                              labs(fill="Expr [cut z-score log1p CPM]") +
                              # theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank()) +
                              guides(colour = guide_colorbar(barwidth = 0.3,
                                                             barheight = 2))



  tbl_cluster = metadata %>%
                  select(cellID, !!(leiden_res)) %>%
                  mutate(leiden_res = get(leiden_res)) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID)))

  tbl_cluster_frac = tbl_cluster %>%
                       group_by(leiden_res) %>%
                       summarise(n_cluster = n()) %>%
                       mutate(frac_cluster = n_cluster/sum(n_cluster)) %>%
                       arrange(leiden_res) %>%
                       mutate(x_max = cumsum(frac_cluster))

  tbl_cluster_frac = tbl_cluster_frac %>%
                     mutate(x_min = c(0, head(tbl_cluster_frac$x_max, -1))) %>%
                     mutate(x_mid = x_min + (x_max - x_min)/2)

   plt_cluster = ggplot(tbl_cluster_frac, aes(fill=leiden_res)) +
                       geom_rect(aes(xmin = x_min, xmax = x_max), ymin=-Inf, ymax=Inf) +
                       theme_publication +
                       # scale_x_discrete(expand = c(0,0)) +
                       # theme(axis.text.x=element_blank(),
                       #       axis.ticks.x=element_blank(),
                       #       axis.text.y=element_blank(),
                       #       axis.ticks.y=element_blank()) +
                       theme(legend.position="none") +
                       # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                       scale_x_continuous(name=NULL, breaks=tbl_cluster_frac$x_mid, labels=(tbl_cluster_frac$leiden_res), expand = c(0,0))
                       # , guide = guide_axis(angle=45)) +
                       labs(x=NULL)


  plt_panel = (plt_DE_heatmap / plt_cluster) + plot_layout(heights = c(1, 0.02))

return(plt_panel)
}


plot_DE_dotplot <- function(fileID, leiden_res){
  
  de_wilcoxon = read_csv(paste0("output/", fileID, "/DE-diffxpy/DE_wilcoxon_", fileID, "_DE_", leiden_res, "_top100.csv.gz")) %>%
                                  slice(1:20)
  
  DE_genes = de_wilcoxon %>%
                tidyr::pivot_longer(everything(), names_to="clusterID", values_to="geneID") %>%
                arrange(clusterID) %>%
                pull(geneID) %>%
                unique()  
  
  tbl_expr = get_adata_raw_expr(DE_genes)
  
  
  tbl_DE_expr = metadata %>%
                  select(cellID, !!(leiden_res)) %>%
                  mutate(clusterID = get(leiden_res)) %>%
                  left_join(tbl_expr, by="cellID") %>%
                  tidyr::pivot_longer(starts_with("expr_"), names_to="geneID", values_to="expr") %>%
                  mutate(geneID = gsub("expr_", "", geneID)) %>%
                  group_by(geneID) %>%
                  mutate(expr_scaled = (expr - min(expr))/(max(expr) - min(expr))) %>%
                  mutate(expr_scaled_cut = ifelse(expr_scaled > 2.5, 2.5, ifelse(expr_scaled < -2.5, -2.5, expr_scaled))) %>%                  
                  ungroup() %>%
                  group_by(geneID, clusterID) %>%
                  summarise(frac_expressed = sum(expr > 0)/length(expr),
                            expr_average = mean(expr_scaled)) %>%
                  ungroup() %>%
                  mutate(geneID = factor(geneID, levels = rev(DE_genes)))
                  
  
  # expr_Gene
  plt_DE_dotplot = ggplot(tbl_DE_expr, aes(x=clusterID, y=geneID, size=frac_expressed, colour=expr_average)) +
                              geom_point() +
                              scale_size_continuous(range = c(0.1, 3)) +
                              theme_publication +
                              scale_colour_viridis() +
                              labs(x=paste0(leiden_res), y=NULL)
                              

return(plt_DE_dotplot)
}




plt_heatmap_scenic <- function(tbl_scenic_auc_z_score, tbl_scenic_rss_group, metadata){

  topreg = c()

  for (groupID in sort(tbl_scenic_rss_group$group)){
    top_group = tbl_scenic_rss_group %>%
                    filter(group == groupID) %>%
                    melt(id.vars="group") %>%
                    as_tibble() %>%
                    arrange(desc(value)) %>%
                    head(10) %>%
                    pull(variable) %>%
                    as.vector()
    topreg = c(topreg, top_group)
  }

  topreg = unique(topreg)

  tbl_reg_scores = tbl_scenic_auc_z_score %>%
                      select(cellID, !!topreg) %>%
                      melt(id.vars = "cellID") %>%
                      mutate(variable = factor(variable, levels = topreg)) %>%
                      left_join(metadata %>% select(cellID, group), by="cellID") %>%
                      arrange(group) %>%
                      mutate(cellID = factor(cellID, levels=unique(cellID))) %>%
                      mutate(value_cut = if_else(value > 2, 2, if_else(value < -2, -2, value)))

  tbl_group = metadata %>%
                  select(cellID, group) %>%
                  arrange(group) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID)))

  tbl_group_frac = tbl_group %>%
                        group_by(group) %>%
                        summarise(n_group = n()) %>%
                        mutate(frac_group = n_group/sum(n_group)) %>%
                        arrange(group) %>%
                        mutate(group = factor(group, levels=group)) %>%
                        mutate(x_max = cumsum(frac_group))

  tbl_group_frac = tbl_group_frac %>%
                      mutate(x_min = c(0, head(tbl_group_frac$x_max, -1))) %>%
                      mutate(x_mid = x_min + (x_max - x_min)/2) %>%
                      mutate(group_label_short = group)


  plt_group = ggplot(tbl_group_frac, aes(fill=group)) +
                      geom_rect(aes(xmin = x_min, xmax = x_max), ymin=-Inf, ymax=Inf) +
                      theme_publication +
                      theme(legend.position="none") +
                      # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                      scale_x_continuous(name=NULL, breaks=tbl_group_frac$x_mid, labels=(tbl_group_frac$group_label_short), expand = c(0,0),
                      guide = guide_axis(angle=45)) +
                      labs(x=NULL)
                      
                      
                      
  plt_cluster_marker_expr = ggplot(tbl_reg_scores, aes(x=cellID, y=variable, fill=value_cut)) +
                              geom_tile() +
                              # scale_fill_distiller(palette = 'RdYlBu', limits=c(-2, 6), na.value='chartreuse') +
                              scale_fill_distiller(palette = 'RdYlBu') +
                              theme_publication +
                              labs(x=NULL, y=NULL) +
                              scale_x_discrete(expand = c(0,0)) +
                              scale_y_discrete(expand = c(0,0)) +
                              labs(fill="AUC [z-score]") +
                              # theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank()) +
                              guides(fill = guide_colorbar(barwidth = 0.5,
                                                           barheight = 3)) +
                               theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) 

  plt_test = (plt_cluster_marker_expr / plt_group) + plot_layout(heights = c(1, 0.02))

return(plt_test)
}




plot_DE_dpt_heatmap <- function(fileID, pseudotime_DE){
    
  DE_genes = pseudotime_DE
  tbl_expr = get_adata_raw_expr(pseudotime_DE)
  
  cmat = as.matrix(tbl_expr %>% select(-cellID))
  rownames(cmat) = tbl_expr %>% pull(cellID)

  c <- cor(cmat, method="spearman")
  d <- as.dist(1-c)  # because correlation is a similarity measure and
  mat_cluster_cols <- hclust(d)
  
  
  geneID_order = colnames(cmat)[mat_cluster_cols$order] %>%
                  gsub("expr_", "", .)

  
  
  tbl_DE_expr = metadata %>%
                  select(cellID, !!(leiden_res), dpt_pseudotime) %>%
                  mutate(clusterID = get(leiden_res)) %>%
                  left_join(tbl_expr, by="cellID") %>%
                  tidyr::pivot_longer(starts_with("expr_"), names_to="geneID", values_to="expr") %>%
                  mutate(geneID = gsub("expr_", "", geneID)) %>%
                  group_by(geneID) %>%
                  mutate(expr_scaled = (expr - min(expr))/(max(expr) - min(expr))) %>%
                  mutate(expr_zscore = scale(expr)) %>%
                  mutate(expr_scaled_zscore = scale(expr_scaled)) %>%
                  mutate(expr_zscore_cut = ifelse(expr_zscore > 2.5, 2.5, ifelse(expr_zscore < -2.5, -2.5, expr_zscore))) %>%
                  ungroup() %>%
                  arrange(clusterID) %>%
                  arrange(dpt_pseudotime) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID))) %>%
                  mutate(geneID = factor(geneID, levels = rev(geneID_order)))


  plt_DE_heatmap = ggplot(tbl_DE_expr, aes(x=cellID, y=geneID, fill=expr_zscore_cut)) +
                              geom_tile() +
                              scale_fill_distiller(palette = 'RdYlBu') +
                              theme_publication +
                              labs(x=NULL, y=NULL) +
                              scale_x_discrete(expand = c(0,0)) +
                              scale_y_discrete(expand = c(0,0)) +
                              labs(fill="Expr [cut z-score log1p CPM]") +
                              # theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank()) +
                              guides(colour = guide_colorbar(barwidth = 0.3,
                                                             barheight = 2))



  tbl_cluster = metadata %>%
                  select(cellID, !!(leiden_res), dpt_pseudotime) %>%
                  mutate(leiden_res = get(leiden_res)) %>%
                  arrange(dpt_pseudotime) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID)))

   plt_cluster = ggplot(tbl_cluster, aes(x=cellID, y="", fill=leiden_res)) +
                       geom_tile() +
                       theme_publication +
                       scale_x_discrete(expand = c(0,0)) +
                       theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank()) +
                       labs(x=NULL, y=NULL) +
                       theme(legend.position="none") +
                       # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                       # , guide                               labs(x=NULL, y=NULL) +
# = guide_axis(angle=45)) +
                       labs(x=NULL)

   plt_cluster_dpt = ggplot(tbl_cluster, aes(x=cellID, y="", fill=dpt_pseudotime)) +
                       geom_tile() +
                       theme_publication +
                       scale_fill_viridis() +
                       scale_x_discrete(expand = c(0,0)) +
                       theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank()) +
                       labs(x=NULL, y=NULL) +
                       theme(legend.position="none") +
                       # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                       # , guide                               labs(x=NULL, y=NULL) +
# = guide_axis(angle=45)) +
                       labs(x=NULL)


  plt_panel = (plt_DE_heatmap / plt_cluster / plt_cluster_dpt) + plot_layout(heights = c(1, 0.02, 0.02))

return(plt_panel)
}






plot_DE_monocle2_heatmap <- function(fileID, pseudotime_DE){
    
  DE_genes = pseudotime_DE
  tbl_expr = get_adata_raw_expr(pseudotime_DE)
  
  cmat = as.matrix(tbl_expr %>% select(-cellID))
  rownames(cmat) = tbl_expr %>% pull(cellID)

  c <- cor(cmat, method="spearman")
  d <- as.dist(1-c)  # because correlation is a similarity measure and
  mat_cluster_cols <- hclust(d)
  
  
  geneID_order = colnames(cmat)[mat_cluster_cols$order] %>%
                  gsub("expr_", "", .)

  
  
  tbl_DE_expr = metadata %>%
                  select(cellID, !!(leiden_res), monocle_pseudotime) %>%
                  mutate(clusterID = get(leiden_res)) %>%
                  left_join(tbl_expr, by="cellID") %>%
                  tidyr::pivot_longer(starts_with("expr_"), names_to="geneID", values_to="expr") %>%
                  mutate(geneID = gsub("expr_", "", geneID)) %>%
                  group_by(geneID) %>%
                  mutate(expr_scaled = (expr - min(expr))/(max(expr) - min(expr))) %>%
                  mutate(expr_zscore = scale(expr)) %>%
                  mutate(expr_scaled_zscore = scale(expr_scaled)) %>%
                  mutate(expr_zscore_cut = ifelse(expr_zscore > 2.5, 2.5, ifelse(expr_zscore < -2.5, -2.5, expr_zscore))) %>%
                  ungroup() %>%
                  arrange(clusterID) %>%
                  arrange(monocle_pseudotime) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID))) %>%
                  mutate(geneID = factor(geneID, levels = rev(geneID_order)))


  plt_DE_heatmap = ggplot(tbl_DE_expr, aes(x=cellID, y=geneID, fill=expr_zscore_cut)) +
                              geom_tile() +
                              scale_fill_distiller(palette = 'RdYlBu') +
                              theme_publication +
                              labs(x=NULL, y=NULL) +
                              scale_x_discrete(expand = c(0,0)) +
                              scale_y_discrete(expand = c(0,0)) +
                              labs(fill="Expr [cut z-score log1p CPM]") +
                              # theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank()) +
                              guides(colour = guide_colorbar(barwidth = 0.3,
                                                             barheight = 2))



  tbl_cluster = metadata %>%
                  select(cellID, !!(leiden_res), monocle_pseudotime) %>%
                  mutate(leiden_res = get(leiden_res)) %>%
                  arrange(monocle_pseudotime) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID)))

   plt_cluster = ggplot(tbl_cluster, aes(x=cellID, y="", fill=leiden_res)) +
                       geom_tile() +
                       theme_publication +
                       scale_x_discrete(expand = c(0,0)) +
                       theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank()) +
                       labs(x=NULL, y=NULL) +
                       theme(legend.position="none") +
                       # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                       # , guide                               labs(x=NULL, y=NULL) +
# = guide_axis(angle=45)) +
                       labs(x=NULL)

   plt_cluster_dpt = ggplot(tbl_cluster, aes(x=cellID, y="", fill=monocle_pseudotime)) +
                       geom_tile() +
                       theme_publication +
                       scale_fill_viridis() +
                       scale_x_discrete(expand = c(0,0)) +
                       theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank()) +
                       labs(x=NULL, y=NULL) +
                       theme(legend.position="none") +
                       # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                       # , guide                               labs(x=NULL, y=NULL) +
# = guide_axis(angle=45)) +
                       labs(x=NULL)


  plt_panel = (plt_DE_heatmap / plt_cluster / plt_cluster_dpt) + plot_layout(heights = c(1, 0.02, 0.02))

return(plt_panel)
}


plot_UMAP_expr <- function(geneID){
  
  tbl_expr = get_adata_raw_expr(geneID)
  
  metadata = metadata %>%
              left_join(tbl_expr, by="cellID") %>%
              mutate(expr_Gene = get(paste0("expr_", geneID)))
              
  # expr_Gene
  plt_UMAP_expr = ggplot(slice(metadata, sample(1:n())), aes(x=UMAP_dim1, y=UMAP_dim2, z=expr_Gene)) +
                              stat_summary_hex(bins=100) +
                              theme_publication +
                              scale_fill_viridis() +
                              labs(x="UMAP 1", y="UMAP 2") +
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank()) +
                              labs(fill=paste0(geneID, " [log1p]"))

  return(plt_UMAP_expr)
}




plot_DE_cluster_heatmap <- function(DE_genes, cluster_num_subset=c(), guide_angle=45){

  tbl_expr = get_adata_raw_expr(DE_genes)
  
  if (length(cluster_num_subset) != 0){
    metadata = metadata %>%
                            mutate(clusterID_num = get(leiden_res)) %>%
                            filter(clusterID_num %in% cluster_num_subset)
  }

  tbl_DE_expr = metadata %>%
                  select(cellID, !!(leiden_res_annot), !!(leiden_res)) %>%
                  mutate(clusterID_num = get(leiden_res)) %>%
                  mutate(clusterID = get(leiden_res_annot)) %>%
                  mutate(clusterID = factor(clusterID, levels=annot_order)) %>%
                  left_join(tbl_expr, by="cellID") %>%
                  tidyr::pivot_longer(starts_with("expr_"), names_to="geneID", values_to="expr") %>%
                  mutate(geneID = gsub("expr_", "", geneID)) %>%
                  group_by(geneID) %>%
                  mutate(expr_scaled = (expr - min(expr))/(max(expr) - min(expr))) %>%
                  mutate(expr_zscore = scale(expr)) %>%
                  mutate(expr_scaled_zscore = scale(expr_scaled)) %>%
                  mutate(expr_zscore_cut = ifelse(expr_zscore > 1, 1, ifelse(expr_zscore < -1, -1, expr_zscore))) %>%
                  ungroup() %>%
                  arrange(clusterID) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID))) %>%
                  mutate(geneID = factor(geneID, levels = DE_genes))


  plt_DE_heatmap = ggplot(tbl_DE_expr, aes(x=cellID, y=geneID, fill=expr_zscore_cut)) +
                              # geom_tile() +
                              geom_tile_rast(raster.dpi=300) +
                              scale_fill_distiller(palette = 'RdYlBu') +
                              # scale_fill_gradient2(low = "magenta", mid = "black", high = "yellow") +
                              theme_publication +
                              labs(x=NULL, y="Genes") +
                              scale_x_discrete(expand = c(0,0)) +
                              scale_y_discrete(expand = c(0,0)) +
                              labs(fill="Expr [z-score]") +
                              # theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))
                              theme(axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank()) +
                              guides(fill = guide_colorbar(barwidth = 0.3,
                                                             barheight = 2))+
                              theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))



  tbl_cluster = metadata %>%
                  select(cellID, !!(leiden_res), !!(leiden_res_annot)) %>%
                  mutate(leiden_res = get(leiden_res)) %>%
                  mutate(leiden_res_annot = get(leiden_res_annot)) %>%
                  mutate(cellID = factor(cellID, levels=unique(cellID)))

  tbl_cluster_frac = tbl_cluster %>%
                       group_by(leiden_res_annot) %>%
                       summarise(n_cluster = n()) %>%
                       mutate(frac_cluster = n_cluster/sum(n_cluster)) %>%
                       mutate(leiden_res_annot = factor(leiden_res_annot, levels=annot_order)) %>%
                       arrange(leiden_res_annot) %>%
                       mutate(x_max = cumsum(frac_cluster))

  tbl_cluster_frac = tbl_cluster_frac %>%
                     mutate(x_min = c(0, head(tbl_cluster_frac$x_max, -1))) %>%
                     mutate(x_mid = x_min + (x_max - x_min)/2)

   plt_cluster = ggplot(tbl_cluster_frac, aes(fill=leiden_res_annot)) +
                       geom_rect(aes(xmin = x_min, xmax = x_max), ymin=-Inf, ymax=Inf) +
                       theme_publication +
                       scale_fill_manual(values=colour_map, na.value="#CCCCCC") +
                       theme(legend.position="none") +
                       # theme(plot.margin=unit(c(0, 0, 0, 0), "cm")) +
                       scale_x_continuous(name=NULL, breaks=tbl_cluster_frac$x_mid, labels=(tbl_cluster_frac$leiden_res_annot), expand = c(0,0),
                              guide = guide_axis(angle=guide_angle)) +
                       labs(x=NULL)


  plt_panel = (plt_DE_heatmap / plt_cluster) + plot_layout(heights = c(1, 0.02))

return(plt_panel)
}
