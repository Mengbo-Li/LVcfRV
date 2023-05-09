## Color annotations for the heatmap
colorRange <- colorRampPalette(rev(c(brewer.pal(11, "PRGn")[2:6], brewer.pal(11, "PiYG")[7:10])))(100)
colorRange2 <- colorRampPalette(rev(brewer.pal(8, "Blues")[8:1]))(10)
# add.flag <- function(pheatmap, kept.labels, repel.degree) {
#   
#   # repel.degree = number within [0, 1], which controls how much 
#   #                space to allocate for repelling labels.
#   ## repel.degree = 0: spread out labels over existing range of kept labels
#   ## repel.degree = 1: spread out labels over the full y-axis
#   heatmap <- pheatmap$gtable
#   new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
#   # keep only labels in kept.labels, replace the rest with ""
#   new.label$label <- ifelse(new.label$label %in% kept.labels, 
#                             new.label$label, "")
#   # calculate evenly spaced out y-axis positions
#   repelled.y <- function(d, d.select, k = repel.degree){
#     # d = vector of distances for labels
#     # d.select = vector of T/F for which labels are significant
#     
#     # recursive function to get current label positions
#     # (note the unit is "npc" for all components of each distance)
#     strip.npc <- function(dd){
#       if(!"unit.arithmetic" %in% class(dd)) {
#         return(as.numeric(dd))
#       }
#       d1 <- strip.npc(dd$arg1)
#       d2 <- strip.npc(dd$arg2)
#       fn <- dd$fname
#       return(lazyeval::lazy_eval(paste(d1, fn, d2)))
#     }
#     full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
#     selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
#     return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
#                     to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
#                     length.out = sum(d.select)), 
#                 "npc"))
#   }
#   new.y.positions <- repelled.y(new.label$y, d.select = new.label$label != "")
#   new.flag <- segmentsGrob(x0 = new.label$x,
#                            x1 = new.label$x + unit(0.15, "npc"),
#                            y0 = new.label$y[new.label$label != ""],
#                            y1 = new.y.positions)
#   # shift position for selected labels
#   new.label$x <- new.label$x + unit(0.2, "npc")
#   new.label$y[new.label$label != ""] <- new.y.positions
#   # add flag to heatmap
#   heatmap <- gtable::gtable_add_grob(x = heatmap, grobs = new.flag, t = 4, l = 4)
#   # replace label positions in heatmap
#   heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
#   # plot result
#   grid.newpage()
#   grid.draw(heatmap)
#   # return a copy of the heatmap invisibly
#   invisible(heatmap)
# }
# ### Get order of genes after clustering
# genesInHeatOrder_v <- donorHeat$tree_col$labels[donorHeat$tree_row$order]
# printGeneLabels <- which(genesInHeatOrder_v %in% donorTop$Symbol[1:10])
# printGeneLabels <- genesInHeatOrder_v[printGeneLabels]
# add.flag(donorHeat, kept.labels = printGeneLabels, repel.degree = -0.4)


pathPval <- function(k, pCutoff = 0.05) {
  topkegg <- mutate(k, Pathway = str_c(rownames(k), Pathway, sep = ":")) %>% 
    mutate(Pathway = str_sub(Pathway, start = 6)) %>% 
    filter(P.Up <= pCutoff | P.Down <= pCutoff) %>% 
    mutate(class = ifelse(P.Up <= pCutoff, "Up", "Down")) %>% 
    dplyr::select(Pathway, P.Up, P.Down, class) %>% 
    pivot_longer(2:3) %>% 
    filter(value <= pCutoff) %>% 
    dplyr::select(-name) %>% 
    mutate(value = -log10(value)) %>% 
    mutate(class = factor(class, levels = c("Up", "Down"))) %>% 
    arrange(class, value) %>% 
    mutate(Pathway = factor(Pathway, levels = Pathway))
  ggplot(topkegg, aes(x = Pathway, y = value, fill = class)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.4) + 
    coord_flip() + 
    scale_fill_manual(values = brewer.pal(8, "Paired")[c(5, 3)]) + 
    labs(x = "", y = "-log10(P-value)", fill = "Direction") + 
    theme_pubr()
}


# goPval <- function(go, pCutoff = 0.1) {
#   topGO <- mutate(go, Term = str_c(rownames(go), Term, sep = ":")) %>% 
#     filter(adj.P.Up <= pCutoff | adj.P.Down <= pCutoff) %>% 
#     mutate(class = ifelse(adj.P.Up <= pCutoff, "Up", "Down")) %>% 
#     dplyr::select(Term, Ont, adj.P.Up, adj.P.Down, class) %>% 
#     pivot_longer(3:4) %>% 
#     filter(value <= pCutoff) %>% 
#     dplyr::select(-name) %>% 
#     mutate(value = -log10(value)) %>% 
#     mutate(class = factor(class, levels = c("Up", "Down"))) %>% 
#     arrange(Ont, class, value) %>% 
#     mutate(Term = factor(Term, levels = Term))
#   ggplot(topGO, aes(x = Term, y = value, fill = class)) +
#     geom_col() +
#     geom_hline(yintercept = 0) +
#     geom_hline(yintercept = -log10(0.05), color = "white", linetype = "dashed", size = 0.4) + 
#     coord_flip() + 
#     facet_wrap(~Ont, nrow = 1) + 
#     scale_fill_manual(values = brewer.pal(8, "Paired")[c(5, 3)]) + 
#     labs(x = "", y = "-log10(P-value)", fill = "Direction") + 
#     theme_pubr()
# }


pathPval2 <- function(k, pCutoff = 0.05) {
  topkegg <- filter(k, P.Up <= pCutoff | P.Down <= pCutoff) %>% 
    mutate(class = ifelse(P.Up <= pCutoff, "Up", "Down")) %>% 
    dplyr::select(Pathway, P.Up, P.Down, class) %>% 
    pivot_longer(2:3) %>% 
    filter(value <= pCutoff) %>% 
    dplyr::select(-name) %>% 
    mutate(value = -log10(value)) %>% 
    mutate(class = factor(class, levels = c("Up", "Down"))) %>% 
    arrange(class, value) %>% 
    mutate(Pathway = factor(Pathway, levels = Pathway))
  ggplot(topkegg, aes(x = Pathway, y = value, fill = class)) +
    geom_col() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.4) + 
    coord_flip() + 
    scale_fill_manual(values = brewer.pal(8, "Paired")[c(5, 3)]) + 
    labs(x = "", y = "-log10(P-value)", fill = "Direction") + 
    theme_pubr()
}





geneOverlap <- function(k, theToptable, up.Pcut, down.Pcut) {
  topkegg <- mutate(k, PathwayID = rownames(k), .before = 1) %>% 
    filter(P.Up <= 0.05 | P.Down <= 0.05)
  GK <- getGeneKEGGLinks(species.KEGG = "hsa")
  PathwayNames <- getKEGGPathwayNames(species.KEGG = "hsa", remove.qualifier = TRUE)
  ups <- filter(topkegg, P.Up <= up.Pcut)
  upPath_genes <- lapply(ups$PathwayID, function(pi) {
    ID <- GK$GeneID[GK$Pathway == pi]
    genes <- filter(theToptable, theToptable$GeneID %in% ID, adj.P.Val <= 0.05)
    up_genes <- genes$Symbol[genes$logFC >= 0]
    # down_genes <- genes$Symbol[genes$logFC < 0]
    data.frame(PathwayID = pi, gene = up_genes)
  }) %>% do.call(rbind, .) %>% left_join(PathwayNames, by = "PathwayID") %>% arrange(gene)
  downs <- filter(topkegg, P.Down <= down.Pcut)
  downPath_genes <- lapply(downs$PathwayID, function(pi) {
    ID <- GK$GeneID[GK$Pathway == pi]
    genes <- filter(theToptable, theToptable$GeneID %in% ID, adj.P.Val <= 0.05)
    # up_genes <- genes$Symbol[genes$logFC >= 0]
    down_genes <- genes$Symbol[genes$logFC < 0]
    data.frame(PathwayID = pi, gene = down_genes)
  }) %>% do.call(rbind, .) %>% left_join(PathwayNames, by = "PathwayID") %>% arrange(gene)
  
  upGenePlot <- ggplot(upPath_genes, aes(axis1 = PathwayID, axis2 = gene, fill = Description)) + 
    geom_alluvium(alpha = 0.7) + 
    geom_stratum(alpha = 1, color = "white", show.legend = FALSE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.8) +
    scale_fill_manual(values = brewer.pal(8, "Pastel1")) + 
    scale_color_manual(values = brewer.pal(8, "Pastel1")) + 
    theme_void()
  downGenePlot <- ggplot(downPath_genes, aes(axis1 = PathwayID, axis2 = gene, fill = Description)) + 
    geom_alluvium(alpha = 0.7) + 
    geom_stratum(alpha = 1, color = "white", show.legend = FALSE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.8) +
    scale_fill_manual(values = brewer.pal(8, "Pastel1")) + 
    scale_color_manual(values = brewer.pal(8, "Pastel1")) + 
    theme_void()
  list(upPath_genes = upPath_genes, upGenePlot = upGenePlot, downPath_genes = downPath_genes, downGenePlot = downGenePlot)
}



getPathGenes <- function(k, theToptable, up.Pcut, down.Pcut) {
  topkegg <- mutate(k, PathwayID = rownames(k), .before = 1) %>% 
    filter(P.Up <= 0.05 | P.Down <= 0.05)
  GK <- getGeneKEGGLinks(species.KEGG = "hsa")
  PathwayNames <- getKEGGPathwayNames(species.KEGG = "hsa", remove.qualifier = TRUE)
  ups <- filter(topkegg, P.Up <= up.Pcut) %>% arrange(P.Up)
  upPath_genes <- lapply(ups$PathwayID, function(pi) {
    ID <- GK$GeneID[GK$Pathway == pi]
    genes <- filter(theToptable, theToptable$GeneID %in% ID, adj.P.Val <= 0.05)
    up_genes <- genes$Symbol[genes$logFC >= 0]
    down_genes <- genes$Symbol[genes$logFC < 0]
    if (length(down_genes) > 0) data.frame(PathwayID = pi, gene = up_genes, direction = "Up") %>% 
      add_row(data.frame(PathwayID = pi, gene = down_genes, direction = "Down"))
    else data.frame(PathwayID = pi, gene = up_genes, direction = "Up")
  }) %>% do.call(rbind, .) %>% left_join(PathwayNames, by = "PathwayID") %>% 
    arrange(PathwayID) %>% mutate(Path.direction = "Up-reg")
  downs <- filter(topkegg, P.Down <= down.Pcut) %>% arrange(P.Down)
  downPath_genes <- lapply(downs$PathwayID, function(pi) {
    ID <- GK$GeneID[GK$Pathway == pi]
    genes <- filter(theToptable, theToptable$GeneID %in% ID, adj.P.Val <= 0.05)
    up_genes <- genes$Symbol[genes$logFC >= 0]
    down_genes <- genes$Symbol[genes$logFC < 0]
    if (length(up_genes) > 0) data.frame(PathwayID = pi, gene = down_genes, direction = "Down") %>% 
      add_row(data.frame(PathwayID = pi, gene = up_genes, direction = "Up"))
    else data.frame(PathwayID = pi, gene = down_genes, direction = "Down")
  }) %>% do.call(rbind, .) %>% left_join(PathwayNames, by = "PathwayID") %>%
    arrange(PathwayID) %>% mutate(Path.direction = "Down-reg")
  pathDE <- add_row(upPath_genes, downPath_genes) %>% 
    dplyr::rename(Symbol = gene, Gene.direction = direction) %>% 
    dplyr::select(1, 4, 2, 3, 5)
  pathDE
}


ggMDS <- function(dat, meta, label) {
  mds <- plotMDS(t(dat), plot = FALSE)
  data.frame(smpName = rownames(dat), x = mds$x, y = mds$y) %>% 
    left_join(meta, by = "smpName") %>% 
    mutate(label = label)
}
ggMDS2 <- function(dat, meta, label) {
  mds <- plotMDS(t(dat), plot = FALSE)
  data.frame(sample = rownames(dat), x = mds$x, y = mds$y) %>% 
    left_join(meta, by = "sample") %>% 
    mutate(label = label)
}

# get first p dimensions in MDS
getMDS <- function(data) {
  mds <- plotMDS(data, plot = FALSE)
  mds$distance.matrix.squared
  # mds$eigen.vectors
}

# knn graph
library(igraph)
library(mstknnclust)
make.knn.graph <- function(D, edgeBtwCut = 0.9){
  # calculate euclidean distances between cells
  dist <- as.matrix(dist(D))
  # # make a list of edges to k nearest neighbors for each cell
  # edges <- mat.or.vec(0, 2)
  # for (i in 1:nrow(dist)){
  #   # find closes neighbours
  #   matches <- setdiff(order(dist[i, ], decreasing = FALSE)[1:(k+1)], i)
  #   # add edges in both directions
  #   edges <- rbind(edges, cbind(rep(i, k), matches))  
  #   edges <- rbind(edges, cbind(matches, rep(i, k)))
  # }
  # # create a graph from the edgelist
  # graph <- graph_from_edgelist(edges, directed = FALSE)
  # # V(graph)$frame.color <- NA
  graph <- mst.knn(dist)
  graph <- graph$network
  # membership <- igraph::clusters(graph)$membership
  # filter nodes without edge
  V(graph)$label <- V(graph)$name <- rownames(D)[as.numeric(V(graph)$name)]
  # degree of graph
  graphdegree <- igraph::degree(graph, v = V(graph), loops = TRUE, normalized = TRUE)
  # graphdegree <- knn(graph, vids = V(graph), weights = NULL)
  # edge betweenness
  edgeBtw <- edge_betweenness(graph, e = E(graph), directed = FALSE, weights = NULL)
  edgeBtw <- as.numeric(scale(edgeBtw, center = FALSE))
  # graph <- subgraph.edges(graph, eids = which(edgeBtw > quantile(edgeBtw, prob = filterEdge)), delete.vertices = TRUE)
  # edgeBtw <- edgeBtw[edgeBtw > median(edgeBtw)]
  col.edge <- c("#C8C8C8", "#808080")[(edgeBtw > quantile(edgeBtw, prob = edgeBtwCut))+1]
  # make a layout for visualizing in 2D
  set.seed(1)
  g.layout <- layout_with_fr(graph)
  return(list(graph = graph, layout = g.layout, graphdegree = graphdegree, edgeBtw = edgeBtw, 
              col.edge = col.edge, dist = dist))        
}


gls.series <- function(M, design = NULL, ndups = 2, spacing = 1, block = NULL, 
                        correlation = NULL, weights = NULL, 
                        cor.trend = FALSE, ###### added for the case where dupcor needs to be run within lmFit
                        trim = 0.5, ###### temporarily added for QR calculation when correlation is a vector, so that QR returned
                        ##### is from a consensus correlation calculated from the input vector, instead of a series of QR's
                        ...) {
  M <- as.matrix(M)
  ngenes <- nrow(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- matrix(1, narrays, 1)
  design <- as.matrix(design)
  if (nrow(design) != narrays) 
    stop("Number of rows of design matrix does not match number of arrays")
  nbeta <- ncol(design)
  coef.names <- colnames(design)
  
  if (is.null(correlation)) {
    ######################
    correlation <- duplicateCorrelation(M, design = design, 
                                        ndups = ndups, spacing = spacing, block = block, 
                                        weights = weights, ...)
    #### For tidiness, there is no smoothing option added here, should be put into duplicateCorrelation later
    if (cor.trend) correlation <- correlation$atanh.correlations
    correlation <- correlation$consensus.correlation
    ######################
  }
  
  # if (abs(correlation) >= 1) 
  if (abs(correlation) >= 1 || sum( abs(correlation) >= atanh(1)) >= 1) 
    stop("correlation is 1 or -1, so the model is degenerate")
  #####################
  
  if (!is.null(weights)) {
    weights[is.na(weights)] <- 0
    weights <- asMatrixWeights(weights, dim(M))
    M[weights < 1e-15] <- NA
    weights[weights < 1e-15] <- NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates (ndups<2)")
      ndups <- 1
      correlation <- 0
    }
    
    #####################
    # cormatrix <- diag(rep_len(correlation, narrays), nrow = narrays, 
    #                   ncol = narrays) %x% array(1, c(ndups, ndups))
    ### correlation[1] instead of correlation to accommodate the genewise correlation case 
    cormatrix <- diag(rep_len(correlation[1], narrays), nrow = narrays, 
                      ncol = narrays) %x% array(1, c(ndups, ndups))
    #####################
    
    if (is.null(spacing)) 
      spacing <- 1
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    if (!is.null(weights)) 
      weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep_len(1, ndups)
    colnames(design) <- coef.names
    ngenes <- nrow(M)
    narrays <- ncol(M)
  }
  else {
    ndups <- spacing <- 1
    block <- as.vector(block)
    if (length(block) != narrays) 
      stop("Length of block does not match number of arrays")
    ub <- unique(block)
    nblocks <- length(ub)
    Z <- matrix(block, narrays, nblocks) == matrix(ub, narrays, 
                                                   nblocks, byrow = TRUE)
    
    #####################
    # cormatrix <- Z %*% (correlation * t(Z))
    ### correlation[1] instead of correlation to accommodate the genewise correlation case 
    ### In the genewise case, this is only a placeholder and no tanh-transformation is applied yet 
    cormatrix <- Z %*% (correlation[1] * t(Z))
    #####################
  }
  diag(cormatrix) <- 1
  stdev.unscaled <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M), 
                                                              coef.names))
  
  #######################
  NoProbeWts <- all(is.finite(M)) && (is.null(weights) ||
                                        !is.null(attr(weights, "arrayweights")))
  ##### this fast calculation only used for consensus correlation + arrayweights or no weights
  ##### and there must not be any NA in M
  ##### substitute the if-condition
  # if (NoProbeWts) {
  if (NoProbeWts && length(correlation) == 1L) { 
    ########################
    
    V <- cormatrix
    if (!is.null(weights)) {
      wrs <- 1/sqrt(weights[1, ])
      V <- wrs * t(wrs * t(V))
    }
    cholV <- chol(V)
    y <- backsolve(cholV, t(M), transpose = TRUE)
    dimnames(y) <- rev(dimnames(M))
    X <- backsolve(cholV, design, transpose = TRUE)
    dimnames(X) <- dimnames(design)
    fit <- lm.fit(X, y)
    if (fit$df.residual > 0) {
      if (is.matrix(fit$effects)) 
        fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank), 
                                               , drop = FALSE]^2))
      else fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
    }
    else fit$sigma <- rep_len(NA_real_, ngenes)
    fit$fitted.values <- fit$residuals <- fit$effects <- NULL
    fit$coefficients <- t(fit$coefficients)
    fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
    est <- fit$qr$pivot[1:fit$qr$rank]
    dimnames(fit$cov.coefficients) <- list(coef.names[est], 
                                           coef.names[est])
    stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)), 
                                    ngenes, fit$qr$rank, byrow = TRUE)
    fit$stdev.unscaled <- stdev.unscaled
    fit$df.residual <- rep_len(fit$df.residual, length.out = ngenes)
    dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
    fit$pivot <- fit$qr$pivot
    fit$ndups <- ndups
    fit$spacing <- spacing
    fit$block <- block
    fit$correlation <- correlation
    return(fit)
  }
  
  
  
  ########################
  ####### when correlation is genewise and/or M has observational/gene weights, need to loop through all genes
  ####### observational weights (same dim as M) or gene weights (length ngenes vector) are the same in implementation
  ####### as gene weights were already converted into the matrix format
  ####### in the notes below just mention them both as the observational weights for simplicity
  beta <- stdev.unscaled
  sigma <- rep_len(NA_real_, ngenes)
  df.residual <- rep_len(0, ngenes)
  
  ########################
  ####### for consensus cor + observation weights - no need to update `V` for each gene
  ####### Or for consensus + arrayweights/no weights, but there are missing values in Y
  if (length(correlation) == 1L) {
    
    for (i in 1:ngenes) {
      y <- drop(M[i, ])
      o <- is.finite(y)
      y <- y[o]
      n <- length(y)
      if (n > 0) {
        X <- design[o, , drop = FALSE]
        V <- cormatrix[o, o]
        
        #### observational weights
        if (!is.null(weights)) {
          wrs <- 1/sqrt(drop(weights[i, o]))
          V <- wrs * t(wrs * t(V))
        }
        ####### weights end  ------
        
        cholV <- chol(V)
        y <- backsolve(cholV, y, transpose = TRUE)
        if (all(X == 0)) {
          df.residual[i] <- n
          sigma[i] <- sqrt(array(1/n, c(1, n)) %*% y^2)
        }
        else {
          X <- backsolve(cholV, X, transpose = TRUE)
          out <- lm.fit(X, y)
          est <- !is.na(out$coefficients)
          beta[i, ] <- out$coefficients
          stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, 
                                                       size = out$rank)))
          df.residual[i] <- out$df.residual
          if (df.residual[i] > 0) 
            sigma[i] <- sqrt(array(1/out$df.residual, 
                                   c(1, n)) %*% out$residuals^2)
        }
      }
    }
    
  } 
  
  ###### lastly when using genewise correlations, with each type of weights: array, observational or none
  ###### also allow NAs in y
  ###### need to update V and weights in every iteration
  if (length(correlation) > 1L) {
    
    ####################################
    # transform atanh-correlations back to correlation values
    correlation.raw <- tanh(correlation)
    ####################################
    
    for (i in 1:ngenes) {
      y <- drop(M[i, ])
      o <- is.finite(y)
      y <- y[o]
      n <- length(y)
      if (n > 0) {
        X <- design[o, , drop = FALSE]
        V <- cormatrix[o, o]
        
        ######## substitute genewise correlations
        V[V != 0] <- correlation.raw[i]
        diag(V) <- 1
        
        ####### weights start  ------
        if (!is.null(weights)) {
          wrs <- 1/sqrt(drop(weights[i, o]))
          V <- wrs * t(wrs * t(V))
        }
        ####### weights end  ------
        
        ############ should all be solvable with lower bound on genewise correlations
        cholV <- tryCatch(chol(V), error = function(e) NA)
        y <- backsolve(cholV, y, transpose = TRUE)
        if (is.na(y[1])) { 
          warning(paste0("Model cannot be solved: Gene ", rownames(M)[i]))
          next
        }
        if (all(X == 0)) {
          df.residual[i] <- n
          sigma[i] <- sqrt(array(1/n, c(1, n)) %*% y^2)
        }
        else {
          X <- backsolve(cholV, X, transpose = TRUE)
          out <- lm.fit(X, y)
          est <- !is.na(out$coefficients)
          beta[i, ] <- out$coefficients
          stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, 
                                                       size = out$rank)))
          df.residual[i] <- out$df.residual
          if (df.residual[i] > 0) 
            sigma[i] <- sqrt(array(1/out$df.residual, 
                                   c(1, n)) %*% out$residuals^2)
        }
      }
    }
    
    ################ for now, QR is returned with the consensus correlation
    ######## or do we want to return QR for each V? => the output list will be gigantic
    ######## NOTE: in this implementation, as the consensus correlation is calculated from the input vector, 
    ######## the trimmed-mean (consensus cor) is calculated on the smoothed values if the input correlations are 
    ######## already smoothed. (Check w Gordon)
    cormatrix[cormatrix != 0] <- tanh(mean(correlation, trim = trim, na.rm = TRUE))
    diag(cormatrix) <- 1
    ##########################
  }
  
  cholV <- chol(cormatrix)
  QR <- qr(backsolve(cholV, design, transpose = TRUE))
  cov.coef <- chol2inv(QR$qr, size = QR$rank)
  est <- QR$pivot[1:QR$rank]
  dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
  list(coefficients = beta, stdev.unscaled = stdev.unscaled, 
       sigma = sigma, df.residual = df.residual, ndups = ndups, 
       spacing = spacing, block = block, correlation = correlation, 
       cov.coefficients = cov.coef, pivot = QR$pivot, rank = QR$rank)
}









## add to lmFit 
## this is a temporary modification on lmFit to replace gls.series by gls.series

lmFit <- function (object, design = NULL, ndups = 1, spacing = 1, block = NULL, 
                    correlation, weights = NULL, method = "ls", ...) 
{
  y <- getEAWP(object)
  if (!nrow(y$exprs)) 
    stop("expression matrix has zero rows")
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
    if (nrow(design) != ncol(y$exprs)) 
      stop("row dimension of design doesn't match column dimension of data object")
  }
  ne <- nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups <- y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing <- y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights <- y$weights
  method <- match.arg(method, c("ls", "robust"))
  if (ndups > 1) {
    if (!is.null(y$probes)) 
      y$probes <- uniquegenelist(y$probes, ndups = ndups, 
                                 spacing = spacing)
    if (!is.null(y$Amean)) 
      y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), 
                                     ndups = ndups, spacing = spacing), na.rm = TRUE)
  }
  if (method == "robust") 
    fit <- mrlm(y$exprs, design = design, ndups = ndups, 
                spacing = spacing, weights = weights, ...)
  else if (ndups < 2 && is.null(block)) 
    fit <- lm.series(y$exprs, design = design, ndups = ndups, 
                     spacing = spacing, weights = weights)
  else {
    if (missing(correlation)) 
      stop("the correlation must be set, see duplicateCorrelation")
    fit <- gls.series(y$exprs, design = design, ndups = ndups, 
                       spacing = spacing, block = block, correlation = correlation, 
                       weights = weights, ...)
  }
  if (NCOL(fit$coefficients) > 1) {
    n <- rowSums(is.na(fit$coefficients))
    n <- sum(n > 0 & n < NCOL(fit$coefficients))
    if (n > 0) 
      warning("Partial NA coefficients for ", n, " probe(s)", 
              call. = FALSE)
  }
  fit$genes <- y$probes
  fit$Amean <- y$Amean
  fit$method <- method
  fit$design <- design
  new("MArrayLM", fit)
}

#  DUPS.R
#  Functions to handle duplicate spots or blocking

unwrapdups <- function(M,ndups=2,spacing=1) {
  #	Unwrap M matrix for a series of experiments so that all spots for a given gene are in one row
  #	Gordon Smyth
  #	18 Jan 2002. Last revised 2 Nov 2002.
  
  if(ndups==1) return(M)
  M <- as.matrix(M)
  nspots <- dim(M)[1]
  nslides <- dim(M)[2]
  ngroups <- nspots / ndups / spacing
  dim(M) <- c(spacing,ndups,ngroups,nslides)
  M <- aperm(M,perm=c(1,3,2,4))
  dim(M) <- c(spacing*ngroups,ndups*nslides)
  M
}

uniquegenelist <- function(genelist,ndups=2,spacing=1) {
  #	Eliminate entries in genelist for duplicate spots
  #	Gordon Smyth
  #	2 Nov 2002.  Last revised 10 Jan 2005
  
  if(ndups <= 1) return(genelist)
  i <- drop(unwrapdups(1:NROW(genelist),ndups=ndups,spacing=spacing)[,1])
  if(is.null(dim(genelist)))
    return(genelist[i])
  else
    return(genelist[i,,drop=FALSE])
}

duplicateCorrelation2 <- function(object,design=NULL,ndups=2L,spacing=1L,block=NULL,trim=0.15,weights=NULL)
  #	Estimate the correlation between duplicates given a series of arrays
  #	Gordon Smyth
  #	25 Apr 2002. Last revised 15 Feb 2021.
{
  #	Extract components from y
  y <- getEAWP(object)
  M <- y$exprs
  ngenes <- nrow(M)
  narrays <- ncol(M)
  
  #	Check design matrix
  if(is.null(design)) design <- y$design
  if(is.null(design))
    design <- matrix(1,ncol(y$exprs),1)
  else {
    design <- as.matrix(design)
    if(!identical(mode(design),"numeric")) stop("design must be a numeric matrix")
  }
  if(!identical(nrow(design),narrays)) stop("Number of rows of design matrix does not match number of arrays")
  nbeta <- ncol(design)
  
  #	Check whether design and block are of full rank
  QR <- qr(design)
  if(QR$rank < nbeta) message("Note: design matrix not of full rank (",nbeta-QR$rank," coef not estimable).")
  if(!is.null(block)) {
    MaxBlockSize <- max(table(block))
    if(identical(MaxBlockSize,1L)) {
      warning("Blocks all of size 1: setting intrablock correlation to zero.")
      return( list(consensus.correlation=0,cor=0,atanh.correlations=rep_len(0,nrow(M))) )
    }
    design.block <- model.matrix(~factor(block))
    design.block <- design.block[,-1,drop=FALSE]
    QtBlock <- qr.qty(QR,design.block)
    if(max(abs(QtBlock[-(1:QR$rank),])) < 1e-8) {
      warning("Block factor already encoded in the design matrix: setting intrablock correlation to zero.")
      return( list(consensus.correlation=0,cor=0,atanh.correlations=rep_len(0,nrow(M))) )
    }
  }
  
  #	Check weights
  if(is.null(weights)) weights <- y$weights
  if(!is.null(weights)) {
    weights <- asMatrixWeights(weights,dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  
  #	Setup spacing or blocking arguments
  if(is.null(block)) {
    #		If present, use ndups and spacing stored in object
    if(!is.null(y$printer$ndups)) ndups <- y$printer$ndups
    if(!is.null(y$printer$spacing)) spacing <- y$printer$spacing
    if(ndups<2L) {
      warning("No duplicates: setting correlation between duplicates to zero.")
      return( list(consensus.correlation=0,cor=0,atanh.correlations=rep_len(0,nrow(M))) )
    }
    if(is.character(spacing)) {
      if(spacing=="columns") spacing <- 1
      if(spacing=="rows") spacing <- y$printer$nspot.c
      if(spacing=="topbottom") spacing <- nrow(M)/2
    }
    Array <- rep(1:narrays,each=ndups)
  } else {
    ndups <- 1L
    nspacing <- 1L
    Array <- block
  }
  
  #	Unwrap data to get all data for a gene in one row
  if(is.null(block)) {
    M <- unwrapdups(M,ndups=ndups,spacing=spacing)
    ngenes <- nrow(M)
    if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
    design <- design %x% rep_len(1,ndups)
  }
  
  #	Compute genewise correlations
  if(!requireNamespace("statmod",quietly=TRUE)) stop("statmod package required but is not installed (or can't be loaded)")
  nafun <- function(e) NA
  rho <- unlist(parallel::mclapply(1:ngenes, mc.cores = 12, mc.preschedule = FALSE, function(i) {
    y <- drop(M[i,])
    o <- is.finite(y)
    A <- factor(Array[o])
    nobs <- sum(o)
    nblocks <- length(levels(A))
    if(nobs>(nbeta+2) && nblocks>1 && nblocks<nobs-1) {
      y <- y[o]
      X <- design[o,,drop=FALSE]
      Z <- model.matrix(~0+A)
      if(!is.null(weights)) {
        w <- drop(weights[i,])[o]
        s <- tryCatch(statmod::mixedModel2Fit(y,X,Z,w,only.varcomp=TRUE,maxit=20)$varcomp,error=nafun)
      } else
        s <- tryCatch(statmod::mixedModel2Fit(y,X,Z,only.varcomp=TRUE,maxit=20)$varcomp,error=nafun)
      if(!is.na(s[1])) { return(s[2]/sum(s)) } else { return(NA)}
    } else return(NA)
  }))
  
  #	Keep correlations away from limits to ensure correlation matrix is positive-definite
  rhomax <- 0.99
  if(is.null(block))
    rhomin <- 1/(1-ndups) + 0.01
  else
    rhomin <- 1/(1-MaxBlockSize) + 0.01
  if(min(rho, na.rm = TRUE) < rhomin) rho[rho < rhomin] <- rhomin
  if(max(rho, na.rm = TRUE) > rhomax) rho[rho > rhomax] <- rhomax
  
  arho <- atanh(rho)
  mrho <- tanh(mean(arho,trim=trim,na.rm=TRUE))
  list(consensus.correlation=mrho,cor=mrho,atanh.correlations=arho)
}

avedups <- function(x,ndups,spacing,weights) UseMethod("avedups")

avedups.default <- function(x,ndups=2,spacing=1,weights=NULL)
  #	Average over duplicate spots, for matrices or vectors
  #	Gordon Smyth
  #	6 Apr 2006.
{
  if(ndups==1) return(x)
  if(is.null(x)) return(NULL)
  x <- as.matrix(x)
  nspots <- dim(x)[1]
  nslides <- dim(x)[2]
  rn <- rownames(x)
  cn <- colnames(x)
  ngroups <- nspots / ndups / spacing
  dim(x) <- c(spacing,ndups,ngroups*nslides)
  x <- aperm(x,perm=c(2,1,3))
  if(mode(x)=="character")
    x <- x[1,,]
  else {
    if(is.null(weights))
      x <- colMeans(x,na.rm=TRUE)
    else {
      weights <- as.matrix(weights)
      dim(weights) <- c(spacing,ndups,ngroups*nslides)
      weights <- aperm(weights,perm=c(2,1,3))
      weights[is.na(weights) | is.na(x)] <- 0
      weights[weights<0] <- 0
      x <- colSums(weights*x,na.rm=TRUE)/colSums(weights)
    }
  }
  dim(x) <- c(spacing*ngroups,nslides)
  colnames(x) <- cn
  rownames(x) <- avedups(rn,ndups=ndups,spacing=spacing)
  x
}

avedups.MAList <- function(x,ndups=x$printer$ndups,spacing=x$printer$spacing,weights=x$weights)
  #	Average over duplicate spots for MAList objects
  #	Gordon Smyth
  #	6 Apr 2006.
{
  if(is.null(ndups) || is.null(spacing)) stop("Must specify ndups and spacing")
  y <- x
  y$M <- avedups(x$M,ndups=ndups,spacing=spacing,weights=weights)
  y$A <- avedups(x$A,ndups=ndups,spacing=spacing,weights=weights)
  other <- names(x$other)
  for (a in other) object$other[[a]] <- avedups(object$other[[a]],ndups=ndups,spacing=spacing,weights=weights)
  y$weights <- avedups(x$weights,ndups=ndups,spacing=spacing)
  y$genes <- uniquegenelist(x$genes,ndups=ndups,spacing=spacing)
  y$printer <- NULL
  y
}

avedups.EList <- function(x,ndups=x$printer$ndups,spacing=x$printer$spacing,weights=x$weights)
  #	Average over duplicate spots for EList objects
  #	Gordon Smyth
  #	2 Apr 2010.
{
  if(is.null(ndups) || is.null(spacing)) stop("Must specify ndups and spacing")
  y <- x
  y$E <- avedups(x$E,ndups=ndups,spacing=spacing,weights=weights)
  other <- names(x$other)
  for (a in other) object$other[[a]] <- avedups(object$other[[a]],ndups=ndups,spacing=spacing,weights=weights)
  y$weights <- avedups(x$weights,ndups=ndups,spacing=spacing)
  y$genes <- uniquegenelist(x$genes,ndups=ndups,spacing=spacing)
  y$printer <- NULL
  y
}

avereps <- function(x,...)
  #	4 June 2008
  UseMethod("avereps")

avereps.default <- function(x,ID=rownames(x),...)
  #	Average over irregular replicate spots, for matrices or vectors
  #	Gordon Smyth
  #	Created 3 June 2008.  Last modified 1 Dec 2010.
  #	Revised 19 Aug 2009 following suggestions from Axel Klenk.
  #	Revised 28 March 2010 following suggestion from Michael Lawrence.
{
  if(is.null(x)) return(NULL)
  x <- as.matrix(x)
  if(is.null(ID)) stop("No probe IDs")
  ID <- as.character(ID)
  if(mode(x)=="character") {
    d <- duplicated(ID)
    if(!any(d)) return(x)
    y <- x[!d,,drop=FALSE]
    return(y)
  }
  ID <- factor(ID,levels=unique(ID))
  #	rowsum(x,ID,reorder=FALSE,na.rm=TRUE)/as.vector(table(ID))
  y <- rowsum(x,ID,reorder=FALSE,na.rm=TRUE)
  n <- rowsum(1L-is.na(x),ID,reorder=FALSE)
  y/n
}

avereps.MAList <- function(x,ID=NULL,...)
  #	Average over irregular replicate spots for MAList objects
  #	Gordon Smyth
  #	3 June 2008.  Last modified 8 Sep 2010.
{
  if(is.null(ID)) {
    ID <- x$genes$ID
    if(is.null(ID)) ID <- rownames(x)
    if(is.null(ID)) stop("Cannot find probe IDs")
  }
  y <- x
  y$M <- avereps(x$M,ID=ID)
  y$A <- avereps(x$A,ID=ID)
  other <- names(x$other)
  for (a in other) y$other[[a]] <- avereps(x$other[[a]],ID=ID)
  y$weights <- avereps(x$weights,ID=ID)
  y$genes <- x$genes[!duplicated(ID),]
  y$printer <- NULL
  y
}

avereps.EList <- function(x,ID=NULL,...)
  #	Average over irregular replicate probes for EList objects
  #	Gordon Smyth
  #	2 April 2010.  Last modified 20 May 2011.
{
  if(is.null(ID)) {
    ID <- x$genes$ID
    if(is.null(ID)) ID <- rownames(x)
    if(is.null(ID)) stop("Cannot find probe IDs")
  }
  y <- x
  y$E <- avereps(x$E,ID=ID)
  other <- names(x$other)
  for (a in other) y$other[[a]] <- avereps(x$other[[a]],ID=ID)
  y$weights <- avereps(x$weights,ID=ID)
  y$genes <- x$genes[!duplicated(ID),]
  y$printer <- NULL
  y
}

avereps.RGList <- function(x,ID=NULL,...)
  #	Warn users that averaging should not be applied prior to normalization
  #	Gordon Smyth
  #	2 December 2013.
{
  stop("avereps should not be applied to an RGList object")
  invisible()
}


avereps.EListRaw <- function(x,ID=NULL,...)
  #	Warn users that averaging should not be applied prior to normalization
  #	Gordon Smyth
  #	2 December 2013.
{
  stop("avereps should not be applied to an EListRaw object")
  invisible()
}
