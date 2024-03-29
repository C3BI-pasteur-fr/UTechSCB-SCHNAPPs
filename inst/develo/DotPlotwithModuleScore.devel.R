require(cowplot)
require(Seurat)
require(ggplot2)
require(tidyr)
require(cowplot)
#' DotPlotwithModuleScore
#' 
#' adapted from Seurat::DotPlot but adds the module score per group
#' 
#' @param features list of character vectors of gene names in "symbol" annotation (SCHNAPPs specific)
#' 
#' @export DotPlotwithModuleScore
#' 

DotPlotwithModuleScore <- function (object, assay = NULL, features,
                                    featureDat = featureDat, 
                                    cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                                    idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, clusters="sampleNames",
                                    scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  
  # object = seurDat
  # assay="RNA"
  # features = features
  # cols = col
  # col.min = col.min
  # col.max = col.max
  # dot.min = dot.min
  # dot.scale = dot.scale
  # idents = NULL
  # group.by = NULL
  # split.by = NULL
  # cluster.idents = FALSE
  # scale = TRUE
  # scale.by = scale.by
  # scale.min = NA
  # scale.max = NA
  
  # save(file = "~/SCHNAPPsDebug/DotPlotwithModuleScoreF.RData", list = c(ls()))
  # cp = load("~/SCHNAPPsDebug/DotPlotwithModuleScoreF.RData")
  require(tidyr)
  gmtd = features
  
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  # browser()
  color.byNoQ <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  color.by <- enquo(color.byNoQ)
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  
  if(is.null(data.plot$feature.groups)) {
    data.plot$feature.groups = "1"
    feature.groups = data.plot$feature.groups
    names(feature.groups) = data.plot$feature.groups
  }
  
  s=AddModuleScore(
    object = object,
    features = gmtd,
    ctrl = min(100, data.plot$feature.groups %>% table() %>% max()),
    name = 'Combined'
  )
  ams = lapply(Idents(s) %>% levels(), FUN= function(id){
    lapply(seq(length(gmtd)), FUN = function(cIdx){
      subset(s, idents = id)[[paste0("Combined",cIdx)]][,1]%>% mean()
    })
  }) %>% unlist() %>% matrix(ncol=length(gmtd),byrow=T) %>% data.frame()
  if(is.list(gmtd))  colnames(ams) = names(gmtd)
  if(is.character(gmtd)) {
    colnames(ams) = gmtd
  }
  ams$id = Idents(s) %>% levels() 
  amsLong = gather(ams, feature.groups, avg.exp.scaled, 1:(ncol(ams)-1) )
  # colnames(amsLong) = c("feature.groups", "avg.exp.scaled")
  amsLong$features.plot = "XXXXXXXM-Score"
  
  meanAvgExpSca = data.plot[,c("feature.groups", "avg.exp", "id")] %>% tidyr::pivot_wider( names_from = c(feature.groups), values_from = avg.exp, values_fn = mean) %>% 
    tidyr::pivot_longer(values_to = "avg.exp", cols = names(gmtd), names_to = "feature.groups")
  
  amsLong = merge(meanAvgExpSca, amsLong, by=c("id", "feature.groups"))
  
  pct.expSca = data.plot[,c("feature.groups", "pct.exp", "id")] %>% tidyr::pivot_wider( names_from = c(feature.groups), values_from = pct.exp, values_fn = mean) %>% 
    tidyr::pivot_longer(values_to = "pct.exp", cols = names(gmtd), names_to = "feature.groups")
  
  amsLong = merge(pct.expSca, amsLong, by=c("id", "feature.groups"))
  # rbind(amsLong, data.plot)
  data = rbind(amsLong, data.plot)
  data$features.plot = factor(data$features.plot, levels = c(unique(features),"M-Score"))
  plot <- ggplot(data, 
                 mapping = aes(x = features.plot, y = id)) + 
    # color.byNoQ
  geom_point(mapping = aes(size = pct.exp,color = .data[[color.byNoQ]])) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  
  if (!is.null(x = feature.groups)) {
    # Warning: The `facets` argument of `facet_grid()` is deprecated as of ggplot2 2.2.0.
    # ℹ Please use the `rows` argument instead.
    # ℹ The deprecated feature was likely used in the Seurat package.
    # Please report the issue at <https://github.com/satijalab/seurat/issues>.
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + 
      theme(panel.spacing = unit(x = 1, units = "lines"), strip.background = element_blank())
  }
  # plot + scale_color_distiller(palette = cols)
  
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  
  labFun <- function(breakval){
    out = featureDat[breakval,"symbol"]
    out[is.na(out)] = "M-Score"
    return(out)
  }
  ylabFun <- function(val){stringr::str_replace(val,"SeuratProject_","bla")}
  plot= plot + scale_x_discrete(labels = labFun) +
    scale_y_discrete(labels = ylabFun) +
    ylab(clusters) +
    theme(axis.text.x = element_text(angle = 25, vjust = 0.5),
          strip.text.x = element_text(angle=5))
  # ggplotly(plot)
  return(plot)
}

