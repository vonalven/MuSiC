# Utility Function
#
#
# Author: Jiaxin Fan, Xuran Wang
# Adapted by: Giacomo von Alvensleben, EPFL
###################################################
#' @title MuSiC2_Deconvolution
#'
#' @description This function is used to deconvolve bulk RNA-seq data using single-cell reference generated under a different condition.
#' @param bulk.control.mtx Matrix of expression for bulk data, control group
#' @param bulk.case.mtx Matrix of expression for bulk data, case group
#' @param sc.sce SingleCellExperiment for single cell data
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character, the phenoData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by the single cell dataset;
#' @param method can be 'toast' or 't_stats'. See MuSiC2 github documentation for more info;
#' @param expr_low numeric, cutoff for defining lowly expressed genes in bulk data. Genes with mean expression across samples in bulk data < expr_low will be excluded from cell-type-specific DE gene detection. Default is 20;
#' @param prop_r numeric, cutoff for defining rare cell types. Cell types with mean proportion across samples in bulk data < prop_r will be characterized as rare cell types. Otherwise, will be characterized as common cell types. Default is 0.1;
#' @param eps_c numeric, convergence cutoff for common cell types. The cell type proportion estimate is converged if absolute relative change of proportion estimates for the current iteration against the previous iteration < eps_c. Default is 0.05;
#' @param eps_r numeric, convergence cutoff for rare cell types. The cell type proportion estimate is converged if absolute change of proportion estimates for the current iteration against the previous iteration < eps_r. Default is 0.01;
#' @param n_resample numeric, number of resamples used for detecting cell-type-specific DE genes. Only relevant for t_stats mode. Default is 40;
#' @param sample_prop numeric, proportion of samples to be randomly sampled without replacement under each condition at each resampling. Only relevant for t_stats mode. Default is 0.5;
#' @param cutoff_expr numeric, cutoff for defining lowly expressed genes over resamples. Genes with average cell-type-specific expression calculated over all resamples in the lower cutoff_expr quantile are excluded from cell-type-specific DE gene detection. Only relevant for toast mode. Default is 0.05;
#' @param cutoff_c numeric, cutoff for defining cell-type-specific DE genes for common cell types. Genes with the value of statistic, defined as the absolute value of the ratio of the mean and standard deviation of the log fold change over all resamples, in the upper cutoff_c quantile are considered as cell-type-specific DE genes. Default is 0.05;
#' @param cutoff_r numeric, cutoff for defining cell-type-specific DE genes for rare cell types. Genes with the value of statistic, defined as the absolute value of the ratio of the mean and standard deviation of the log fold change over all resamples, in the upper cutoff_r quantile are considered as cell-type-specific DE genes. Default is 0.01;
#' @param maxiter numeric, maximum number of iterations. Default is 200;
#' @param markers vector or list of gene names. Default as NULL, i.e., use all genes that provided by both bulk and single cell datasets;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param cap numeric, cutoff on maximum number of genes removed for each cell type. For each cell type, at the maximum, genes with FDR adjusted p-value within the lower cap quantile are removed. Default is 0.3;
#' @param centered logic, subtract avg of Y and D;
#' @param normalize logic, divide Y and D by their standard deviation;
#' @param nb_cores number of cores, for parallel computation. Only relevant for t_stats mode;
#' @return If MuSiC2 converges, return:
#' \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {n.iter: numeric, number of iterations.}
#'    \item {DE.genes: vector, cell-type-specific DE genes being removed.}
#'    }
#'  Or if MuSiC2 does not converge, return:
#'  \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {id.not.converge: vector, sample ids that failed to converge.}
#'    }
#' @details
#' As detailed in the MuSiC2 github documentation, MuSiC2 can be run in 2 modes to identify cell-type-specific DE genes:
#' \itemize{
#'    \item {t-statistics: the original method described in the manuscript.}
#'    \item {TOAST (Li and Wu (2019)) with P-value cutoffs for defining cell-type-specific DE genes. MuSiC2 with TOAST runs faster and converges faster than MuSiC2 with T statistics as it does not require resamplings. Therefore, based on our simulation and real data studies, we recommend using TOAST (Li and Wu (2019)) with P-value cutoffs for deconvolving bulk samples with healthy single cell reference, and ad hoc T statistics for deconvolving bulk samples with diseased single cell reference.}
#'    }
#' When the difference between the diseased and healthy bulk samples is small, or the sample sizes are small, TOAST (Li and Wu (2019)) may not be able to detect any cell-type-specific DE genes, and we recommend using MuSiC2 with T statistics under this case.
#' @author Jiaxin Fan, Xuran Wang
#' @author Adapted by: Giacomo von Alvensleben, EPFL
#' @seealso
#' \code{\link{music_prop}}
#' @export

music2_prop <- function(bulk.control.mtx,
                        bulk.case.mtx,
                        sc.sce,
                        clusters,
                        samples,
                        select.ct,
                        method      = "toast",
                        expr_low    = 20, 
                        prop_r      = 0.1, 
                        eps_c       = 0.05, 
                        eps_r       = 0.01,
                        n_resample  = 40, 
                        sample_prop = 0.5,
                        cutoff_expr = 0.05, 
                        cutoff_c    = 0.05, 
                        cutoff_r    = 0.01, 
                        maxiter     = 200, 
                        markers     = NULL, 
                        cell_size   = NULL, 
                        ct.cov      = FALSE,
                        nu          = 1e-04,
                        eps         = 0.01,
                        cap         = 0.3,
                        centered    = FALSE,
                        normalize   = FALSE,
                        # random_seed = 123,
                        nb_cores    = 1){
  
  # set.seed(random_seed)
  
  gene.bulk <- intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if(length(gene.bulk) < 0.1 * min(nrow(bulk.control.mtx), nrow(bulk.case.mtx))){
    stop("Not enough genes for bulk data! Please check gene annotations.")
  }
  gene_all <- intersect(gene.bulk, rownames(sc.sce))
  if(length(gene_all) < 0.2 * min(length(gene.bulk), nrow(sc.sce))){
    stop("Not enough genes between bulk and single-cell data! Please check gene annotations.")
  }
  if(!(method %in% c("toast", "t_stats"))){
    stop("method must be 'toast' or 't_stats'.")
  }
  if(is.null(select.ct)){
    select.ct <- unique(as.character(SingleCellExperiment::colData(sc.sce)[[clusters]]))
  }
  # remove spaces and non-alphanumeric characters from cell types to allow creating design matrices with TOAST
  ct.map <- data.frame(input_ct     = select.ct,
                       formatted_ct = gsub(" ", "_", gsub("[^[:alnum:] ]", "", select.ct)), 
                       stringsAsFactors = F)
  sc.coldata             <- SingleCellExperiment::colData(sc.sce)
  sc.coldata[[clusters]] <- plyr::mapvalues(sc.coldata[[clusters]], from = ct.map$input_ct, to = ct.map$formatted_ct)
  colData(sc.sce)        <- sc.coldata
  select.ct              <- unique(as.character(SingleCellExperiment::colData(sc.sce)[[clusters]]))

  if(!is.null(cell_size)){
    cell_size[[1]] <- plyr::mapvalues(as.character(cell_size[[1]]), from = ct.map$input_ct, to = ct.map$formatted_ct)
  }
  
  bulk.mtx    <- cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, ])
  bulk.mtx    <- bulk.mtx[gene_all, ]
  sc.iter.sce <- sc.sce[gene_all, ]
  
  # remove lowly expressed genes from DE analysis: i.e., gene with average expression < expr_low
  cat("> Removing lowly expressed genes (avg gene expression <", expr_low, ")...\n")
  expr      <- apply(bulk.mtx, 1, mean)
  exp_genel <- names(expr[expr >= expr_low])
  
  # Analyse separately based on their case/control status
  bulk.control <- bulk.mtx[, colnames(bulk.control.mtx)]
  bulk.case    <- bulk.mtx[, colnames(bulk.case.mtx)]
  
  # Step 1: cell type deconvolution, set initial value
  # estimate cell type proportion for controls using music.
  # tutorial for basic music see https://xuranw.github.io/MuSiC/articles/MuSiC.html
  cat("> Estimating cell type proportions from control samples using all scRNA-seq genes (N =", nrow(sc.sce), ") and bulk RNA-seq genes (N =", nrow(bulk.control), ")...\n")
  prop_control <- music_prop(bulk.mtx  = bulk.control, 
                             sc.sce    = sc.sce,
                             clusters  = clusters, 
                             samples   = samples, 
                             select.ct = select.ct,
                             markers   = markers, 
                             cell_size = cell_size, 
                             ct.cov    = ct.cov, 
                             iter.max  = maxiter,
                             nu        = nu, 
                             eps       = eps, 
                             centered  = centered, 
                             normalize = normalize, 
                             verbose   = F)
  prop_control <- prop_control$Est.prop.weighted
  
  cat("> Estimating cell type proportions from case samples using all scRNA-seq genes (N =", nrow(sc.sce), ") and bulk RNA-seq genes (N =", nrow(bulk.case.mtx), ")...\n")
  prop_case_fix <- NULL
  prop_case_ini <- music_prop(bulk.mtx  = bulk.case.mtx, 
                              sc.sce    = sc.sce,
                              clusters  = clusters, 
                              samples   = samples, 
                              select.ct = select.ct,
                              markers   = markers, 
                              cell_size = cell_size, 
                              ct.cov    = ct.cov, 
                              iter.max  = maxiter,
                              nu        = nu, 
                              eps       = eps, 
                              centered  = centered, 
                              normalize = normalize, 
                              verbose   = F)
  prop_case_ini <- prop_case_ini$Est.prop.weighted
  # do a copy
  prop_case_ini.copy <- prop_case_ini
  
  
  # prop_all <- music_prop(bulk.mtx  = cbind(bulk.control.mtx, bulk.case.mtx), 
  #                        sc.sce    = sc.sce,
  #                        clusters  = clusters, 
  #                        samples   = samples, 
  #                        select.ct = select.ct,
  #                        markers   = markers, 
  #                        cell_size = cell_size, 
  #                        ct.cov    = ct.cov, 
  #                        iter.max  = maxiter,
  #                        nu        = nu, 
  #                        eps       = eps, 
  #                        centered  = centered, 
  #                        normalize = normalize, 
  #                        verbose   = F)
  # prop_all <- prop_all$Est.prop.weighted
  # prop_all.2 <- rbind(prop_control, prop_case_ini)
  # identical(prop_all, prop_all.2)
  
  prop_CASE <- prop_case_ini
  prop_all  <- rbind(prop_control, prop_CASE)
  
  # start iteration
  cat("> Starting iterative analysis...\n")
  iter    <- 1
  ncell   <- length(select.ct)
  id_conv <- NULL
  
  if(method == "toast"){
    
    Pheno           <- data.frame(condition = factor(c(rep("control", ncol(bulk.control.mtx)), rep("case", ncol(bulk.case.mtx))), levels = c("control", "case")))
    rownames(Pheno) <- c(colnames(bulk.control.mtx), colnames(bulk.case.mtx))
    
    # step 2: identify cell-type-specific DE genes using TOAST
    while(iter <= maxiter){
      cat("## Iteration", iter, "\n")
      
      # first, log transformed the bulk expression
      Y_raw      <- log1p(bulk.mtx)
      design     <- Pheno
      Prop       <- prop_all[rownames(Pheno),]
      stopifnot(identical(rownames(design), rownames(Prop)))
      
      Design_out   <- makeDesign(design, Prop)
      fitted_model <- fitModel(Design_out, Y_raw)
      # run TOAST to detect DE between conditions
      res_table    <- csTest(fitted_model, coef = "condition", verbose = F)
      # lapply(res_table, function(x) base::range(x$fdr))
      
      # average cell type proportion
      mex <- apply(prop_all, 2, mean)
      lr  <- NULL
      for(celltype in select.ct){
        # celltype <- select.ct[1]
        m           <- mex[celltype]
        DE          <- res_table[[celltype]]
        pval        <- DE$fdr
        names(pval) <- rownames(DE)
        pval        <- pval[names(pval) %in% exp_genel]
        # select DE genes 
        if(m >= prop_r){
          lr <- c(lr, names(pval[pval <= cutoff_c & pval <= quantile(pval, prob = cap)]))
        }else{
          lr <- c(lr, names(pval[pval <= cutoff_r & pval <= quantile(pval, prob = cap)]))
        }
      }
      lr <- unique(lr)
      
      # step 3: update sc gene list
      # remove identified DE genes from sc rna-seq data
      cat("> ", length(lr), " DE genes identified...\n")
      l           <- setdiff(gene_all, lr)
      sc.iter.sce <- sc.sce[l, ]
      
      # step 1: update cell type proportion based on new gene list
      if(length(id_conv) > 0){
        case_sample <- bulk.case[ , !(colnames(bulk.case) %in% id_conv), ]
      } else{
        case_sample <- bulk.case
      }
      
      cat("> Estimating cell type proportions in case samples using fitered scRNA-seq genes (N =", nrow(sc.iter.sce), ") and bulk RNA-seq genes (N =", nrow(case_sample), ")...\n")
      prop_case <- music_prop(bulk.mtx  = case_sample, 
                              sc.sce    = sc.iter.sce,
                              clusters  = clusters, 
                              samples   = samples, 
                              select.ct = select.ct,
                              markers   = markers, 
                              cell_size = cell_size, 
                              ct.cov    = ct.cov, 
                              iter.max  = maxiter,
                              nu        = nu, 
                              eps       = eps, 
                              centered  = centered, 
                              normalize = normalize, 
                              verbose   = F)
      prop_case <- prop_case$Est.prop.weighted
      prop_CASE <- rbind(prop_case, prop_case_fix)
      
      if(length(id_conv) == 1){
        rownames(prop_CASE) <- c(rownames(prop_case), id_conv)
      }
      prop_all <- rbind(prop_control, prop_CASE)
      
      # check convergence, by cell type
      prop_case <- prop_case[rownames(prop_case_ini), ]
      pc        <- abs(prop_case - prop_case_ini)
      conv      <- pc
      conv[,]   <- 1
      # use difference if rare cell type
      conv[prop_case_ini <= prop_r] <- ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
      # use percent change if common cell type
      pc[prop_case_ini > prop_r]    <- pc[prop_case_ini > prop_r] / prop_case_ini[prop_case_ini > prop_r]
      conv[prop_case_ini > prop_r]  <- ifelse(pc[prop_case_ini > prop_r] < eps_c, 0, 1)
      convf                         <- apply(conv, 1, function(x) all(x==0))
      
      # if an id converged, not updating anymore
      all_converge  <- FALSE
      id_conv       <- c(id_conv, names(convf[convf == TRUE]))
      prop_case_ini <- prop_CASE[!rownames(prop_CASE) %in% id_conv,]
      prop_case_fix <- prop_CASE[rownames(prop_CASE) %in% id_conv,]
      
      # if all converged or if only one subjects not converging --> music2 converged
      if(is.vector(prop_case_ini)){
        all_converge <- TRUE
        break
      }else if(nrow(prop_case_ini) == 0){
        all_converge <- TRUE
        break
      }
      iter <- iter + 1
    }
    
    cat("## Done!")
    
    
    # quantify differences with initial estimate and print a message
    pad_strings <- function(strings, length){
      sapply(strings, function(x) {
        paste0(x, strrep(" ", length - nchar(x)))
      })
    }
    prop_all.control <- prop_all[colnames(bulk.control.mtx), ]
    prop_all.case    <- prop_all[colnames(bulk.case.mtx), ]
    stopifnot(identical(colnames(prop_all.control), colnames(prop_control)))
    stopifnot(identical(colnames(prop_all.case), colnames(prop_case_ini.copy)))
    stopifnot(identical(colnames(prop_all.control), colnames(prop_all.case)))
    
    avg_prop_control.diff <- (apply(prop_all.control, 2, mean) - apply(prop_control, 2, mean)) * 100
    avg_prop_case.diff    <- (apply(prop_all.case, 2, mean) - apply(prop_case_ini.copy, 2, mean)) * 100
    diff.df               <- data.frame(celltype          = names(avg_prop_control.diff),
                                        diff_perc_control = unname(avg_prop_control.diff),
                                        diff_perc_case    = unname(avg_prop_case.diff),
                                        stringsAsFactors = F)
    diff.df$original_celltype <- plyr::mapvalues(diff.df$celltype, from = ct.map$formatted_ct, to = ct.map$input_ct)
    
    diff.df$label <- paste0("\t\t(", 1:nrow(diff.df), ") ",
                            pad_strings(diff.df$original_celltype, max(nchar(diff.df$original_celltype))), 
                            "\t\tcontrol: ", ifelse(diff.df$diff_perc_control > 0, "+", ""), format(round(diff.df$diff_perc_control, 3), nsmall = 3), "% (expected 0%);",
                            "\t\tcase: ", ifelse(diff.df$diff_perc_case > 0, "+", ""), format(round(diff.df$diff_perc_case, 3), nsmall = 3), "%")
    cat("\n> The following changes have been detected compared to the initial estimates:\n", paste(diff.df$label, collapse = "\n"), "\n\n")
    cat("> Original set of cell types:\n", paste(paste0("\t\t", sort(ct.map$input_ct)), collapse = "\n"), "\n\n")
    
    colnames(prop_all) <- plyr::mapvalues(colnames(prop_all), from = ct.map$formatted_ct, to = ct.map$input_ct)
    
    if(all_converge){
      return(list("Est.prop"    = prop_all,
                  "convergence" = TRUE,
                  "n.iter"      = iter,
                  "DE.genes"    = lr))
    } else{
      return(list("Est.prop"        = prop_all,
                  "convergence"     = FALSE,
                  "id.not.converge" = rownames(prop_case_ini)))
    }
    
    
  } else if(method == "t_stats"){
    
    while(iter <= maxiter){
      cat("## Iteration", iter, "\n")
      
      # step 2: identify cell-type-specific DE genes
      # calculate mean/sd of log fold change as an indicator of differential expression
      
      cat("> Indentifying cell-type-specific DE genes...\n")
      stats.list <- parallel::mclapply(mc.cores = nb_cores, X = 1:n_resample, FUN = function(i){
        # i <- 1
        # set.seed(random_seed)
        id_h      <- sample(colnames(bulk.control), round(ncol(bulk.control) * sample_prop))
        control_s <- bulk.control[exp_genel, colnames(bulk.control) %in% id_h]
        prop_h    <- prop_control[colnames(control_s), ]
        
        mod0 <- apply(control_s, 1, function(x){
          mod <- nnls(prop_h, x)
          if(mod$mode == 1){
            return(mod$x)
          }else{
            return(rep(0, ncell))
          }
        })
        # MOD0 <- MOD0 + mod0
        
        # set.seed(random_seed)
        id_d   <- sample(colnames(bulk.case), round(ncol(bulk.case) * sample_prop))
        case_s <- bulk.case[exp_genel, colnames(bulk.case) %in% id_d]
        prop_d <- prop_CASE[colnames(case_s), ]
        
        mod1 <- apply(case_s, 1, function(x){
          mod <- nnls(prop_d,x)
          if(mod$mode == 1){
            return(mod$x)
          }else{
            return(rep(0, ncell))
          }
        })
        # MOD1  <- MOD1 + mod1
        # LOGFC <- rbind(LOGFC, log1p(mod1) - log1p(mod0))
        list(mod0 = mod0,
             mod1 = mod1,
             lfc  = log1p(mod1) - log1p(mod0))
      })
      MOD0  <- Reduce("+", lapply(stats.list, function(x) x$mod0))
      MOD1  <- Reduce("+", lapply(stats.list, function(x) x$mod1))
      LOGFC <- do.call(rbind, lapply(stats.list, function(x) x$lfc))
      
      rcv <- parallel::mclapply(mc.cores = nb_cores, X = 1:ncell, FUN = function(i){
        s   <- LOGFC[seq(from = i, to = nrow(LOGFC), by = ncell), ]
        apply(s, 2,function(x) ifelse(mean(x) == 0, 0, mean(x) / sd(x)))
      })
      rcv <- do.call(rbind, rcv)
      
      abs_rcv_logfc  <- abs(rcv)
      MOD0           <- MOD0 / n_resample
      MOD1           <- MOD1 / n_resample
      rownames(MOD0) <- rownames(MOD1) <- rownames(abs_rcv_logfc) <- select.ct
      
      # average cell type proportion
      mex <- apply(prop_all, 2, mean)
      lr  <- lapply(select.ct, function(celltype){
        m  <- mex[celltype]
        rh <- MOD0[celltype,]
        rd <- MOD1[celltype,]
        # for genes with average expression within lower cutoff_expr for both conditions are removed from cell-type-specific DE genes detection
        llr <- unique(intersect(names(rd[rd <= quantile(rd, prob = cutoff_expr)]),
                                names(rh[rh <= quantile(rh, prob = cutoff_expr)])))
        x <- abs_rcv_logfc[celltype, ]
        x <- x[!(names(x) %in% llr)]
        # select genes with large mean/cv of log fc as DE
        if(m >= prop_r){
          names(x[x >= quantile(x, prob = 1 - cutoff_c)])
        }else{
          names(x[x >= quantile(x, prob = 1 - cutoff_r)])
        }
      })
      lr <- unique(unlist(lr))
      
      # step 3: update sc gene list
      # remove identified DE genes from sc rna-seq data
      cat("> ", length(lr), " DE genes identified...\n")
      l           <- setdiff(gene_all, lr)
      sc.iter.sce <- sc.sce[l, ]
      
      # step 1: update cell type proportion based on new gene list
      if(length(id_conv) > 0){
        case_sample <- bulk.case[ ,(!colnames(bulk.case) %in% id_conv)]
      }else{
        case_sample <- bulk.case
      }
      
      cat("> Estimating cell type proportions in case samples using fitered scRNA-seq genes (N =", nrow(sc.iter.sce), ") and bulk RNA-seq genes (N =", nrow(case_sample), ")...\n")
      prop_case <- music_prop(bulk.mtx  = case_sample, 
                              sc.sce    = sc.iter.sce,
                              clusters  = clusters, 
                              samples   = samples, 
                              select.ct = select.ct,
                              markers   = markers, 
                              cell_size = cell_size, 
                              ct.cov    = ct.cov, 
                              iter.max  = maxiter,
                              nu        = nu, 
                              eps       = eps, 
                              centered  = centered, 
                              normalize = normalize, 
                              verbose   = F)
      prop_case <- prop_case$Est.prop.weighted
      prop_CASE <- rbind(prop_case, prop_case_fix)
      
      if(length(id_conv) == 1){
        rownames(prop_CASE) <- c(rownames(prop_case), id_conv)
      }
      prop_all <- rbind(prop_control,prop_CASE)
      
      # check convergence, by cell type
      prop_case <- prop_case[rownames(prop_case_ini),]
      pc        <- abs(prop_case-prop_case_ini)
      conv      <- pc
      conv[,]   <- 1
      # use difference if rare cell type
      conv[prop_case_ini <= prop_r] <- ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
      # use percent change if common cell type
      pc[prop_case_ini > prop_r]    <- pc[prop_case_ini > prop_r] / prop_case_ini[prop_case_ini > prop_r]
      conv[prop_case_ini > prop_r]  <- ifelse(pc[prop_case_ini > prop_r] < eps_c, 0, 1)
      
      convf  <- apply(conv, 1, function(x) all(x==0))
      
      # if an id converged, not updating anymore
      all_converge  <- FALSE
      id_conv       <- c(id_conv, names(convf[convf == TRUE]))
      prop_case_ini <- prop_CASE[!rownames(prop_CASE) %in% id_conv, ]
      prop_case_fix <- prop_CASE[rownames(prop_CASE) %in% id_conv, ]
      
      # if all converged or if only one subjects not converging--> music2 converged
      if(is.vector(prop_case_ini)){
        all_converge <- TRUE
        break
      }else if(nrow(prop_case_ini) == 0){
        all_converge <- TRUE
        break
      }
      iter <- iter + 1
    }
    
    cat("## Done!")
    
    # quantify differences with initial estimate and print a message
    pad_strings <- function(strings, length){
      sapply(strings, function(x) {
        paste0(x, strrep(" ", length - nchar(x)))
      })
    }
    prop_all.control <- prop_all[colnames(bulk.control.mtx), ]
    prop_all.case    <- prop_all[colnames(bulk.case.mtx), ]
    print(colnames(prop_all.control))
    print(colnames(prop_control))
    print(colnames(prop_all.case))
    print(colnames(prop_case_ini.copy))
    stopifnot(identical(colnames(prop_all.control), colnames(prop_control)))
    stopifnot(identical(colnames(prop_all.case), colnames(prop_case_ini.copy)))
    stopifnot(identical(colnames(prop_all.control), colnames(prop_all.case)))
    
    avg_prop_control.diff <- (apply(prop_all.control, 2, mean) - apply(prop_control, 2, mean)) * 100
    avg_prop_case.diff    <- (apply(prop_all.case, 2, mean) - apply(prop_case_ini.copy, 2, mean)) * 100
    diff.df               <- data.frame(celltype          = names(avg_prop_control.diff),
                                        diff_perc_control = unname(avg_prop_control.diff),
                                        diff_perc_case    = unname(avg_prop_case.diff),
                                        stringsAsFactors = F)
    diff.df$original_celltype <- plyr::mapvalues(diff.df$celltype, from = ct.map$formatted_ct, to = ct.map$input_ct)
    
    diff.df$label <- paste0("\t\t(", 1:nrow(diff.df), ") ",
                            pad_strings(diff.df$original_celltype, max(nchar(diff.df$original_celltype))), 
                            "\t\tcontrol: ", ifelse(diff.df$diff_perc_control > 0, "+", ""), format(round(diff.df$diff_perc_control, 3), nsmall = 3), "% (expected 0%);",
                            "\t\tcase: ", ifelse(diff.df$diff_perc_case > 0, "+", ""), format(round(diff.df$diff_perc_case, 3), nsmall = 3), "%")
    cat("\n> The following changes have been detected compared to the initial estimates:\n", paste(diff.df$label, collapse = "\n"), "\n\n")
    cat("> Original set of cell types:\n", paste(paste0("\t\t", sort(ct.map$input_ct)), collapse = "\n"), "\n\n")
    
    
    colnames(prop_all) <- plyr::mapvalues(colnames(prop_all), from = ct.map$formatted_ct, to = ct.map$input_ct)
    
    if(all_converge){
      return(list("Est.prop"    = prop_all,
                  "convergence" = TRUE,
                  "n.iter"      = iter,
                  "DE.genes"    = lr))
    } else{
      return(list("Est.prop"        = prop_all,
                  "convergence"     = FALSE,
                  "id.not.converge" = rownames(prop_case_ini)))
    }
    
  }
  
}



