#' Plot log info detailed
#' @description plot output of optim 
#' @usage plotOptimInfo(coco.object, ...)
#' @param coco.object an optimized coco.object
#' @param ... arguments for par()
#' @returns Outputs a sequence of plots detailing parameters during the 
#' optimization routine
#' @author Federico Blasi
#' @seealso [cocoOptim()]
plotOptimInfo <- function(coco.object, ...){
  
  if(length(coco.object@output) == 0){stop('did not find an output to work with.')}
  
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  
  index_pars <- base::grepl('par', base::colnames(coco.object@output$loginfo))
  index_grad <- base::grepl('gr', base::colnames(coco.object@output$loginfo))
  x_grid <- 1:base::dim(coco.object@output$loginfo)[1]

  graphics::par(mfrow=c(1,2), ...)
  
  tmp_dm <- cocons::getDesignMatrix(coco.object@model.list, data = coco.object@data)
  
  tmp_names <- unlist(tmp_dm$par.pos[unlist(lapply(tmp_dm$par.pos, FUN = is.logical))])
  
  names_to_display <- rep(colnames(getDesignMatrix(coco.object@model.list,
                                                   data = coco.object@data)$model.matrix), 6)[tmp_names]
  
  tmp_names <- names(tmp_names[which(tmp_names == TRUE)])
  
  if(coco.object@info$optim.type == "reml" ||
     coco.object@info$optim.type == "pml"){
    
    shift_mean <- length(base::attr(stats::terms(coco.object@model.list$mean), "term.labels")) + 
      base::attr(stats::terms(coco.object@model.list$mean), "intercept")
    
  } else{ shift_mean <- 0}

  for(ii in 1:sum(index_pars)){
    graphics::plot(x = x_grid, y = coco.object@output$loginfo[, index_pars][, ii], type = 'l', 
                   main = base::paste0(tmp_names[ii + shift_mean]," ", names_to_display[ii + shift_mean]),
                   xlab = "Iters", ylab = " ")
    graphics::abline(h = coco.object@info$boundaries$theta_lower[ii + shift_mean], col = 'blue')
    graphics::abline(h = coco.object@info$boundaries$theta_upper[ii + shift_mean], col = 'red')
    
    graphics::plot(x = x_grid, y = coco.object@output$loginfo[,index_grad][,ii], type = 'l',
                   main = 'gradient', xlab = "Iters", ylab = " ")
  }
}
