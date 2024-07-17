#' Plot log info detailed
#' @description plot output of optim 
#' @usage plotOptimInfo(svcov.object, ...)
#' @param svcov.object an optimized svcov.object
#' @param ... arguments for par()
#' @returns Outputs a sequence of plots detailing parameters during the 
#' optimization routine
#' @author Federico Blasi
plotOptimInfo <- function(svcov.object, ...){
  
  if(length(svcov.object@output) == 0){stop('did not find an output to work with.')}
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  index_pars <- base::grepl('par', base::colnames(svcov.object@output$loginfo))
  index_grad <- base::grepl('gr', base::colnames(svcov.object@output$loginfo))
  x_grid <- 1:base::dim(svcov.object@output$loginfo)[1]
  
  #graphics::par(mfrow = grDevices::n2mfrow(prod(dim(svcov.object@output$loginfo))/256), 
  #              oma = c(3, 2, 2, 2), mar = c(2.5,1,1,1))
  
  graphics::par(mfrow=c(1,2), ...)
  
  tmp_dm <- svcov::getDesignMatrix(svcov.object@model.list, data = svcov.object@data)
  
  tmp_names <- unlist(tmp_dm$par.pos[unlist(lapply(tmp_dm$par.pos, FUN = is.logical))])
  
  names_to_display <- rep(colnames(getDesignMatrix(svcov.object@model.list,
                                                   data = svcov.object@data)$model.matrix),6)[
                                                     tmp_names]
  
  tmp_names <- names(tmp_names[which(tmp_names == TRUE)])
  
  for(ii in 1:base::length(svcov.object@info$boundaries$theta_init)){
    graphics::plot(x = x_grid, y = svcov.object@output$loginfo[, index_pars][, ii], type = 'l', 
                   main = base::paste0(tmp_names[ii]," ",names_to_display[ii]),
                   xlab = "Iters", ylab = " ")
    graphics::abline(h = svcov.object@info$boundaries$theta_lower[ii], col = 'blue')
    graphics::abline(h = svcov.object@info$boundaries$theta_upper[ii], col = 'red')
    
    graphics::plot(x = x_grid, y = svcov.object@output$loginfo[,index_grad][,ii], type = 'l',
                   main = 'gradient', xlab = "Iters", ylab = " ")
  }
}

#' Plot variable
#' @description plot output of optim 
#' @usage plotOptimInfo(svcov.object, ...)
#' @param svcov.object an optimized svcov.object
#' @param ... arguments for par()
#' @returns Outputs a sequence of plots detailing parameters during the 
#' optimization routine
#' @author Federico Blasi
plotVar <- function(locs, vartoplot, ...){
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  par(mfrow = c(1, 1), oma = master_oma, mar = master_mars)
  # smooth
  layout(matrix(c(1,1,1,1,1,1,1,1,1,1,2), nrow = 1))
  
  plot(x@locs[,1], x@locs[,2],
       col = tim.colors(128)[cut(vartoplot, 
                                 breaks = seq(min(vartoplot),max(vartoplot),length.out=128), labels = F)], 
       cex=0.3,xlab = 'Longitude', ylab = 'Latitude', pch = 20, cex.axis = 1.5, cex.lab = 1.5)
  map(add=TRUE, resolution = 0)
  par(mar = c(4.5, 0, 0.5, 2.5))
  
  tmp_z <- matrix(1:100, nrow =1)
  tmp_x <- 1
  tmp_y <- seq(min(vartoplot), max(vartoplot), len = 100) 
  image(tmp_x, tmp_y, tmp_z, col = tim.colors(128), axes = FALSE, xlab = "", ylab = "")
  axis(4,lwd = 2,cex.axis = 1.5)

}



