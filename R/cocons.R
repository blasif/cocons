#' coco Class
#' @description Creates a coco object
#' @usage coco(type, data, locs, z, model.list, info, output = list())
#' @param type One of two available types "dense" or "sparse". See description.
#' @param data A data.frame with covariates information, where colnames(data) matches model.list specification
#' @param locs a matrix with locs matching data
#' @param z A vector with response values
#' @param model.list A list specyfing a model for each aspect of the spatial structure.
#' @param info A list specifying characteristics of the coco object
#' @param output empty or an output from optimparallel output, including as well boundaries 
#' @returns Creates a coco object
#' @author Federico Blasi
#' @seealso [spam::cov.wend1()]
#' 
coco <- function(type,
                  data,
                  locs,
                  z,
                  model.list,
                  info = list(),
                  output = list()) {
  
  .coco.check.type(type)
  .coco.check.data(data)
  .coco.check.locs(locs)
  
  if(is.vector(z)){
    z <- matrix(z, ncol = 1)
  }
  
  .coco.check.z(z)
  
  .coco.check.model.list(model.list,
                          data)
  .coco.check.info(type = type, 
                    info = info,
                    model.list = model.list)
  
  .coco.check.output(output)
  
  if(is.null(info$lambda)){
    info$lambda <- 0
  }
  
  if(T){
    if(!is.formula(model.list$smooth)){
      info$smooth_limits <- c(model.list$smooth[1],
                              model.list$smooth[1])
    }
  }
  
  methods::new(
    Class = "coco",
    type = type,
    data = data,
    locs = locs,
    z = z,
    model.list = model.list,
    info = info,
    output = output
  )
  
}
