#' @title GaussSmoothArray
#' @description An internal function named GaussSmoothArray. Original from AnalyzeFMRI package
#' @param x The array to be smoothed.
#' @param voxdim The dimensions of the volume elements (voxel) that make up the array.
#' @param ksize The dimensions (in number of voxels) of the 3D discrete smoothing kernel used to smooth the array.
#' @param sigma The covariance matrix of the 3D Gaussian smoothing kernel. This matrix doesn't have to be non-singular; zero on the diagonal of sigma indicate no smoothing in that direction.
#' @param mask A 3D 0-1 mask that delimits where the smoothing occurs.
#' @param var.norm Logical flag indicating whether to normalize the variance of the smoothed array.
#' 
#' @export


GaussSmoothArray <- function(x, voxdim = c(1, 1, 1), ksize = 5, sigma = diag(3, 3), mask = NULL, var.norm= FALSE )
{
  filtmat <- GaussSmoothKernel(voxdim, ksize, sigma)
  
  if(!is.array(x)) return("x should be an array")
  if(length(dim(x)) != 3 && length(dim(x)) != 4) return("array x should be 3D or 4D")
  tmp <- FALSE
  if(length(dim(x)) == 3) {
    x <- array(x, dim = c(dim(x), 1))
    tmp <- TRUE
  }
  if(is.null(mask)) mask <- array(1, dim = dim(x)[1:3])
  
  if(var.norm){
    d <- .Fortran("gaussfilter2",
                  as.double(x),
                  as.integer(dim(x)[1]),
                  as.integer(dim(x)[2]),
                  as.integer(dim(x)[3]),
                  as.integer(dim(x)[4]),
                  as.double(filtmat),
                  as.integer(ksize),
                  as.double(mask),
                  double(length(x)),
                  PACKAGE = "TCIU")
    c1 <- array(d[[9]], dim = dim(x))
    if(tmp) c1 <- c1[, , , 1]
  }
  
  else
  {   d <- .Fortran("gaussfilter1",
                    as.double(x),
                    as.integer(dim(x)[1]),
                    as.integer(dim(x)[2]),
                    as.integer(dim(x)[3]),
                    as.integer(dim(x)[4]),
                    as.double(filtmat),
                    as.integer(ksize),
                    as.double(mask),
                    as.double(x),
                    PACKAGE = "TCIU")
  c1 <- array(d[[9]], dim = dim(x))
  if(tmp) c1 <- c1[, , , 1]
  }
  
  return(c1)
}
