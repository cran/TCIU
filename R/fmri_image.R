#' @title interactive graph object of the fMRI image
#' @description fMRI image visualization method, based on package \code{plotly}.
#'
#' @param fmridata a 4D array contains information for the fMRI spacetime image. The data should only contain the magnitude for the fMRI image.
#' @param option The default is 'manually'. If choose 'auto', then this function will lead you to key in the space (x,y,z) parameters and time (time) parameter for this function to generate graphs.
#' @param voxel_location a 3D array indicating the spatial location of the brain. If option is auto, set the voxel_location as NULL.
#' @param time time location for the voxel
#'
#' @details
#' The function \code{fmri_image} is used to create images for front view, side view, and top view of the fMRI image.
#' When providing the 4D array of the fMRI spacetime image and input the x,y,z position of the voxel, 
#' three views of the fMRI image and the time series image of the voxel will be shown.
#'
#' @author SOCR team <\url{http://socr.umich.edu/people/}>
#'
#' @return an interactive graph object of the fMRI image created by \code{plotly}
#' @export
#'
#' @import plotly
#' @importFrom forecast auto.arima
#' @importFrom ggplot2 ggplot geom_hline geom_vline xlim ylim
#' 
#' @examples
#' fmri_generate = fmri_simulate_func(dim_data = c(64, 64, 40), mask = mask)
#' fmri_image(fmri_generate$fmri_data, option='manually', voxel_location = c(40,22,33), time = 4)
#' 

fmri_image =  function(fmridata,
                       option = "manually",
                       voxel_location = NULL,
                       time = NULL) {
  xdim=dim(fmridata)[1]
  ydim=dim(fmridata)[2]
  zdim=dim(fmridata)[3]
  tdim=dim(fmridata)[4]
  if (option == "auto") {
    message("Please key in x, y, z, and time in sequence")
    x = 0
    y = 0
    z = 0
    t = 0
    while (x > xdim | x < 1 | x%%1 != 0) {
      x <- scan(n = 1)
      if (x > xdim | x < 1 | x%%1 != 0) {
        message("input x is not in range or is not an integer! Please retype!")
      }
    }
    while (y > ydim | y < 1 | y%%1 != 0) {
      y <- scan(n = 1)
      if (y > ydim | y < 1 | y%%1 != 0) {
        message("input y is not in range or is not an integer! Please retype!")
      }
    }
    while (z > zdim | z < 1 | z%%1 != 0) {
      z <- scan(n = 1)
      if (z > zdim | z < 1 | z%%1 != 0) {
        message("input z is not in range or is not an integer! Please retype!")
      }
    }
    while (t > tdim | t < 1 | t%%1 != 0) {
      t <- scan(n = 1)
      if (t > tdim | t < 1 | t%%1 != 0) {
        message("input time is not in range or is not an integer! Please retype!")
      }
    }
  } else if (option == "manually") {
    x = voxel_location[1]
    y = voxel_location[2]
    z = voxel_location[3]
    t = time
  }
  # time series
  try1 <- fmridata[x, y, z, ]
  try1 <- detrend(try1, bp = seq(21, tdim, by = 20))
  tstry1 <- ts(try1)
  trymod1 <- auto.arima(tstry1)
  ksmth <- ksmooth(c(1:tdim), tstry1, kernel = "normal", bandwidth = 5)
  smth <- smooth(tstry1)
  ksth <- data.frame(tstry1, smth, ksmth$y)
  TScore <- GTSplot(ksth, ts_name = c("original", "smooth", "ksmooth"), COLO = c("3399FF", "66FF33", "FF6666"))
  
  # z slice
  
  zfMRI <- t(fmridata[, , z, t])
  zslice <- ggplotly(ggplot() + geom_hline(yintercept = y, color = "red") + geom_vline(xintercept = x, color = "red") + xlim(0, 
                                                                                                                             xdim) + ylim(0, ydim)) %>% add_contour(z = ~zfMRI)
  
  # x slice
  
  xfMRI <- t(fmridata[x, , , t])
  xslice <- ggplotly(ggplot() + geom_hline(yintercept = z, color = "red") + geom_vline(xintercept = y, color = "red") + xlim(0, 
                                                                                                                             ydim) + ylim(0, zdim)) %>% add_contour(z = ~xfMRI)
  
  # y slice
  yfMRI <- t(fmridata[, y, , t])
  yslice <- ggplotly(ggplot() + geom_hline(yintercept = z, color = "red") + geom_vline(xintercept = x, color = "red") + xlim(0, 
                                                                                                                             xdim) + ylim(0, zdim)) %>% add_contour(z = ~yfMRI)
  
  # combine all 4 plots:
  result <- subplot(TScore, zslice, xslice, yslice, nrows = 2)
  return(result)
}