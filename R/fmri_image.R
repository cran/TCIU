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

fmri_image = function(fmridata,
                      option = "manually",
                      voxel_location = NULL,
                      time = NULL) {
  
  # Check input dimensions
  if (length(dim(fmridata)) != 4) {
    stop("fmridata must be a 4D array")
  }
  
  xdim = dim(fmridata)[1]
  ydim = dim(fmridata)[2]
  zdim = dim(fmridata)[3]
  tdim = dim(fmridata)[4]
  
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
    # Check manual inputs
    if (is.null(voxel_location) || length(voxel_location) != 3) {
      stop("voxel_location must be a vector of length 3 (x, y, z)")
    }
    if (is.null(time) || length(time) != 1) {
      stop("time must be a single integer")
    }
    
    x = voxel_location[1]
    y = voxel_location[2]
    z = voxel_location[3]
    t = time
    
    # Validate manual inputs
    if (x < 1 || x > xdim || x%%1 != 0) {
      stop(sprintf("x must be an integer between 1 and %d", xdim))
    }
    if (y < 1 || y > ydim || y%%1 != 0) {
      stop(sprintf("y must be an integer between 1 and %d", ydim))
    }
    if (z < 1 || z > zdim || z%%1 != 0) {
      stop(sprintf("z must be an integer between 1 and %d", zdim))
    }
    if (t < 1 || t > tdim || t%%1 != 0) {
      stop(sprintf("time must be an integer between 1 and %d", tdim))
    }
  }
  
  # Time series data processing
  try1 <- fmridata[x, y, z, ]
  
  # Detrend if pracma is available
  if (requireNamespace("pracma", quietly = TRUE)) {
    try1 <- pracma::detrend(try1, bp = seq(21, tdim, by = 20))
  } else {
    warning("Package 'pracma' not available. Skipping detrending.")
  }
  
  tstry1 <- ts(try1)
  
  # Smoothing
  if (requireNamespace("stats", quietly = TRUE)) {
    ksmth <- stats::ksmooth(1:tdim, tstry1, kernel = "normal", bandwidth = 5)
    smth <- stats::smooth(tstry1)
  } else {
    stop("stats package is required for smoothing")
  }
  
  # Create time series plot directly with plotly (replacing GTSplot)
  time_df <- data.frame(
    Time = 1:tdim,
    Original = as.numeric(tstry1),
    Smooth = as.numeric(smth),
    KSmooth = ksmth$y
  )
  
  # Create interactive time series plot
  TScore <- plot_ly(time_df) %>%
    add_trace(x = ~Time, y = ~Original, type = 'scatter', mode = 'lines',
              name = 'Original', line = list(color = '#3399FF')) %>%
    add_trace(x = ~Time, y = ~Smooth, type = 'scatter', mode = 'lines',
              name = 'Smooth', line = list(color = '#66FF33')) %>%
    add_trace(x = ~Time, y = ~KSmooth, type = 'scatter', mode = 'lines',
              name = 'KSmooth', line = list(color = '#FF6666')) %>%
    layout(
      title = paste0("Voxel (", x, ",", y, ",", z, ") Time Series"),
      xaxis = list(title = "Time"),
      yaxis = list(title = "Signal"),
      hovermode = 'x unified'
    )
  
  # Create slice plots
  # z slice (axial)
  zfMRI <- t(fmridata[, , z, t])
  zslice <- plot_ly(z = ~zfMRI, type = "contour") %>%
    layout(
      title = paste0("Z = ", z, " (Axial View)"),
      xaxis = list(title = "X", range = c(0, xdim)),
      yaxis = list(title = "Y", range = c(0, ydim)),
      shapes = list(
        list(
          type = "line",
          x0 = x, x1 = x, y0 = 0, y1 = ydim,
          line = list(color = "red", width = 2)
        ),
        list(
          type = "line",
          x0 = 0, x1 = xdim, y0 = y, y1 = y,
          line = list(color = "red", width = 2)
        )
      )
    )
  
  # x slice (sagittal)
  xfMRI <- t(fmridata[x, , , t])
  xslice <- plot_ly(z = ~xfMRI, type = "contour") %>%
    layout(
      title = paste0("X = ", x, " (Sagittal View)"),
      xaxis = list(title = "Y", range = c(0, ydim)),
      yaxis = list(title = "Z", range = c(0, zdim)),
      shapes = list(
        list(
          type = "line",
          x0 = y, x1 = y, y0 = 0, y1 = zdim,
          line = list(color = "red", width = 2)
        ),
        list(
          type = "line",
          x0 = 0, x1 = ydim, y0 = z, y1 = z,
          line = list(color = "red", width = 2)
        )
      )
    )
  
  # y slice (coronal)
  yfMRI <- t(fmridata[, y, , t])
  yslice <- plot_ly(z = ~yfMRI, type = "contour") %>%
    layout(
      title = paste0("Y = ", y, " (Coronal View)"),
      xaxis = list(title = "X", range = c(0, xdim)),
      yaxis = list(title = "Z", range = c(0, zdim)),
      shapes = list(
        list(
          type = "line",
          x0 = x, x1 = x, y0 = 0, y1 = zdim,
          line = list(color = "red", width = 2)
        ),
        list(
          type = "line",
          x0 = 0, x1 = xdim, y0 = z, y1 = z,
          line = list(color = "red", width = 2)
        )
      )
    )
  
  # Combine all 4 plots
  result <- subplot(
    TScore, 
    zslice, 
    xslice, 
    yslice,
    nrows = 2,
    titleY = TRUE,
    titleX = TRUE
  ) %>%
    layout(
      title = paste0("fMRI Visualization - Voxel (", x, ",", y, ",", z, "), Time = ", t),
      showlegend = TRUE
    )
  
  return(result)
}