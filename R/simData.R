#' Simulate spatially autocorrelated (SAC) data
#'
#' Use this function to simulate data on three different landscapes; a Gaussian, Bernoulli or zero-inflated Poisson distributed repsonse variable; and four different causes of SAC or reference data, i.e. no SAC. 
#' @usage simData(dataset,
#'         filename = "default",
#'         gridsize = c(50L, 50L),
#'         cvfold = 10L,
#'         cvblock.size = c(10,10),
#'         r.seed = 20151126,
#'         n.predictor = 7L,
#'         f.smooth = list(function() lon,
#'                         function() lat,
#'                         function() (lon - mean(lon))^2,
#'                         function() (lat - mean(lat))^2,
#'                         function() x3^x4 * x4^x3,
#'                         function() x1^x1 * x3^x4,
#'                         function() x2^x1 * x4^x3 * log(x5 + 1)),
#'         f.realistic = list(var = 0.1, scale = 0.1),
#'         f.real = list(resolution = 10L,
#'                       extent = c("5N24E", "7S37E"),
#'                       bio.vars = c("bio1", "bio19", "bio2", "bio12", "bio4", "bio18", "bio3")),
#'         f.response = c("x1", "x4", "x4^2", "x3*x4", "x3"),
#'         par.response = "default", 
#'         f.sac1 = list(corCoef = -3,
#'                       sarFactor = 1),
#'         f.sac2 = "x1", 
#'         f.sac3 = c("^","*"), 
#'         f.sac4 = list(dispersal.max = 0.1, 
#'                       dispersal.shape = 30),
#'         interactive = TRUE)
#' @param dataset Input character of the form \code{"abc"}, with:
#' \describe{
#'    \item{\code{a}}{predictor landscape: \describe{
#'      \item{1}{smooth (linear and non-linear gradients without noise)}
#'      \item{2}{realistic (unconditional Gaussian random fields from exponential covariance models)}
#'      \item{3}{real (Real bio-climatic predictors from \url{http://www.worldclim.org})}
#'      }}
#'    \item{\code{b}}{distribution of the response variable: \describe{
#'      \item{1}{Gaussian}
#'      \item{2}{Bernoulli}
#'      \item{3}{zero-inflated Poisson}
#'      }}
#'    \item{\code{c}}{SAC scenario: \describe{
#'      \item{0}{Reference, i.e. no SAC}
#'      \item{1}{SAC onto response variable}
#'      \item{2}{Omitted predictor}
#'      \item{3}{Wrong functional form, e.g. intentionally miss quadratic term or interaction}
#'      \item{4}{Dispersal}
#'      }}
#' }
#' @param filename The destination file name (character). Defaults to "dataset\code{dataset}", e.g. "dataset110".
#' @param gridsize Vector defining [1] the number of cells in x direction (Longitude), and [2] the number of cells in y direction (Latitude).
#' @param cvfold Number of unique Cross-Validation (CV) IDs to be assigned blockwise to the data. 
#' @param cvblock.size Number of cells in x, y direction in one CV block.
#' @param r.seed Randomisation value to be used in \code{\link[base]{set.seed}} before any stochastic process. Defaults to \code{20151126}
#' @param n.predictor Number of predictors to be simulated. 
#' @param f.smooth If \code{dataset = "1**"}: List of \code{n.predictor} linear and non-linear functions. Can be functions of Longitude (\code{lon}) and/or Latitude (\code{lat}).
#' @param f.realistic If \code{dataset = "2**"}: A list of \code{var} and \code{scale}, which are passed to \code{\link[RandomFields]{RMexp}} to compute the exponential covariance model. If both arguments are of length \code{n.predictor} a new model is computed for every predictor. 
#' @param f.real If \code{dataset = "3**"}: A list comprising: 
#' \itemize{
#'  \item \code{resolution} [minutes of a degree] = 2.5, 5, and 10 (default). Defines the resoultion of the global interpolated climate data from \url{http://www.worldclim.org}; 
#'  \item \code{extent} = numeric vector of two geogr. coordinates (diagonal corners).
#'  \item \code{bio.vars} = character string of length \code{n.predictor} defining which bioclimatic should variables be used.
#'}
#' @param f.response Character string of mathematical terms (based on predictors x1, x2,...) yielding the response varibale.
#' @param par.response Coefficients for the elements in f.response. By default this numeric vector contains an intercept (first element) and beta values for every element in f.response. If the distribution is set to Gaussian, an addtional (last) element is provided to set the standard deviation in \code{\link[stats]{rnorm}}. 
#' Defaults to 
#' \itemize{
#'  \item Gaussian: 0.8 (intercept), 0.2, 0.9, -0.8, -0.6, 0.5, 0.2 (Gaussian error)
#'  \item Bernoulli: 0.2 (intercept), 4.5, -1.2, -1.2, -1.1, 0.9
#'  \item Poisson: 0.2 (intercept), 1.6, 0.9, 0.8, -0.8, 0.5  
#'}
#' \strong{Please note:} Poisson is zero-inflated and therefore requires a list of two numeric parameter vectors. First item is a numeric vector setting the Bernoulli coefficients, second the Poisson coefficients. 
#' @param f.sac1 If \code{dataset = "**1"}: List of \code{corCoef}, a coefficient impacting the correlation structure, and (only if \code{dataset = "*11"}) \code{sarFactor}, a factor determining the magnitude of SAC added to the existing response varibale.
#' @param f.sac2 If \code{dataset = "**2"}: Name of the predictor(s) to be omitted in the model structure.
#' @param f.sac3 If \code{dataset = "**3"}: Character string of "^" and/or "*" to filter (\code{\link[base]{grep}}) and omit respective terms in the model structure. 
#' @param f.sac4 If \code{dataset = "**4"}: A list of \code{dispersal.max} = maximum dispersal factor, and \code{dispersal.shape} = shape factor, the higher the more skewed the exponential curve.
#' @param interactive Defaults to \code{TRUE}. If \code{FALSE} existing files are overwritten and data (if not existiing) downloaded automatically, i.e. without asking.
#' @return No R output. Data and instruction is saved to netCDF file (\code{filename.nc}).
#' @details Designed to run in default mode. 
#'          
#'          The cell values in the case of the \strong{\code{realistic} predictor landscape} are simulated (with \code{\link[RandomFields]{RFsimulate}}) based on a stationary isotropic covariance model. Because the covariance function (\eqn{C(r)=e^{-r}}) largely (by default: only) depends on the distance between two points the resulting predictor landscape becomes spatially autocorrelated. 
#'          (Please note: We modify the convariance function by altering \code{var} (additional variance argument), and \code{scale} (additional scale argument) in \code{f.realistic}.)
#' @seealso \code{\link{extract.ncdf}} which allows you to readily extract data and attributes from the netCDF file.
#'
#' @author Severin Hauenstein <severin.hauenstein@biom.uni-freiburg.de>
#'
#' @note \itemize{
#' \item SAC cause 1, i.e. spatial error onto response (\code{dataset = "**1"}), is computationally burdensome because of the inversion of the distance matrix. This can be severe for large grids, i.e. \code{gridsize > c(100,100)}.
#' \item It may be necessary to change the parameter settings (\code{par.response}) for SAC cause 3. Particularly in the case of a Bernoully distributed response variable the magnitude of SAC is rather sensitive to the \code{par.response} values.
#' }
#' @import ncdf4 RandomFields raster wordspace
#' @export
#' @examples
#' \dontrun{
#'
#' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# 
#' # example structure                           #
#' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#' # landscape, distribution, SAC cause
#' 
#' # simulate data with simData
#' 
#' # extract data with extract.ncdf
#' # lattice::levelplot of the response variable
#' 
#' # build linear model: model structure in attributes of netCDF file
#' 
#' # compute residulas
#' # Uni- and multivariate spatial correlograms with ncf::correlog
#' # plot correlogram to check spatial autocorrelation
#'
#' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#' # 12 datasets with different settings         #
#' #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#' 
#' library(lattice)
#' library(ncf)
#' #---------------------------------------------#
#' # smooth landscape, Gaussian distribution, refrence data
#' 
#' simData("110")
#' 
#' d110 <- extract.ncdf("dataset110.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d110) # levelplot response
#' 
#' fm110 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d110, family = "gaussian")
#' summary(fm110) 
#' res110 <- residuals(fm110) # calculate residuals
#' co110 <- correlog(d110$Lat, d110$Lon, res110, increment=0.02, resamp=1) # check autocorrleation
#' plot(co110$mean.of.class, co110$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset110")
#' 
#' #---------------------------------------------#
#' # smooth landscape, Gaussian distribution, SAC onto response
#' 
#' simData("111")
#' 
#' d111 <- extract.ncdf("dataset111.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d111) # levelplot response
#' 
#' fm111 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d111, family = "gaussian")   
#' summary(fm111) 
#' res111 <- residuals(fm111) # calculate residuals
#' co111 <- correlog(d111$Lat, d111$Lon, res111, increment=0.02, resamp=1) # check autocorrleation
#' plot(co111$mean.of.class, co111$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset111")
#' 
#' #---------------------------------------------#
#' # smooth landscape, Gaussian distribution, omitted predictor
#' 
#' simData("112")
#' 
#' d112 <- extract.ncdf("dataset112.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d112) # levelplot response
#' 
#' fm112 <- glm(y ~ x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d112, family = "gaussian") # omit x1   
#' summary(fm112) 
#' res112 <- residuals(fm112) # calculate residuals
#' co112 <- correlog(d112$Lat, d112$Lon, res112, increment=0.02, resamp=1) # check autocorrleation
#' plot(co112$mean.of.class, co112$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset112")
#' 
#' #---------------------------------------------#
#' # smooth landscape, Bernoulli distribution, refrence data
#' 
#' simData("120")
#' 
#' d120 <- extract.ncdf("dataset120.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d120) # levelplot response
#' 
#' fm120 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d120, family = "gaussian")   
#' summary(fm120) 
#' res120 <- residuals(fm120) # calculate residuals
#' co120 <- correlog(d120$Lat, d120$Lon, res120, increment=0.02, resamp=1) # check autocorrleation
#' plot(co120$mean.of.class, co120$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset120")
#' 
#' #---------------------------------------------#
#' # smooth landscape, Bernoulli distribution, SAC onto response
#' 
#' simData("121")
#' 
#' d121 <- extract.ncdf("dataset121.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d121) # levelplot response
#' 
#' fm121 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d121, family = "gaussian")   
#' summary(fm121) 
#' res121 <- residuals(fm121) # calculate residuals
#' co121 <- correlog(d121$Lat, d121$Lon, res121, increment=0.02, resamp=1) # check autocorrleation
#' plot(co121$mean.of.class, co121$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset121")
#' 
#' #---------------------------------------------#
# smooth landscape, Bernoulli distribution, omitted predictor
#' 
#' simData("122")
#' 
#' d122 <- extract.ncdf("dataset122.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d122) # levelplot response
#' 
#' fm122 <- glm(y ~ x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d122, family = "gaussian") # omit x1   
#' summary(fm122) 
#' res122 <- residuals(fm122) # calculate residuals
#' co122 <- correlog(d122$Lat, d122$Lon, res122, increment=0.02, resamp=1) # check autocorrleation
#' plot(co122$mean.of.class, co122$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset122")
#' 
#' #---------------------------------------------#
#' # real landscape, Gaussian distribution, refrence data
#' 
#' simData("310")
#' 
#' d310 <- extract.ncdf("dataset310.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d310) # levelplot response
#' 
#' fm310 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d310, family = "gaussian")   
#' summary(fm310) 
#' res310 <- residuals(fm310) # calculate residuals
#' co310 <- correlog(d310$Lat, d310$Lon, res310, increment=0.16, resamp=1) # check autocorrleation
#' plot(co310$mean.of.class, co310$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset310")
#' 
#' #---------------------------------------------#
#' # real landscape, Gaussian distribution, SAC onto response
#' 
#' simData("311")
#' 
#' d311 <- extract.ncdf("dataset311.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d311) # levelplot response
#' 
#' fm311 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d311, family = "gaussian")   
#' summary(fm311) 
#' res311 <- residuals(fm311) # calculate residuals
#' co311 <- correlog(d311$Lat, d311$Lon, res311, increment=0.16, resamp=1) # check autocorrleation
#' plot(co311$mean.of.class, co311$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset311")
#' 
#' #---------------------------------------------#
#' # real landscape, Gaussian distribution, omitted predictor
#' 
#' simData("312")
#' 
#' d312 <- extract.ncdf("dataset312.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d312) # levelplot response
#' 
#' fm312 <- glm(y ~ x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d312, family = "gaussian") # omit x1   
#' summary(fm312) 
#' res312 <- residuals(fm312) # calculate residuals
#' co312 <- correlog(d312$Lat, d312$Lon, res312, increment=0.16, resamp=1) # check autocorrleation
#' plot(co312$mean.of.class, co312$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset312")
#' 
#' #---------------------------------------------#
#' # real landscape, Bernoulli distribution, refrence data
#' 
#' simData("320")
#' 
#' d320 <- extract.ncdf("dataset320.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d320) # levelplot response
#' 
#' fm320 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d320, family = "gaussian")   
#' summary(fm320) 
#' res320 <- residuals(fm320) # calculate residuals
#' co320 <- correlog(d320$Lat, d320$Lon, res320, increment=0.16, resamp=1) # check autocorrleation
#' plot(co320$mean.of.class, co320$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset320")
#' 
#' #---------------------------------------------#
#' # real landscape, Bernoulli distribution, SAC onto response
#' 
#' simData("321")
#' 
#' d321 <- extract.ncdf("dataset321.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d321) # levelplot response
#' 
#' fm321 <- glm(y ~ x1 + x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d321, family = "gaussian")   
#' summary(fm321) 
#' res321 <- residuals(fm321) # calculate residuals
#' co321 <- correlog(d321$Lat, d321$Lon, res321, increment=0.16, resamp=1) # check autocorrleation
#' plot(co321$mean.of.class, co321$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset321")
#' 
#' #---------------------------------------------#
#' # real landscape, Bernoulli distribution, omitted predictor
#' 
#' simData("322")
#' 
#' d322 <- extract.ncdf("dataset322.nc")[[2]] # extract data
#' levelplot(y~Lon+Lat,data=d322) # levelplot response
#' 
#' fm322 <- glm(y ~ x4 + I(x4^2) + x3*x4 + x3 + x2 + x5 + x6 + x7, 
#'              data = d322, family = "gaussian") # omit x1   
#' summary(fm322) 
#' res322 <- residuals(fm322) # calculate residuals
#' co322 <- correlog(d322$Lat, d322$Lon, res322, increment=0.16, resamp=1) # check autocorrleation
#' plot(co322$mean.of.class, co322$correlation, type = "o", ylim = c(-1,1), # plot correlogram
#'      ylab="Moran Similarity", xlab="averaged distance class", main = "dataset322")
#' }
simData <- function(dataset, # three numbers character string: landscape, dsitribution, SA scenario
                    # landscape: 1 = smooth, 2 = realistic, 3 = real
                    # distribution: 1 = normal, 2 = binomial, 3 = poisson
                    # scenario: 0 = reference (no SAC), 1 = SA onto response, 2 = omitted predictor, 
                    #3 = wrong functional form, 4 = dispersal
                    filename = "default",
                    gridsize = c(50L, 50L), # vector of length 2: number of cells in 
                    #x and y direction, respectively
                    cvfold = 10L, # 
                    cvblock.size = c(10,10), # number of cells in x, y direction in one block
                    r.seed = 20151126, # control randomisation seed
                    n.predictor = 7L, # Number of predictors
                    f.smooth = list(function() lon, # lon and lat as expanded grid of X and Y values
                                    function() lat, 
                                    function() (lon - mean(lon))^2, 
                                    function() (lat - mean(lat))^2,
                                    function() x3^x4 * x4^x3,
                                    function() x1^x1 * x3^x4,
                                    function() x2^x1 * x4^x3 * log(x5 + 1)),
                    f.realistic = list(var = 0.1, 
                                       scale = 0.1), 
                    f.real = list(resolution = 10L, # 2.5, 5, 10 minutes of a degree
                                  extent = c("5N24E", "7S37E"), # provide in geogr. coordinates: Diagonal corners.
                                  bio.vars = c("bio1", "bio19", "bio2", "bio12", "bio4", "bio18", "bio3")), 
                    f.response = c("x1", "x4", "x4^2", "x3*x4", "x3"),
                    par.response = "default", # dependant on distribution (dataset[2])
                    f.sac1 = list(corCoef = -3,
                                  sarFactor = 1),
                    f.sac2 = "x1", # predictor to omit in model structure
                    f.sac3 = c("^","*"), # functional forms to omitin model structure, here quadratic effects and interactions
                    f.sac4 = list(dispersal.max = 0.1, # maximum dispersal rate
                                  dispersal.shape = 30), # the higher the more skewed
                    interactive = TRUE){ # when FALSE, existing files will be replaced, download will be processed

  # check if dataset correctly specified
  if (!is.character(dataset)) stop("dataset must be a character string of three numbers defining\n landscape [1;3], distribution [1;3] and SAC scenario [0;4].")
  # extract distribution, landscape, SA scenario
  dataset <- unlist(strsplit(dataset, ""))
  
  # check if dataset correctly specified 
  if (length(dataset) != 3) stop("define dataset as three cypher character string:\n landscape, distribution, SAC scenario, 
                                e.g. 110 for smooth landscape, gaussian distribution and reference scenario (no SAC)")
  if (!dataset[1] %in% 1:3) stop("There are three landscapes available (1-3):\n 1 = smooth, 2 = realistic and  3 = real, please pick one of them")
  if (!dataset[2] %in% 1:3) stop("There are three distributions available (1-3):\n 1 = gaussian, 2 = bernoulli and 3 = poisson, please pick one of them")
  if (!dataset[3] %in% 0:4) stop("There are five SAC scenarios available (0-4):\n 0 = reference (no SAC), 1 = SAC onto response, 2 = omitted predictor,\n 3 = wrong functional form and 4 = dispersal, please pick one of them")
  
  # set filename
  # if default
  if (filename == "default") filename <- paste0("dataset",Reduce(paste0,dataset),".nc")
  # if not check whether extension (.nc) needs to appended
  if (length(grep(".nc", filename)) == 0){
    filename <- paste0(filename,".nc")
    } else {
      check.nc <- strsplit(filename, split = "")[[1]]
      if(!Reduce(paste0, check.nc[(length(check.nc)-2):length(check.nc)]) == ".nc"){
        filename <- paste0(filename,".nc")
      }
    }
  
  # check if file already exists. If so, ask if file should be overwritten
  if (file.exists(filename)){
    if (interactive){
      check.overwrite <- keep.asking(Q = "File already exists. Do you want to overwrite it?")
      if(check.overwrite == "n") stop("File already exists, please specify new filename.")
    }
  }
  
  #--------------------------------------------------------------------------#
  # set up list for readme 
  readme <- vector("list", 5L)
  #--------------------------------------------------------------------------#
  #--------------------------------------------------------------------------#
  
  # simulate predictors
  
  # define grid
  # lon lat coords
  Xvec <- seq(0, 1, len = gridsize[1]); Yvec <- seq(0, 1, len = gridsize[2]) 
  grd <- expand.grid(X = Xvec, Y = Yvec)
  lon <- grd$X; lat <- grd$Y
  
  #%%%%%%%%%%%%%%%%%%%%%#
  # 1) smooth landscape #
  if(dataset[1] == 1){
    # there must be as many functions specified in f.smooth as n.predictor
    if(n.predictor != length(f.smooth)) stop("n.predictor must be equivalent to the number of functions\n provided in f.smooth to simulate smooth predictors.")
    # call functions from f.smooth
    for(i in seq(n.predictor)) assign(paste0("x",i), do.call(f.smooth[[i]], args = list()))
    # write functions in readme entry 1
    functional.form <- sapply(seq.int(n.predictor), 
                              function(x) paste0("x",x," = ",deparse(f.smooth[[x]]), ";")[2])
    readme[[1]] <- paste0("Smooth landscape: The ", n.predictor, " predictors are linear and ",
                          "non-linear gradients without noise. ", 
                          "Latitude (lat) and longitude (lon) are the x and y column, ", 
                          "respectively, from an expanded grid of two sequences from 0 to 1.", 
                          "The functional form of the gradients is defined in f.smooth: ", 
                          Reduce(paste, functional.form), "These ", n.predictor, 
                          " predictors are rescaled to [-1,1].")
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%#
  # 2) realistic landscape #
  if(dataset[1] == 2){
    # recycle first value of f.realistic$var if not of length = n.predictor
    if(length(f.realistic$var) != n.predictor){
      f.var <- rep(f.realistic$var[1], n.predictor)
    }else{
      f.var <- f.realistic$var
    }
    # recycle first value of f.realistic$scale if not of length = n.predictor
    if(length(f.realistic$scale) != n.predictor){
      f.scale <- rep(f.realistic$scale[1], n.predictor)
    }else{
      f.scale <- f.realistic$scale
    }
    
    set.seed(r.seed)
    for(i in seq.int(n.predictor)){
      # compute the exponential covariance model 
      expCov <- RMexp(var = f.var[i], scale = f.scale[i])
      # simulate expCov on defined grid
      assign(paste0("x",i), as.vector(RFsimulate(expCov, 
                                                 x = Xvec, y = Yvec, 
                                                 spConform = FALSE)))
    }
    readme[[1]] <- paste0("Realistic landscape: The ", n.predictor, " predictors are ",
                          "simulated unconditional Gaussian random fields from exponential ",
                          "covariance models with variance = ", paste(f.var, collapse=","), 
                          " and scale = ", paste(f.scale, collapse=","),". These ", 
                          n.predictor, " predictors are rescaled to [-1,1].")
  }
  
  #%%%%%%%%%%%%%%%%%%%#
  # 3) real landscape #
  if(dataset[1] == 3){
    # check if data already exists
    if(!file.exists("wc10")){
      if(interactive){
        check.download <- keep.asking(Q = "Do you want to download data from http://www.worldclim.org?")
        if(check.download == "y"){ # check if data should be downloaded
          bio <- getData('worldclim', download = TRUE, var='bio', res=f.real$resolution)
        }else{ # if not, but file does not exist: error
          stop(paste0("If data shall not be downloaded, please provide them in the working directory.
                      Folder must be named: wc",f.real$resolution))
        }
      }else{ # if not interactive and file not existing, do download
        bio <- getData('worldclim', download = TRUE, var='bio', res=f.real$resolution)
      }
    }else{ # if data alread exists, load from default location
      bio <- getData('worldclim', download = FALSE, var='bio', res=f.real$resolution)
    }
    
    # define area of interest as extent obj
    ext <- sapply(1:2, function(x) geo.to.num(f.real$extent[x]))
    # sort, to go from lower left to upper right 
    ext <- extent(c(sort(ext[2,]), sort(ext[1,])))              
    # crop data to extent
    bio <- crop(bio, ext)
    # get gridsize for real landscape
    gridsize <- dim(bio)[1:2]
    
    # define grid for real landscape
    # lon lat coords
    grd <- as.data.frame(coordinates(bio))
    lon <- grd$x; lat <- grd$y
    
    # write real predictors as indiv. vectors x_1 : x_n.predictor
    for(i in seq.int(n.predictor)) assign(paste0("x",i), values(bio[[f.real$bio.vars[i]]]))
    readme[[1]] <- paste0("Real landscape: Real bio-climatic predictors ",
                          "from http://www.worldclim.org for the extent of ",
                          Reduce(paste, paste(c("xmin =", ", xmax =", ", ymin =", ", ymax ="),f.real$extent)),
                          " decimal degree, and a resolution of ", f.real$resolution, " minutes of a degree. ",
                          Reduce(paste, paste0(f.real$bio.vars, ",")), " have been pre-selected ",
                          "based on spearman's rho^2. These ", n.predictor, " predictors are rescaled to [-1,1].")
  }
  
  
  # rescale predictors to [-1,1] # same scale for all landscapes
  allpredictors <- numeric(n.predictor)
  for(i in seq.int(n.predictor)){
    xi <- get(paste0("x",i))
    scale.factor <- (1-(-1))/(max(xi, na.rm = TRUE) - min(xi, na.rm = TRUE))
    assign(paste0("x",i), -1 + (xi - min(xi)) * scale.factor)
    # predictor names to extract nuisances from
    allpredictors[i] <- paste0("x",i)
  }
  
  
  
  
  #--------------------------------------------------------------------------#
  # response
  
  # compute dist matrix for scenario 1 or 4
  if (dataset[3] %in% c(1,4)){
    D <- dist.matrix(as.matrix(grd), method = "euclidean")
    # Dispersal
    if (dataset[3] == 4){
      D.scaled <- (D - min(D, na.rm=T))/(max(D, na.rm=T)-min(D, na.rm=T)) # scale to [0;1]
      # dispersal factor from variably skewed exponential curve
      disp <- f.sac4$dispersal.max / exp(D.scaled * f.sac4$dispersal.shape) 
      diag(disp) <- 0
    }
  }
  
  # 1) normal distribution #
  if (dataset[2] == 1){
    readme[[2]] <- "Gaussian"
    set.seed(r.seed)
    # set default parameters
    if (par.response[1] == "default"){
      par.gaus <- c(0.8,0.2, 0.9, -0.8,-0.6,0.5,0.2)
    } else {
      if (length(par.response)-1 != length(f.response)+1) 
         stop("In case of a Gaussian distribution you must provide
               par.response as numeric vector of length
               length(f.response)+1")
      par.gaus <- par.response
    }
    # compute response variable using specified terms in f.response and parameters
    resp.formula <- Reduce(paste, c(paste(par.gaus[1],"+"),
                                    paste(paste(par.gaus[2:(length(par.gaus)-1)],"*",f.response),
                                          collapse="+", sep = " ")))
    y <- eval(parse(text = resp.formula)) + rnorm(prod(gridsize), 0, par.gaus[length(par.gaus)]) # add Gaussian noise
    if (dataset[3] == 4){ # if dispersal
      y <- sapply(seq_along(y), function(x) y[x] * (1+ sum(disp[ ,x] * y)))
    }
  }
  
  # 1) binomial distribution #
  if (dataset[2] == 2){
    readme[[2]] <- "Bernoulli"
    # set default parameters
    if (par.response[1] == "default"){
      par.bern <- c(0.2, 4.5, 1.2,-1.2,-1.1,0.9)
    } else {
      if (length(par.response)-1 != length(f.response)) 
         stop("In case of a Bernoulli distribution you must provide
               par.response as numeric vector of same length
               as f.response")
      par.bern <- par.response
    }
    # compute logit 
    resp.formula <- Reduce(paste, c(paste(par.bern[1],"+"),
                                    paste(paste(par.bern[2:length(par.bern)],"*",f.response),
                                          collapse="+", sep = " ")))
    pi.log <- eval(parse(text = resp.formula))
    if (dataset[3] %in% c(0,2,3)){
      set.seed(r.seed)
      y <- rbinom(prod(gridsize), size = 1, prob=plogis(pi.log))
    } else {
      if (dataset[3] == 4){ # if dispersal
        disp.prob <- sapply(seq_along(pi.log), function(x) plogis(pi.log[x]) * 
                              prod((1 + disp[ ,x] * plogis(pi.log))))
        disp.prob <- ifelse(disp.prob > 1, 1, disp.prob)
        set.seed(r.seed)
        y <- rbinom(prod(gridsize), size = 1, prob=disp.prob)
      }
    }
  }

  # 1) zero-inflated poisson distribution #
  if (dataset[2] == 3){
    readme[[2]] <- "zero-inflated Poisson"
    # set default parameters
    if (par.response[1] == "default"){
      par.bern <- c(0.2, 4.5, 1.2, -1.2, -1.1, 0.9)
      par.pois <- c(0.2, 1.6, 0.9, -0.8, -0.8, 0.5)
    } else {
      if (length(par.response) != 2 | length(par.response[[1]])-1 != length(f.response) | 
         length(par.response[[2]])-1 != length(f.response)) stop("In case of a Poisson 
               distribution you must provide par.response as a list of two numeric vectors, both 
               of the same length as f.response. The first for the Bernoulli distribution, 
               the second for the then zero-inflated Possion")
      par.bern <- par.response[[1]]
      par.pois <- par.response[[2]]
    }
    set.seed(r.seed)
    # zero-inflation from Bernoulli distribution
    zeroone <- rbinom(prod(gridsize), size = 1, 
                      prob=plogis(eval(parse(text = Reduce(paste, 
                            c(paste(par.bern[1],"+"),
                  paste(paste(par.bern[2:length(par.bern)],"*",
     f.response), collapse="+", sep = " ")))))))
    resp.formula <- Reduce(paste, c(paste(par.pois[1],"+"),
                                    paste(paste(par.pois[2:length(par.pois)],"*",f.response),
                                          collapse="+", sep = " ")))
    pi <- exp(eval(parse(text = resp.formula)) + rnorm(prod(gridsize), sd=0.2)) # add a bit of noise
    set.seed(r.seed)
    # use poisson error for cases where y = 1
    y <- rpois(prod(gridsize), lambda = pi)
    y <- ifelse(zeroone == 0, 0, y)
    if(dataset[3] == 4){ # if dispersal
      y <- sapply(seq_along(y), function(x) as.integer(y[x] * (1+ sum(disp[ ,x] * y))))
    }
  }
  
  # add attribute for coefficients
  readme[[5]] <- resp.formula
  
  #--------------------------------------------------------------------------#
  # nuisance
  to.omit <- unlist(sapply(f.response, function(x) grep(x, allpredictors, fixed = TRUE)))
  nuisance <- allpredictors[-to.omit]
  # SAC causes
  
  # 0) reference scenario (no SAC) #
    if(dataset[3] == 0){
      # readme
      readme[[3]] <- "Reference, i.e. no SAC"
      readme[[4]] <- paste("Model structure to use: y ~",
                            paste(f.response, collapse = " + "), "+",
                            paste(nuisance, collapse = " + "))
    }
  
  # 1) SAC onto response #
  if (dataset[3] == 1){
    readme[[3]] <- paste0("SAC onto response variable. ",
                          "Correlation structure: e^(",f.sac1$corCoef," * distance.matrix). ",
                          "SAC error added to y was multiplied by ",f.sac1$sarFactor)
    readme[[4]] <- paste("Model structure to use: y ~",
                         paste(f.response, collapse = " + "), "+",
                         paste(nuisance, collapse = " + "))
    
    OMEGA <- exp(f.sac1$corCoef * D) # correlation structure
    W <- chol(chol2inv(chol(OMEGA))) # correlation weights
    WInv <- solve(W) # W1inv = inverse of W
    
    set.seed(r.seed); err <- WInv %*% rnorm(dim(D)[1]) # produces correlated normal errors
    err <- err - mean(err) # center error
    
    if (dataset[2] == 1) y <- y + err * f.sac1$sarFactor
    if (dataset[2] == 2){
      pi <- plogis(pi.log)
      set.seed(r.seed)
      # y <- rbinom(prod(gridsize), size = 1, prob=err * sqrt(pi * (1 - pi)) + pi)
      y <- ifelse((err * sqrt(pi * (1 - pi)) + pi) < 0.5, 0, 1)
    }
    if (dataset[2] == 3){
      y <- as.integer(abs(err * sqrt(pi) + pi))
    }
  }

  # 2) omitted predictor #
    if(dataset[3] == 2){
      
      # readme
      readme[[3]] <- "Omitted predictor."
      to.omit <- unlist(sapply(f.sac2, function(x) grep(x, f.response, fixed = TRUE)))
      readme[[4]] <- paste("Model structure to use: y ~",
                           paste(f.response[-to.omit], collapse = " + "), "+",
                           paste(nuisance, collapse = " + "))
    }
  
  # 3) wrong functional form #
    if(dataset[3] == 3){
      
      # readme
      readme[[3]] <- "Wrong functional form: Purposely miss quadratic effect."
      to.omit <- unlist(sapply(f.sac3, function(x) grep(x, f.response, fixed = TRUE)))
      readme[[4]] <- paste("Model structure to use: y ~",
                           paste(f.response[-to.omit], collapse = " + "), "+",
                           paste(nuisance, collapse = " + "))
    }
  
  # 4) dispersal #
    if(dataset[3] == 4){
      
      # readme
      readme[[3]] <- "Dispersal."
      readme[[4]] <- paste("Model structure to use: y ~",
                           paste(f.response, collapse = " + "), "+",
                           paste(nuisance, collapse = " + "))
    }
  
  #--------------------------------------------------------------------------#
  # CV block ID
  cvblock <- matrix(NA, nrow = gridsize[2], ncol = gridsize[1])
  
  xblocks <- gridsize[1] / cvblock.size[1] 
  yblocks <- gridsize[2] / cvblock.size[2]
  
  # first fill full blocks
  block.counter <- 0
  for(xs in seq(floor(xblocks))){
    for(ys in seq(floor(yblocks))){
      block.counter <- block.counter + 1
      cvblock[(((ys-1) * cvblock.size[2] + 1) : (ys * cvblock.size[2])),
              (((xs-1) * cvblock.size[1] + 1) : (xs * cvblock.size[1]))] <- block.counter
    }
  }
  # check half blocks in x direction
  if(xblocks - floor(xblocks) != 0){
    for(ys2 in seq(floor(yblocks))){
      block.counter <- block.counter + 1
      cvblock[(((ys2-1) * cvblock.size[2] + 1) : (ys2 * cvblock.size[2])),
              ((floor(xblocks) * cvblock.size[1]+1) : (xblocks * cvblock.size[1]))] <- block.counter
    }
  }
  # check half blocks in y direction
  if(yblocks - floor(yblocks) != 0){
    for(xs2 in seq(floor(xblocks))){
      block.counter <- block.counter + 1
      cvblock[((floor(yblocks) * cvblock.size[2]+1) : (yblocks * cvblock.size[2])),
              (((xs2-1) * cvblock.size[1] + 1) : (xs2 * cvblock.size[1]))] <- block.counter
    }
  }
  # check half blocks in x and y direction
  if(xblocks - floor(xblocks) != 0 & yblocks - floor(yblocks) != 0){
    lock.counter <- block.counter + 1
    cvblock[is.na(cvblock)] <- block.counter
  } 
    
  # smooth and realistic landscapes require matrix transponation 
  if(dataset[1] == 3) cvblock <- c(cvblock) else
    cvblock <- c(t(cvblock))
  
  # re-sample to cvfold
  if (cvfold > length(unique(as.vector(cvblock)))) warning("You have specified more cross-validation folds (cvfold) than blocks available.\n To change this, increase gridsize, decrease cvblock.size or decrease cvfold.")
  
  set.seed(r.seed)
  # stratified sampling, randomly select cvfold elements from cvblock and assign 1:cvfold
  cvfold.sample <- rep(NA, length(unique(as.vector(cvblock)))) # as.vector needed to get the number of blocks right: CFD
  while (sum(is.na(cvfold.sample)) > 0){
    whichNA <- which(is.na(cvfold.sample))
    index <- sample(whichNA, size = ifelse(length(whichNA) > cvfold, cvfold, length(whichNA)))
    cvfold.sample[index] <- sample(seq.int(cvfold), size = ifelse(length(index) == cvfold, cvfold, length(index)))
  }
#   cvfold.sample <- sample(cvfold, size = length(unique(cvblock)), replace = 
#                             ifelse(cvfold < length(unique(cvblock)), TRUE, FALSE))
  cvblock <- cvfold.sample[cvblock]
#   for(i in unique(cvblock)) cvblock[cvblock == i] <- cvfold.sample[i]
  
  #--------------------------------------------------------------------------#
  #--------------------------------------------------------------------------#
  
  # write netcdf
  dim1 <- ncdim_def("row", "observation", 1:prod(gridsize))
  dim2 <- ncdim_def("column", "variable", 1:(4+n.predictor))

  var <- vector(mode = "list", length = 4L + n.predictor)
  var[[1]] = ncvar_def("CVid","", dim1, NA, longname = "Block ID for cross-validation")
  var[[2]] = ncvar_def("Lat","", dim1, NA, longname = "Latitude")
  var[[3]] = ncvar_def("Lon","", dim1, NA, longname = "Longitude")
  var[[4]] = ncvar_def("y","", dim1, NA, longname = "Response Variable")
  for(i in seq.int(n.predictor)){ 
    var[[i+4]] <- ncvar_def(paste0("x",i),"", dim1, NA, 
                               longname = paste0("Predictor ", i))
  }

  nc <- nc_create(paste0(filename), var)
  
  

  ncvar_put(nc, var[[1]], cvblock)
  ncvar_put(nc, var[[2]], lat)
  ncvar_put(nc, var[[3]], lon)
  ncvar_put(nc, var[[4]], y)
  for(i in seq.int(n.predictor)){ 
    ncvar_put(nc, var[[4+i]], get(paste0("x", i)))
  }
  # write global attributes into ncdf
  # General information
  ncatt_put(nc, 0, attname = "genInfo", 
               attval = paste0("Gridsize = ",gridsize[1], 
                               "x", gridsize[2], " cells;",
                               " Random seed = ", r.seed, ";",
                               " Number of predictors = ", n.predictor, ";",
                               " Block size for CV = ", cvblock.size[1], "x", 
                               cvblock.size[2], " cells"), 
               prec = "text")
  
  # General information: gridsize, random seed, no. of predictors, cv block size
  ncatt_put(nc, 0, attname = "landsc", 
               attval = readme[[1]], 
               prec = "text")
  # distribution
  ncatt_put(nc, 0, attname = "distr", 
               attval = readme[[2]], 
               prec = "text")
  # SAC scenario
  ncatt_put(nc, 0, attname = "sacScen", 
               attval = readme[[3]], 
               prec = "text")
  # lmodel structure
  ncatt_put(nc, 0, attname = "model", 
               attval = readme[[4]], 
               prec = "text")
  # response coefficients
  ncatt_put(nc, 0, attname = "coeffs",
            attval = readme[[5]], 
            prec = "text")
  
  nc_close(nc)
  return(invisible(""))
  #--------------------------------------------------------------------------#
}







