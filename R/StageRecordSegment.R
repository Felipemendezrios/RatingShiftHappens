#' Recession extraction
#'
#' All recession from the stage record are extracted according to several criteria. NA values are not accepted in the stage record.
#'
#' @param H real vector , stage
#' @param uH real vector, uncertainty of the stage
#' @param time vector, time in POSIXct format is mandatory
#' @param filter.H integer value, filter for number of step records. 10 means only one stage sample every 10 is kept
#' @param chi real value, maximum stage rise between two recession data
#' @param delta.t.min integer value, minimum days between two recession data
#' @param delta.t.max integer value, maximum days between two recession data
#' @param tgood integer value, minimum length of the recession in days
#' @param Nmin.rec integer value, minimum number of data in a recession curve
#' @param tburn.rec numeric, burn factor, >=0 and <1. 0.4 means the first 40 percent of the recession samples are discarded.
#'
#' @return list, all recession satisfying the specified criteria
#' @export
#'
#' @examples
#'
#' recessions=Extraction_recession(H=ArdecheRiverMeyrasStage$H,
#'                                 uH=0.5,
#'                                 time=ArdecheRiverMeyrasStage$Date,
#'                                 chi=1.5,
#'                                 tgood=30)
#' plot_rec_extracted(Rec_extracted = recessions)
Extraction_recession <- function(H,
                                 uH,
                                 time,
                                 filter.H=1,
                                 chi=stats::quantile(H,probs = 0.95,na.rm = TRUE),
                                 delta.t.min=0,
                                 delta.t.max=10,
                                 tgood=20,
                                 Nmin.rec=10,
                                 tburn.rec=0.2){

  if(any(is.na(H) | is.na(uH) | is.na(time)))stop('NA value in input data. Be sure to remove them before running function')
  if(is.null(check_vector_lengths(H,time)))stop('The input data have not the same length')
  if(filter.H  < 1 )stop('filter.H must be strictly positive and gretear than or equal to 1')
  if(filter.H != round(filter.H))stop('filter.H must be a integer')
  if(delta.t.min != round(delta.t.min))stop('delta.t.min must be a integer')
  if(delta.t.min < 0)stop('delta.t.min must be strictly positive')
  if(delta.t.max != round(delta.t.max))stop('delta.t.max must be a integer')
  if(delta.t.max < 0)stop('delta.t.max must be strictly positive')
  if(tgood != round(tgood))stop('tgood must be a integer')
  if(tgood < 0)stop('tgood must be strictly positive')
  if(Nmin.rec < 0)stop('Nmin.rec must be strictly positive')
  if(!lubridate::is.POSIXct(time))stop('time must be POSIXct format')
  if(!dplyr::between(tburn.rec,0,1))stop('tburn.rec must be in [0;1)')

  stage.record = data.frame(h=H, uH=uH,t=time)
  # Order by time
  stage.record=stage.record[order(stage.record$t),]

  # check for NA:
  if (any(is.na(stage.record$h))){
    stage.record = stage.record[-which(is.na(stage.record$h)),]
  }

  # filter stage record by user specifications (filter.H):
  stage.record.post = stage.record[seq(1, nrow(stage.record), filter.H), ]

  # # Convert time format
  # stage.record.post$t=as.numeric(stage.record.post$t)

  # find all local decreasing points:
  min_loc.h = localmin( t = stage.record.post$t,
                        h = stage.record.post$h,
                        uH = stage.record.post$uH )

  data_min = data.frame(tmin=min_loc.h$t_local_min,
                        hmin=min_loc.h$h_local_min,
                        uHmin=min_loc.h$uH_local_min)

  # Find all recession curves :

  # Filter N°1 : delta.t.min
  # Remove data which do not respect the specify time step
  diff_tmin=diff(data_min$tmin)
  units(diff_tmin) <- 'days'

  if(any(diff_tmin<delta.t.min,na.rm = TRUE)){
    data_min_filter_1=data_min[-which(diff_tmin<delta.t.min),]
  }else{
    data_min_filter_1=data_min
  }

  # Filter N°2 : hrec
  # Finding the all the minimum of the minimum (when hmin(i)< hmin(i-1))
  # Filter N°3 : chi and delta.t.max
  # chi is a threshold that defines the end of one recession and the beginning of another
  # delta.t.max helps to handle long time distances between  points

  # initialize:
  hrec = uHrec = data_rec = c()
  trec =  as.POSIXct(character(0))

  hrec[1] = data_min_filter_1$hmin[1]
  trec[1] = data_min_filter_1$tmin[1]
  uHrec[1] = data_min_filter_1$uHmin[1]
  j = 2
  m = 1
  k = 1

  j_max = length(data_min_filter_1$tmin)

  while (j <= j_max) {

    diff_tmin_delta_max=difftime(data_min_filter_1$tmin[j],trec[m],units='days')
    if ((data_min_filter_1$hmin[j] <= hrec[m]) & (diff_tmin_delta_max < delta.t.max)) {
      # ok, recession point
      m = m + 1
      hrec[m] = data_min_filter_1$hmin[j]
      trec[m] = data_min_filter_1$tmin[j]
      uHrec[m] = data_min_filter_1$uHmin[j]
    } else {

      if (((abs(data_min_filter_1$hmin[j] - hrec[m])) >= chi ) | diff_tmin_delta_max >= delta.t.max) { # stop the recession
        # Filter N°4 : tgood
        # Select only recessions with length > tgood
        # Filter N°5 : Nmin
        # Select only recessions with a number of point > Nmin :

        if(((difftime(rev(trec)[1],trec[1],units='days')) >= tgood) & (length(trec) >= Nmin.rec)){
          data_rec_unburned=data.frame(date=trec,
                                       hrec=hrec,
                                       uHrec=uHrec)
          # Apply burn of the first part of the recession to reduce data information
          remove_lines=round(nrow(data_rec_unburned)*tburn.rec)
          data_rec_burned=data_rec_unburned[-c(1:remove_lines),]
          data_rec_burned_plot=cbind(time_rec=as.numeric(
                                       difftime(data_rec_burned$date,
                                                data_rec_burned$date[1], units = "days")),
                                     data_rec_burned,
                                     indx=k)


          data_rec = rbind(data_rec,data_rec_burned_plot)
          k = k + 1
        }

        #Reinitialize variables
        hrec = uHrec = c()
        trec =  as.POSIXct(character(0))
        hrec[1] = data_min_filter_1$hmin[j]
        trec[1] = data_min_filter_1$tmin[j]
        uHrec[1] = data_min_filter_1$uHmin[j]
        m = 1
      }
    }
    j = j+1
  }
  return(Rec_extracted = data_rec)
}

#' Model and segmentation of recession
#'
#' Modelling recession using a catalog of fit models available to apply segmentation procedure.
#' To get the catalog of models available using `GetCatalog()`, all functions specifying Recession in the fit are concerned.
#' For more information about recession model and their default values, please go to the fit function, e.g. `?fitRecession_M3`
#'
#' @param time_rec real vector, recession duration relative to the first data detected during the recession
#' @param hrec real vector, stage value of the recessions
#' @param uhrec real vector, uncertainty of stage value of the recessions
#' @param indx  integer, factor used to gather the data of a same recession
#' @param nSmax integer, maximum number of segments to assess
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#' @param funk the function for estimating the recession
#' @param ... optional arguments to funk
#'
#' @return to complete!
#' @export
#'
#' @examples
#' recessions=Extraction_recession(H=ArdecheRiverMeyrasStage$H,
#'                                 uH=0.5,
#'                                 time=ArdecheRiverMeyrasStage$Date,
#'                                 chi=1.5,
#'                                 tgood=30)
#'
#' D=RBaM::dataset(X=recessions['time_rec'],
#'                 Y=recessions['hrec'],
#'                 Yu=recessions['uHrec'],
#'                 VAR.indx=recessions['indx'],
#'                 data.dir=file.path(tempdir(),'BaM','Recession'))
#'
#' # Give prior knowledge about recession-specific parameters and keep prior on stable parameters
#' alpha1.object = RBaM::parameter_VAR(name='alpha1',
#'                                     index='indx',
#'                                     d=D,
#'                                     init=rep(100,max(recessions$indx)),
#'                                     prior.dist=rep('Uniform',max(recessions$indx)),
#'                                     prior.par=rep(list(c(0,1000)),max(recessions$indx)))
#'
#' alpha2.object = RBaM::parameter_VAR(name='alpha2',
#'                                     index='indx',
#'                                     d=D,
#'                                     init=rep(50,max(recessions$indx)),
#'                                     prior.dist=rep('Uniform',max(recessions$indx)),
#'                                     prior.par=rep(list(c(0,100)),max(recessions$indx)))
#'
#' beta.object = RBaM::parameter_VAR(name='beta',
#'                                     index='indx',
#'                                     d=D,
#'                                     init=rep(52.4,max(recessions$indx)),
#'                                     prior.dist=rep('Uniform',max(recessions$indx)),
#'                                     prior.par=rep(list(c(-97.6,102.4)),max(recessions$indx)))
#'
#' model_rec=ModelAndSegmentation.recession.regression(time_rec=recessions$time_rec,
#'                                                     hrec=recessions$hrec,
#'                                                     uhrec=recessions$uHrec,
#'                                                     indx=recessions$indx,
#'                                                     alpha1.object=alpha1.object,
#'                                                     alpha2.object=alpha2.object,
#'                                                     beta.object=beta.object,
#'                                                     Dataset.object=D,
#'                                                     nCyclesrec=100,
#'                                                     burnrec=0.5,
#'                                                     nSlim=10)
#'
ModelAndSegmentation.recession.regression <- function(time_rec,
                                                      hrec,
                                                      uhrec,
                                                      indx,
                                                      nSmax=2,
                                                      nMin= 1,
                                                      nCycles=100,
                                                      burn=0.5,
                                                      nSlim=max(nCycles/10,1),
                                                      temp.folder=file.path(tempdir(),'BaM'),
                                                      funk=fitRecession_M3,...){
  # Check information given in input
  if(time_rec<0)stop('time_rec must be positive')
  if(uhrec<0)stop('uhrec must be positive')
  if(indx<0)stop('indx must be positive')
  if(indx[1]!=1)stop('Fist number of indx must be one to start the sequence')
  if(any(diff(unique(indx))!=1))stop('indx must be a sequence of consecutive numbers')

return(fitRecession_M3(...))

}
