#' Recession extraction
#'
#' All recession from the stage record are extracted according to several criteria. NA values are not accepted in the stage record.
#'
#' @param H real vector , stage (m)
#' @param uH real vector, uncertainty of the stage
#' @param time vector, time in POSIXct format is mandatory
#' @param filter.H integer value, filter for number of step records. 10 means only one stage sample every 10 is kept
#' @param chi real value, maximum stage rise between two recession data (m)
#' @param delta.t.min integer value, minimum days between two recession data
#' @param delta.t.max integer value, maximum days between two recession data
#' @param tgood integer value, minimum length of the recession in days
#' @param Nmin.rec integer value, minimum number of data in a recession curve
#' @param tburn.rec numeric, burn factor, >=0 and <1. 0.4 means the first 40 percent of the recession samples are discarded.
#'
#' @return List with the following components :
#' \enumerate{
#'    \item Rec_extracted: data frame, all information about stage-recession
#'    \item summary_rec: data frame, summary with time recession and number of stage data by recession
#' }
#' @export
#'
#' @examples
#'
#' recessions=Extraction_recession(H=ArdecheRiverMeyrasStage$H,
#'                                 uH=0.05,
#'                                 time=ArdecheRiverMeyrasStage$Date,
#'                                 chi=0.5,
#'                                 tgood=30)
#'
#' # Check information of recession extracted
#' recessions$Rec_extracted
#'
#' # Check summary table of recession extracted
#' recessions$summary_rec
#'
#' # Plot all recession extracted
#' plot_rec_extracted(Rec_extracted = recessions$Rec_extracted)
Extraction_recession <- function(H,
                                 uH,
                                 time,
                                 filter.H=1,
                                 chi=stats::quantile(H,probs = 0.95,na.rm = TRUE),
                                 delta.t.min=0,
                                 delta.t.max=20,
                                 tgood=10,
                                 Nmin.rec=10,
                                 tburn.rec=0){

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

  # Find all recession curves :

  # initialize:
  hrec = uHrec = data_rec = indx = c()
  trec =  as.POSIXct(character(0))

  hrec[1]  = stage.record.post$h[1]
  trec[1]  = stage.record.post$t[1]
  uHrec[1] = stage.record.post$uH[1]
  indx[1]  = 1
  j = 2  # run all recession stages one by one
  m = 1  # stage accepted (descending point) for a specific recession
  k = 1 # Number of the recession

  j_max = length(stage.record.post$t)

  while (j <= j_max) {

    diff_t=as.numeric(difftime(stage.record.post$t[j],
                    trec[m],
                    units='days'))

    # Check 1: minimum distance between recession data (delta.t.min)
    # minimum distance in case of too
    # many data to process and maximum distance to avoid long periods of no data which may create
    # discontinuous recessions
    if(diff_t < delta.t.min){
      j <- j + 1
      next
    }

    # Check 2: maximum distance between recession data (delta.t.max)
    # Consider different recessions if gap in time is higher than delta.t.max.
    # Handle the missing data in the stage record.
    if(diff_t > delta.t.max){
      # change number of recession (k)
      k = k+1
      m=m+1

      hrec[m]  = stage.record.post$h[j]
      trec[m]  = stage.record.post$t[j]
      uHrec[m] = stage.record.post$uH[j]
      indx[m]  = k
      j <- j + 1
      next
    }

    # Check 3: Delta h is negative ->  potential stage-recession data to save
    # track delta h from two data : stage value - last stage-recession data stored
    diff_h=stage.record.post$h[j]-hrec[m]
    if(diff_h < 0){
      m=m+1

      hrec[m]  = stage.record.post$h[j]
      trec[m]  = stage.record.post$t[j]
      uHrec[m] = stage.record.post$uH[j]
      indx[m]  = k

    }else{ # Case when delta h is positive

      # Check 4: recession is already finished ? (chi)
      # chi is a threshold that defines the end of one recession and the beginning of another
      if(diff_h >= chi){
        m=m+1
        k = k+1

        hrec[m]  = stage.record.post$h[j]
        trec[m]  = stage.record.post$t[j]
        uHrec[m] = stage.record.post$uH[j]
        indx[m]  = k

      }
      # if difference is not higher than chi, do not take into account this data
      # and keep going. So j is increasing and m remains the same value
    }
    j <- j + 1
  }

  data_df=data.frame(hrec=hrec,
                     date=trec,
                     uHrec=uHrec,
                     indx=indx)

  # Summary table :
  summary_temp <- data.frame(data_df %>%
                               group_by(indx) %>%
                               summarise(
                                 first_time = first(date),  # First date for each group
                                 last_time = last(date),    # Last date for each group
                                 diff_time_extremes = as.numeric(difftime(last_time, first_time, units = "days")),  # Difference in hours
                                 count_values = n()  # Total number of values for each group
                                 ))
  # Check 4 : tgood
  # Select only recessions with length > tgood
  rm_indx_tgood = which(summary_temp$diff_time_extremes < tgood)
  # Check 5 : Nmin
  # Select only recessions with a number of point > Nmin :
  rm_indx_Nmin = which(summary_temp$count_values < Nmin.rec)

  rm_indx=unique(c(rm_indx_tgood,rm_indx_Nmin))

  data_df_filtered <- data_df[!data_df$indx %in% rm_indx, ]

  # Re-assign the indx of the recession
  # Same data frame because a colonne will be modify and the use for burning step
  data_rec_unburned=data_df_filtered

  # Step 1: Get unique values of 'indx' and sort them
  unique_indx <- unique(data_df_filtered$indx)

  # Step 2: Create a mapping from old 'indx' values to new sequential values
  indx_map <- data.frame(old = unique_indx, new = seq_along(unique_indx))

  # Step 3: Replace 'indx' values in the filtered data frame with new sequential values
  data_rec_unburned$indx <- indx_map$new[match(data_df_filtered$indx, indx_map$old)]


  if(tburn.rec!=0){
    # Apply burn of the first part of the recession to reduce data information
    data_rec_burned <- data.frame(data_rec_unburned %>%
                             group_by(indx) %>%
                             slice(-(1:floor(n() * tburn.rec))))

  }else{
    data_rec_burned=data_rec_unburned
  }

  data_rec = data.frame(data_rec_burned %>%
                           group_by(indx) %>%
                           mutate(
                             time_rec = as.numeric(difftime(date, first(date), units = "days")),  # Difference in hours
                            ))

  # Summary table
  summary_rec <- data.frame(data_rec_burned %>%
                               group_by(indx) %>%
                               summarise(
                                 diff_time_extremes = as.numeric(difftime( last(date), first(date), units = "days")),  # Difference in hours
                                 count_data = n()  # Total number of values for each group
                               ))



  return(list(Rec_extracted = data_rec,
              summary_rec=summary_rec))
}

#' Model and segmentation of recession
#'
#' Modelling recession using a catalog of fit models available to apply segmentation procedure.
#' Please see details for understaing how to use this function
#'
#' @param time_rec real vector, recession duration relative to the first data detected during the recession
#' @param daterec vector, time in POSIXct format
#' @param hrec real vector, stage value of the recessions
#' @param uHrec real vector, uncertainty of stage
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
#' @return  List with the following components :
#' \enumerate{
#'   \item summary_results : list with this sub components :
#'   \itemize{
#'
#'      \item summary: list, summarize the information to present to the user
#'      \itemize{
#'          \item data: data frame, all data of (H,Q and uQ) with their respective periods after segmentation
#'          \item shift: data frame, all detected shift time
#'      }
#'      \item plot : list, data formatted to use as input for some plot functions
#'      \item res: list, provide all the information of the periods from tree structure
#'      \itemize{
#'          \item tau: real vector, estimated shift times
#'          \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'          \item mcmc: data frame, MCMC simulation
#'          \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'          \item DIC: real, DIC estimation
#'          \item nS: integer, optimal number of segments following DIC criterion
#'        }
#'      \item tree: data frame, provide tree structure
#'      \item origin.date: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#'   }
#'   \item residuals_all_data : list, provide all information of recession data observed and simulated
#'   \itemize{
#'      \item time_rec: real value, time of the recession relative to first data of each recession
#'      \item hrec: real value, stage observed
#'      \item H_rec_sim: real value, stage simulated
#'      \item H_res: real value, residual between stage observed and simulated
#'      \item uhrec: real value, uncertainty on stage observed
#'      \item uH_sim: real value, uncertainty on stage simulated
#'   }
#'   \item parameters : all information about the parameters estimated depending of the fit model used
#' }
#' @export
#' @importFrom stringr str_detect
#' @importFrom utils read.table
#' @import dplyr
#' @details
#' To get the catalog of models available using `GetCatalog()`, all functions specifying Recession in the fit are concerned.
#' For more information about recession model and their default values, please go to the fit function, e.g. `?fitRecession_M3`
#' If prior knowledge about recession parameters is given, be sure to give dataset too.
#' The direction file path for this dataset must be `file.path(tempdir(),'BaM','Recession')` as shown in the example
#' Please respect the name of the parameters if default values will not be used. e.g. beta paremeter's will be "beta"
#'
#' @examples
#' recessions_extracted=Extraction_recession(H=ArdecheRiverMeyrasStage$H,
#'                                           uH=0.05,
#'                                           time=ArdecheRiverMeyrasStage$Date,
#'                                           chi=0.5,
#'                                           tgood=30)
#'
#' recessions = recessions_extracted$Rec_extracted
#'
#' # Plot all recession extracted
#' plot_rec_extracted(Rec_extracted = recessions)
#'
#' D=RBaM::dataset(X=recessions['time_rec'],
#'                 Y=recessions['hrec'],
#'                 Yu=recessions['uHrec'],
#'                 VAR.indx=recessions['indx'],
#'                 data.dir=file.path(tempdir(),'BaM','Recession'))
#'
#' # Choose fit recession model
#' fit=fitRecession_M3
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
#'                                                     daterec=recessions$date,
#'                                                     hrec=recessions$hrec,
#'                                                     uHrec=recessions$uHrec,
#'                                                     indx=recessions$indx,
#'                                                     alpha1.object=alpha1.object,
#'                                                     alpha2.object=alpha2.object,
#'                                                     beta.object=beta.object,
#'                                                     Dataset.object=D,
#'                                                     nCyclesrec=1,
#'                                                     burnrec=0.1,
#'                                                     nSlimrec=5,
#'                                                     funk=fit)
#'
#' # Plot recession segmentation
#' plot_modelAndSegm_recession(model_rec=model_rec,
#'                             spec_recession=c(2,16,28),
#'                             fit=fitRecession_M3,
#'                             equation_rec='Recession_M3_Equation')
#'
#' plotSegmentation(summary = model_rec$summary_results$summary,
#'                  plot_summary = model_rec$summary_results$plot )
ModelAndSegmentation.recession.regression <- function(time_rec,
                                                      daterec,
                                                      hrec,
                                                      uHrec,
                                                      indx,
                                                      nSmax=3,
                                                      nMin= 1,
                                                      nCycles=100,
                                                      burn=0.5,
                                                      nSlim=max(nCycles/10,1),
                                                      temp.folder=file.path(tempdir(),'BaM'),
                                                      funk=fitRecession_M3,...){
  # Check information given in input
  if(any(time_rec<0))stop('time_rec must be positive')
  if(any(uHrec<0))stop('uHrec must be positive')
  if(any(indx<0))stop('indx must be positive')
  if(any(!lubridate::is.POSIXct(daterec)))stop('daterec must be in Posixct format')
  if(indx[1]!=1)stop('Fist number of indx must be one to start the sequence')
  if(any(diff(unique(indx))!=1))stop('indx must be a sequence of consecutive numbers')
  if(any(is.na(time_rec)) | any(is.na(hrec)) | any(is.na(uHrec)) | any(is.na(indx))){
    stop('Missing values not allowed in time, stage, uncertainty or index')
  }
  check <- check_vector_lengths(time_rec,hrec,uHrec,indx)
  if(is.null(check)){
    stop('time,stage, uncertainty or index do not have the same length')
  }

  DF.order <- data.frame(time_rec=time_rec,
                         daterec=daterec,
                         hrec=hrec,
                         uHrec=uHrec,
                         indx=indx)

  # Recession modelling using funk function user-defined:
  rec_model=funk(time_rec=DF.order$time_rec,
                 hrec=DF.order$hrec,
                 uHrec=DF.order$uHrec,
                 indx=DF.order$indx,
                 ...)

  # Residuals of recession observed and simulated
  residualsData <- rec_model$ResultsResiduals
  # Asymptotic stage:
  h_asymptotic <- rec_model$parameters$beta
  # date of last recession data recorded
  date_asymptotic_df= data.frame(DF.order %>%
                                group_by(indx)%>%
                                summarize(
                                  date_asymp=max(daterec)
                                ))
  date_asymptotic=date_asymptotic_df$date_asymp
  # Read summary
  summary.MCMC   <- utils::read.table(file= file.path(temp.folder,'Recession',"Results_Summary.txt"),header=TRUE)
  # Focus in MaxPost and Standard deviation of beta
  # Discuss which uncertainty will we take into account
  u_asymptotic=c(t(summary.MCMC[which(rownames(summary.MCMC)==c('St.Dev.')),
                                which(stringr::str_detect(colnames(summary.MCMC),pattern='beta'))]))
  # Run recursive segmentation
  results=recursive.segmentation(obs=h_asymptotic,
                                 time=date_asymptotic,
                                 u=u_asymptotic,
                                 nSmax=nSmax,
                                 nMin=nMin,
                                 nCycles=nCycles,
                                 burn=burn,
                                 nSlim=nSlim,
                                 temp.folder=temp.folder)

  return(list(summary_results=results,
              residuals_all_data=residualsData,
              parameters=rec_model$parameters))

}
