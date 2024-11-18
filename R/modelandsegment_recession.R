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
#' recessions=Extract_Recessions(H=ArdecheRiverStage$H,
#'                                 uH=0.05,
#'                                 time=ArdecheRiverStage$Date,
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
#' PlotExtract_Recessions(Rec_extracted = recessions$Rec_extracted)
Extract_Recessions <- function(H,
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

#' Criticize recession modeling
#'
#' Weighted RMSE is used as indicator to qualify the performance of recession modeling.
#'
#' @param raw_data data frame, all stage record
#' @param data.object object, model object output by the RBaM function
#' @param residuals data frame, residual file obtained using RBaM function
#' @param summary_MCMC data frame, MCMC summary file obtained using RBaM function
#' @param rmse_threshold real positive value, RMSE threshold to qualify a modeled recession as good
#' @param min.weight real positive value, weight given to the first part of the data recession
#' @param max.weight real positive value, weight given to the last part of the data recession
#'
#' @return list with the following components :
#' \enumerate{
#'    \item df.to.segm: data frame, description below:
#'      \itemize{
#'          \item indx: integer value, recession curve index
#'          \item b_estimated: real value, asymptotic height estimation
#'          \item ub_estimated: real value, uncertainty in asymptotic height estimation expressed as a standard deviation
#'          \item date: date, time in POSIXct format. time at which the recession begins
#'      }
#'    \item rec_data.plot.h.dt: data frame, description below:
#'      \itemize{
#'          \item time_rec: real positive value, time recession
#'          \item datarec: date, time in POSIXct format. time of each stage-recession
#'          \item hrec: real value, stage-recession observed
#'          \item uHrec: real value, uncertainty in stage-recession observed expressed as a standard deviation
#'          \item indx: integer value, common index for all stage-recession data of a same recession curve
#'          \item status: string, after using RMSE weighted, recession modeled will be accepted or rejected
#'      }
#'    \item residuals.all.info: data frame, description below:
#'      \itemize{
#'          \item indx: real value, recession curve index
#'          \item X1_obs: real value, duration of recession observed
#'          \item Y1_obs: real value, stage-recession observed
#'          \item Y1_sim: real value, stage-recession simulated
#'          \item uH_obs: real positive value, uncertainty on the stage-recession observed
#'          \item weight: real positive value, weight assigned depending on the recession duration
#'          \item value_rmse: real positive value, RMSE calculation
#'          \item rmse_threshold: real positive value, threshold user-defined to accept or reject some recession simulations
#'          \item status: string, after using RMSE weighted, recession modeled will be accepted or rejected
#'      }
#' }
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_extract
Criticize_Recessions <- function(raw_data,
                                 data.object,
                                 residuals,
                                 summary_MCMC,
                                 rmse_threshold = 0.3,
                                 min.weight = 0.1,
                                 max.weight = 1){


  if(max.weight < min.weight)stop('max.weight must be higher than min.weight')

  rec_data=data.object$data

  residuals.indx=data.frame(residuals['X1_obs'],
                            residuals['Y1_obs'],
                            residuals['Y1_sim'],
                            uH_obs=rec_data$uHrec,
                            indx=rec_data$indx)

  residuals.indx$weight = NA

  # Mean max value of all recession
  max.time.rec.local =  data.frame(rec_data %>%
                                     dplyr::group_by(indx)%>%
                                     dplyr::summarize(max=max(time_rec)))

  mean.max.val <- mean(max.time.rec.local$max)

  Ncurves=max(rec_data$indx)

  # Weighting for each recession data item to give greater importance to the long duration of the recession
  for(i in 1:Ncurves){
    rec.duration = rec_data$time_rec[which(rec_data$indx==i)]

    min.val <- min(rec.duration)
    max.val <- max(rec.duration)
    if(max.val < mean.max.val){
      max.val = mean.max.val
    }
    # Calculate normalized weights
    weights <- min.weight + (rec.duration - min.val) / (max.val - min.val) * (max.weight - min.weight)

    residuals.indx$weight[which(residuals.indx$indx==i)] <- weights
  }

  # RMSE is in the same units that input data -> meters
  rmse_curve=data.frame(residuals.indx %>%
                          dplyr::group_by(indx) %>%
                          dplyr::summarize(value_rmse=RMSE_Weighted(Y1_obs,Y1_sim,weight)))

  # criticize recession modeled curve by curve


  table_comparison = cbind(rmse_curve,rmse_threshold)

  recession_accepted = table_comparison$indx[which(table_comparison$value_rmse <= table_comparison$rmse_threshold)]

  rejected.curves = unique(rec_data$indx)[!unique(rec_data$indx) %in% recession_accepted]

  # Detect curve to reject in summary_MCMC
  if(length(rejected.curves)!=0){
    summary_MCMC=summary_MCMC[,-match(paste0('beta_k_',rejected.curves),colnames(summary_MCMC))]
  }

  # residuals.all.info helps to plot all recession in a same plot (simulated and observed and a status accepted or jerected)
  # Paste rmse calculation
  residuals.all.info.p=merge(residuals.indx,rmse_curve,by='indx')

  residuals.all.info=data.frame(residuals.all.info.p,
                                rmse_threshold=rmse_threshold,
                                status=NA)

  residuals.all.info$status[residuals.all.info$indx %in% rejected.curves] = "Rejected"
  residuals.all.info$status[residuals.all.info$indx %in% recession_accepted] = "Accepted"

  # Focus in asymptotic height (b) estimated :
  indx_MAP=which(c('MaxPost')==rownames(summary_MCMC))
  indx_st_dev=which(c('St.Dev.')==rownames(summary_MCMC))

  table.asymp.height=(summary_MCMC[c(indx_MAP,indx_st_dev),grep(paste0('beta_k_'),colnames(summary_MCMC))])


  b_obs=data.frame(table.asymp.height ['MaxPost',] %>%
                     tidyr::pivot_longer(
                       cols = dplyr::starts_with("beta_k"),   # Pivot all columns
                       names_to = "b",                        # Create 'b' column for each beta variable
                       values_to = "b_estimated"                    # Store the values in 'value'
                     )) %>%
    dplyr::mutate(indx = stringr::str_extract(b, "\\d+$") %>% as.numeric()) %>%
    dplyr::select(-b)

  ub_obs=data.frame(table.asymp.height ['St.Dev.',] %>%
                      tidyr::pivot_longer(
                        cols = dplyr::starts_with("beta_k"),   # Pivot all columns
                        names_to = "b",                        # Create 'b' column for each beta variable
                        values_to = "ub_estimated"                    # Store the values in 'value'
                      )) %>%
    dplyr::mutate(indx = stringr::str_extract(b, "\\d+$") %>% as.numeric()) %>%
    dplyr::select(-b)

  # Get recession data h(t) for plotting
  rec_data.plot.h.dt = data.frame(raw_data, status=NA)
  rec_data.plot.h.dt$status[raw_data$indx %in% rejected.curves] = "Rejected"
  rec_data.plot.h.dt$status[raw_data$indx %in% recession_accepted] = "Accepted"

  # Get data frame to segment
  df.to.segm.p=merge(b_obs,ub_obs)

  begin.date.record=data.frame(raw_data %>%
                                 group_by(indx) %>%
                                 summarize(date=min(daterec)))

  df.to.segm=merge(df.to.segm.p,begin.date.record)

  return(list(df.to.segm=df.to.segm,
              rec_data.plot.h.dt=rec_data.plot.h.dt,
              residuals.all.info=residuals.all.info))
}


#' Model and segmentation of recession
#'
#' Modelling recession using a catalog of fit models available to apply segmentation procedure.
#' Please see details for understanding how to use this function
#'
#' @param time_rec real vector, recession duration relative to the first data detected during the recession
#' @param daterec vector, time in POSIXct format mandatory
#' @param hrec real vector, stage value of the recessions
#' @param uHrec real vector, uncertainty of stage
#' @param indx  integer, factor used to gather the data of a same recession
#' @param nSmax integer, maximum number of segments to assess
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations of the segmentation
#' @param temp.folder.Recession directory, temporary directory to write computations of the recession
#' @param funk string, model for estimating the recession
#' @param ... optional arguments to funk
#' @param RMSE_args list, optional arguments to pass in `Criticize_Recessions` function for setting weighted RMSE
#'
#' @return  List with the following components:
#' \enumerate{
#'   \item summary.rec.extracted : list with this sub components:
#'      \itemize{
#'          \item data: data frame, all data of recession (time, date, H and uH) with their respective period after segmentation and status if the recession model was accepted or rejected
#'          \item shift: data frame, all detected shift time
#'          }
#'   \item summary.results.segm: list with this sub components:
#'      \itemize{
#'          \item data: data frame, asymptotic height estimation by curve accepted with their respective period after segmentation
#'          \item shift: data frame, all detected shift time
#'      }
#'   \item summary.residual: data frame, the residual information between the MAP simulation and the observation with weighted RMSE to criticise the recession estimate and finally the acceptance or rejection status.
#'   \item plots: list, data formatted to use as input for some plot functions
#'   \item res: list, provide all the information of the periods from tree structure
#'      \itemize{
#'          \item tau: real vector, estimated shift times
#'          \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'          \item mcmc: data frame, MCMC simulation
#'          \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'          \item DIC: real, DIC estimation
#'          \item nS: integer, optimal number of segments following DIC criterion
#'      }
#'    \item tree: data frame, provide tree structure
#'    \item origin.date: positive real or date, date describing origin of the segmentation for a sample
#' }
#' @export
#' @import dplyr
#' @importFrom stringr str_detect
#' @importFrom utils read.table
#' @importFrom RBaM dataset
#'
#' @details
#' To get the catalog of available models, including equation and estimation, please see: `?GetCatalog_Recession`.
#' The funk input must be one of the available models in: `names(GetCatalog_Recession())`
#'
#' Be careful when specifying temporal folder for computations.
#' temp.folder must contain temp.folder.Recession as the default value:
#' temp.folder= `file.path(tempdir(),'BaM')`
#' temp.folder.Recession= `file.path(tempdir(),'BaM','Recession')`
#'
#' @examples
#' recessions_extracted=Extract_Recessions(H=ArdecheRiverStage$H,
#'                                           uH=0.05,
#'                                           time=ArdecheRiverStage$Date,
#'                                           chi=0.5,
#'                                           tgood=30)
#'
#' recessions = recessions_extracted$Rec_extracted
#'
#' # Plot all recession extracted
#' PlotExtract_Recessions(Rec_extracted = recessions)
#'
#' # Choose fit recession model
#' fit='BR1'
#'
#' model_rec=ModelAndSegmentation_Recessions(time_rec=recessions$time_rec,
#'                                                     daterec=recessions$date,
#'                                                     hrec=recessions$hrec,
#'                                                     uHrec=recessions$uHrec,
#'                                                     indx=recessions$indx,
#'                                                     funk=fit,
#'                                                     RMSE_args=list(rmse_threshold =0.3))
#'
#' # Plot recession segmentation user-defined
#' PlotSimObs_Recessions (model_rec=model_rec,
#'                        spec_recession=c(2,16,28,48),
#'                        recession_rejected=TRUE,
#'                        all_recession=FALSE)
#'
#' # Plot the shift times detected
#' PlotSegmentation_Recessions(model_rec=model_rec)
#'
#' # Plot using all data of stage record
#' PlotSegmentation_Recessions_Hdt(time=ArdecheRiverStage$Date,
#'                                 obs=ArdecheRiverStage$H,
#'                                 u=0.05,
#'                                 plot_summary = model_rec$plot)
ModelAndSegmentation_Recessions <- function(time_rec,
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
                                            temp.folder.Recession= file.path(tempdir(),'BaM','Recession'),
                                            funk='BR1',...,
                                            RMSE_args=list()){

  id_recession_model=which(funk==names(GetCatalog_Recession()))
  # Check conditions if they are satisfied
  invisible(check_recession_modeling(time_rec=time_rec,
                                     daterec=daterec,
                                     hrec=hrec,
                                     uHrec=uHrec,
                                     indx=indx,
                                     id_recession_model=id_recession_model))

  # Get equation of the chosen model
  equation_recession_model=GetCatalog_Recession()[[id_recession_model]]$Equation
  # Get estimation function of the chosen model
  recession_model=GetCatalog_Recession()[[id_recession_model]]$funk

  DF.order <- data.frame(time_rec=time_rec,
                         daterec=daterec,
                         hrec=hrec,
                         uHrec=uHrec,
                         indx=indx)
  # dataset object
  D = RBaM::dataset(X=DF.order['time_rec'],
                    Y=DF.order['hrec'],
                    Yu=DF.order['uHrec'],
                    VAR.indx=DF.order['indx'],
                    data.dir= temp.folder.Recession)

  # For detecting and estimating shift times, there are three step:
  # This function is an optional run with default inputs. However, the function
  # could be used separately to modify the inputs according to the user's requirements

  # 1. Recession estimation:
  list.rec.est=recession_model(data.object=D,
                               equation_rec_model=equation_recession_model(),
                               temp.folder.Recession=temp.folder.Recession,
                               ...)

  # 2. Criticize of recession modeling (RMSE weighted)
  list.rec.criticize <- do.call(Criticize_Recessions,
                                c(list(raw_data = DF.order,
                                       data.object = D,
                                       residuals = list.rec.est$residuals,
                                       summary_MCMC = list.rec.est$summary_MCMC),
                                  RMSE_args))

  # 3. Segmentation of asymptotic height
  # Data to segment
  df.to.segm=list.rec.criticize$df.to.segm

  # Run recursive segmentation
  results=Recursive_Segmentation(obs=df.to.segm$b_estimated,
                                 time=df.to.segm$date,
                                 u=df.to.segm$ub_estimated,
                                 nSmax=nSmax,
                                 nMin=nMin,
                                 nCycles=nCycles,
                                 burn=burn,
                                 nSlim=nSlim,
                                 temp.folder=temp.folder)

  # return data recession with period after segmentation
  results$summary$data$indx =  df.to.segm$indx

  period.indx=results$summary$data[,c('period','indx')]

  input.data.with.period=merge( list.rec.criticize$rec_data.plot.h.dt,period.indx, by='indx')
  summary.rec.extracted = list(data = input.data.with.period,
                               shift = results$summary$shift)

  return(list(summary.rec.extracted=summary.rec.extracted,
              summary.results.segm=results$summary,
              summary.residual = list.rec.criticize$residuals.all.info,
              plots=results$plot,
              res=results$res,
              tree=results$tree,
              origin.date=results$origin.date,
              n.rec.max=max(DF.order$indx)))

}
