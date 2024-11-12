
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
criticize_rec_model <- function(raw_data,
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
                          dplyr::summarize(value_rmse=rmse.weighted(Y1_obs,Y1_sim,weight)))

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
