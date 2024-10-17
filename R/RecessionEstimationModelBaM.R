#' Recession model with two exponential and asymptotic
#'
#' @param CalData data frame, calibration data used for calibration
#' @param Ncurves real positive value, number of recession extracted
#' @param M object, model object output by the RBaM function
#' @param MCMC data frame, MCMC results
#' @param t_grid data frame, recession duration
#' @param temp.folder.Recession directory, temporary directory to write computations
#' @param remnant_prior object, remnant error model output by the RBaM function
#'
#' @return list, recession estimation with their uncertainties
Estimation_Recession_M3 <- function(CalData,
                                  Ncurves,
                                  M,
                                  MCMC,
                                  t_grid,
                                  temp.folder.Recession,
                                  remnant_prior){

  allData=c()
  j=1

  for( i in Ncurves){
    alpha1_nonVAR=RBaM::parameter(name=M$par[[1]]$name,
                                  init=M$par[[1]]$init[i],
                                  prior.dist=M$par[[1]]$prior[[i]]$dist,
                                  prior.par=M$par[[1]]$prior[[i]]$par)

    alpha2_nonVAR=RBaM::parameter(name=M$par[[3]]$name,
                                  init=M$par[[3]]$init[i],
                                  prior.dist=M$par[[3]]$prior[[i]]$dist,
                                  prior.par=M$par[[3]]$prior[[i]]$par)

    beta_nonVAR=RBaM::parameter(name=M$par[[5]]$name,
                                init=M$par[[5]]$init[i],
                                prior.dist=M$par[[5]]$prior[[i]]$dist,
                                prior.par=M$par[[5]]$prior[[i]]$par)

    lambda1=RBaM::parameter(name=M$par[[2]]$name,
                            init=M$par[[2]]$init,
                            prior.dist=M$par[[2]]$prior$dist,
                            prior.par=M$par[[2]]$prior$par)

    lambda2=RBaM::parameter(name=M$par[[4]]$name,
                            init=M$par[[4]]$init,
                            prior.dist=M$par[[4]]$prior$dist,
                            prior.par=M$par[[4]]$prior$par)

    M_nonVAR=RBaM::model(ID='Recession_h',
                         nX=1,
                         nY=1,
                         par=list(alpha1_nonVAR,
                                  lambda1,
                                  alpha2_nonVAR,
                                  lambda2,
                                  beta_nonVAR))

    # Columns names to use during prediction
    columns_indx=c(paste0('alpha1_',i),
                   'lambda1',
                   paste0('alpha2_',i),
                   'lambda2',
                   paste0('beta_',i),
                   'Y1_gamma1',
                   'Y1_gamma2',
                   'LogPost')

    parSamples=MCMC[columns_indx]

    # Define a 'prediction' object for total predictive uncertainty
    totalU=RBaM::prediction(X=t_grid, # stage values
                            spagFiles='H_rec_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.Recession, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE, # propagate structural uncertainty ?
                            parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
    paramU=RBaM::prediction(X=t_grid,
                            spagFiles='H_rec_ParamU.spag',
                            data.dir=temp.folder.Recession,
                            doParametric=TRUE,
                            doStructural=FALSE,
                            parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' RC maximizing the posterior density
    maxpost=RBaM::prediction(X=t_grid,
                             spagFiles='H_rec_Maxpost.spag',
                             data.dir=temp.folder.Recession,
                             doParametric=FALSE,
                             doStructural=FALSE,
                             parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    RBaM:: BaM(mod=M_nonVAR,
               data=D,
               workspace = temp.folder.Recession,
               dir.exe = file.path(find.package("RBaM"), "bin"),
               pred=list(totalU,paramU,maxpost), # list of predictions
               remnant = remnant_prior,
               doCalib=FALSE,
               doPred=TRUE)

    # Get total uncertainty from the rating curve
    TotalUenv=read.table(file.path(temp.folder.Recession,'H_rec_TotalU.env'),header = TRUE)

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(TotalUenv))
    TotalUenv[cols_to_convert] <- suppressWarnings(lapply(TotalUenv[cols_to_convert],
                                                          function(x) as.numeric(x)))

    # Add parametric uncertainty
    ParametricUenv=read.table(file.path(temp.folder.Recession,'H_rec_ParamU.env'),header = TRUE)

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(ParametricUenv))
    ParametricUenv[cols_to_convert] <- suppressWarnings(lapply(ParametricUenv[cols_to_convert],
                                                               function(x) as.numeric(x)))

    # Add maxpost rating curve
    MAPREC=read.table(file.path(temp.folder.Recession,'H_rec_Maxpost.spag'))

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(MAPREC))
    MAPREC[cols_to_convert] <- suppressWarnings(lapply(MAPREC[cols_to_convert],
                                                       function(x) as.numeric(x)))

    allData[[j]]=list(CalData=CalData,
                      TotalUenv=TotalUenv,
                      ParametricUenv=ParametricUenv,
                      MAPREC=MAPREC,
                      HgridPlot=t_grid)
    j=j+1
  }
  return(allData)
}

#' Exponential recession model with five parameters
#'
#' @param CalData data frame, calibration data used for calibration
#' @param Ncurves real positive value, number of recession extracted
#' @param M object, model object output by the RBaM function
#' @param MCMC data frame, MCMC results
#' @param t_grid data frame, recession duration
#' @param temp.folder.Recession directory, temporary directory to write computations
#' @param remnant_prior object, remnant error model output by the RBaM function
#'
#' @return list, recession estimation with their uncertainties
Estimation_Recession_BR1 <- function(CalData,
                                    Ncurves,
                                    M,
                                    MCMC,
                                    t_grid,
                                    temp.folder.Recession,
                                    remnant_prior){

  allData=c()
  j=1

  for( i in Ncurves){
    alpha1=RBaM::parameter(name=M$par[[1]]$name,
                            init=M$par[[1]]$init,
                            prior.dist=M$par[[1]]$prior$dist,
                            prior.par=M$par[[1]]$prior$par)

    lambda1=RBaM::parameter(name=M$par[[2]]$name,
                            init=M$par[[2]]$init,
                            prior.dist=M$par[[2]]$prior$dist,
                            prior.par=M$par[[2]]$prior$par)

    lambda2=RBaM::parameter(name=M$par[[3]]$name,
                            init=M$par[[3]]$init,
                            prior.dist=M$par[[3]]$prior$dist,
                            prior.par=M$par[[3]]$prior$par)

    lambda3=RBaM::parameter(name=M$par[[4]]$name,
                            init=M$par[[4]]$init,
                            prior.dist=M$par[[4]]$prior$dist,
                            prior.par=M$par[[4]]$prior$par)

    beta_nonVAR=RBaM::parameter(name=M$par[[5]]$name,
                                init=M$par[[5]]$init[i],
                                prior.dist=M$par[[5]]$prior[[i]]$dist,
                                prior.par=M$par[[5]]$prior[[i]]$par)


    M_nonVAR=RBaM::model(ID='Recession_h',
                         nX=1,
                         nY=1,
                         par=list(alpha1,
                                  lambda1,
                                  lambda2,
                                  lambda3,
                                  beta_nonVAR))

    # Columns names to use during prediction
    columns_indx=c('alpha1',
                   'lambda1',
                   'lambda2',
                   'lambda3',
                   paste0('beta_',i),
                   'Y1_gamma1',
                   'Y1_gamma2',
                   'LogPost')

    parSamples=MCMC[columns_indx]

    # Define a 'prediction' object for total predictive uncertainty
    totalU=RBaM::prediction(X=t_grid, # stage values
                            spagFiles='H_rec_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.Recession, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE, # propagate structural uncertainty ?
                            parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
    paramU=RBaM::prediction(X=t_grid,
                            spagFiles='H_rec_ParamU.spag',
                            data.dir=temp.folder.Recession,
                            doParametric=TRUE,
                            doStructural=FALSE,
                            parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' RC maximizing the posterior density
    maxpost=RBaM::prediction(X=t_grid,
                             spagFiles='H_rec_Maxpost.spag',
                             data.dir=temp.folder.Recession,
                             doParametric=FALSE,
                             doStructural=FALSE,
                             parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    RBaM:: BaM(mod=M_nonVAR,
               data=D,
               workspace = temp.folder.Recession,
               dir.exe = file.path(find.package("RBaM"), "bin"),
               pred=list(totalU,paramU,maxpost), # list of predictions
               remnant = remnant_prior,
               doCalib=FALSE,
               doPred=TRUE)

    # Get total uncertainty from the rating curve
    TotalUenv=read.table(file.path(temp.folder.Recession,'H_rec_TotalU.env'),header = TRUE)

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(TotalUenv))
    TotalUenv[cols_to_convert] <- suppressWarnings(lapply(TotalUenv[cols_to_convert],
                                                          function(x) as.numeric(x)))

    # Add parametric uncertainty
    ParametricUenv=read.table(file.path(temp.folder.Recession,'H_rec_ParamU.env'),header = TRUE)

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(ParametricUenv))
    ParametricUenv[cols_to_convert] <- suppressWarnings(lapply(ParametricUenv[cols_to_convert],
                                                               function(x) as.numeric(x)))

    # Add maxpost rating curve
    MAPREC=read.table(file.path(temp.folder.Recession,'H_rec_Maxpost.spag'))

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(MAPREC))
    MAPREC[cols_to_convert] <- suppressWarnings(lapply(MAPREC[cols_to_convert],
                                                       function(x) as.numeric(x)))

    allData[[j]]=list(CalData=CalData,
                      TotalUenv=TotalUenv,
                      ParametricUenv=ParametricUenv,
                      MAPREC=MAPREC,
                      HgridPlot=t_grid)
    j=j+1
  }
  return(allData)
}

#' Exponential recession model with four parameters
#'
#' @param CalData data frame, calibration data used for calibration
#' @param Ncurves real positive value, number of recession extracted
#' @param M object, model object output by the RBaM function
#' @param MCMC data frame, MCMC results
#' @param t_grid data frame, recession duration
#' @param temp.folder.Recession directory, temporary directory to write computations
#' @param remnant_prior object, remnant error model output by the RBaM function
#'
#' @return list, recession estimation with their uncertainties
Estimation_Recession_BR2 <- function(CalData,
                                     Ncurves,
                                     M,
                                     MCMC,
                                     t_grid,
                                     temp.folder.Recession,
                                     remnant_prior){

  allData=c()
  j=1

  for( i in Ncurves){
    alpha1=RBaM::parameter(name=M$par[[1]]$name,
                           init=M$par[[1]]$init,
                           prior.dist=M$par[[1]]$prior$dist,
                           prior.par=M$par[[1]]$prior$par)

    lambda1=RBaM::parameter(name=M$par[[2]]$name,
                            init=M$par[[2]]$init,
                            prior.dist=M$par[[2]]$prior$dist,
                            prior.par=M$par[[2]]$prior$par)

    c1=RBaM::parameter(name=M$par[[3]]$name,
                       init=M$par[[3]]$init,
                       prior.dist=M$par[[3]]$prior$dist,
                       prior.par=M$par[[3]]$prior$par)

    beta_nonVAR=RBaM::parameter(name=M$par[[4]]$name,
                                init=M$par[[4]]$init[i],
                                prior.dist=M$par[[4]]$prior[[i]]$dist,
                                prior.par=M$par[[4]]$prior[[i]]$par)


    M_nonVAR=RBaM::model(ID='Recession_h',
                         nX=1,
                         nY=1,
                         par=list(alpha1,
                                  lambda1,
                                  c1,
                                  beta_nonVAR))

    # Columns names to use during prediction
    columns_indx=c('alpha1',
                   'lambda1',
                   'c1',
                   paste0('beta_',i),
                   'Y1_gamma1',
                   'Y1_gamma2',
                   'LogPost')

    parSamples=MCMC[columns_indx]

    # Define a 'prediction' object for total predictive uncertainty
    totalU=RBaM::prediction(X=t_grid, # stage values
                            spagFiles='H_rec_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.Recession, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE, # propagate structural uncertainty ?
                            parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
    paramU=RBaM::prediction(X=t_grid,
                            spagFiles='H_rec_ParamU.spag',
                            data.dir=temp.folder.Recession,
                            doParametric=TRUE,
                            doStructural=FALSE,
                            parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' RC maximizing the posterior density
    maxpost=RBaM::prediction(X=t_grid,
                             spagFiles='H_rec_Maxpost.spag',
                             data.dir=temp.folder.Recession,
                             doParametric=FALSE,
                             doStructural=FALSE,
                             parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    RBaM:: BaM(mod=M_nonVAR,
               data=D,
               workspace = temp.folder.Recession,
               dir.exe = file.path(find.package("RBaM"), "bin"),
               pred=list(totalU,paramU,maxpost), # list of predictions
               remnant = remnant_prior,
               doCalib=FALSE,
               doPred=TRUE)

    # Get total uncertainty from the rating curve
    TotalUenv=read.table(file.path(temp.folder.Recession,'H_rec_TotalU.env'),header = TRUE)

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(TotalUenv))
    TotalUenv[cols_to_convert] <- suppressWarnings(lapply(TotalUenv[cols_to_convert],
                                                          function(x) as.numeric(x)))

    # Add parametric uncertainty
    ParametricUenv=read.table(file.path(temp.folder.Recession,'H_rec_ParamU.env'),header = TRUE)

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(ParametricUenv))
    ParametricUenv[cols_to_convert] <- suppressWarnings(lapply(ParametricUenv[cols_to_convert],
                                                               function(x) as.numeric(x)))

    # Add maxpost rating curve
    MAPREC=read.table(file.path(temp.folder.Recession,'H_rec_Maxpost.spag'))

    # handle when value is almost zero (detected as character)
    cols_to_convert=suppressWarnings(character_check(MAPREC))
    MAPREC[cols_to_convert] <- suppressWarnings(lapply(MAPREC[cols_to_convert],
                                                       function(x) as.numeric(x)))

    allData[[j]]=list(CalData=CalData,
                      TotalUenv=TotalUenv,
                      ParametricUenv=ParametricUenv,
                      MAPREC=MAPREC,
                      HgridPlot=t_grid)
    j=j+1
  }
  return(allData)
}

