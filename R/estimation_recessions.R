#' Estimate two exponential recession model with five parameters
#'
#' @param data.object object, model object output by the RBaM function
#' @param equation_rec_model character, equation selected to modelling recessions
#' @param nCyclesrec integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burnrec real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlimrec integer, MCMC slim step
#' @param temp.folder.Recession directory, temporary directory to write computations
#'
#' @return list of following components:
#'  \enumerate{
#'    \item residuals: data frame, description below:
#'       \itemize{
#'          \item X1_obs: real value, duration of recession observed
#'          \item X1_true: real value, unbiased recession duration
#'          \item Y1_obs: real value, stage-recession observed
#'          \item Y1_unbiased: real value, unbiased stage-recession observed
#'          \item Y1_sim: real value, stage-recession simulated
#'          \item Y1_res: real value, residual between observed value and simulated value stage-recession
#'          \item Y1_stdres: real value, residual standar deviation at each stage-recession
#'      }
#'    \item summary_MCMC: data frame, statistical summary of all estimated parameters
#' }
#' @export
#' @importFrom RBaM parameter parameter_VAR xtraModelInfo mcmcOptions mcmcCooking remnantErrorModel BaM
#' @importFrom utils read.table
Estimate_Rec_M3 <- function(data.object,
                                    equation_rec_model,
                                    temp.folder.Recession,
                                    nCyclesrec=100,
                                    burnrec=0.5,
                                    nSlimrec=max(nCyclesrec/10,1)){

  rec.data=data.object$data
  Ncurves=max(rec.data$indx)

  # Define prior information :
  # Parameters
  alpha_1_k=RBaM::parameter_VAR(name='alpha_1_k',
                                index='indx',
                                d=data.object,
                                # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                                init=rep(max(rec.data$hrec)-min(rec.data$hrec),Ncurves), # first guesses
                                prior.dist=rep('FlatPrior+',Ncurves),
                                prior.par =rep(list(NULL), Ncurves)) # prior distributions

  alpha_2_k=RBaM::parameter_VAR(name='alpha_2_k',
                                index='indx',
                                d=data.object,
                                # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                                init=rep(max(rec.data$hrec)-min(rec.data$hrec),Ncurves), # first guesses
                                prior.dist=rep('FlatPrior+',Ncurves),
                                prior.par =rep(list(NULL), Ncurves)) # prior distributions

  # Details : https://en.wikipedia.org/wiki/Exponential_decay
  lambda_1=RBaM::parameter(name='lambda_1',
                           init=log(2)/stats::median(rec.data$time_rec),
                           prior.dist='FlatPrior+')

  lambda_2=RBaM::parameter(name='lambda_2',
                           init=log(2)/stats::median(rec.data$time_rec),
                           prior.dist='FlatPrior+')

  beta_k=RBaM::parameter_VAR(name='beta_k',
                             index='indx',
                             d=data.object,
                             # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                             init=rep(min(rec.data$hrec),Ncurves), # first guesses
                             prior.dist=rep('FlatPrior',Ncurves),
                             prior.par =rep(list(NULL), Ncurves)) # prior distributions

  # use xtraModelInfo to pass the names of the inputs and the formulas
  xtra=RBaM::xtraModelInfo(object=list(inputs=c('t'),formulas=equation_rec_model))
  # model
  mod=RBaM::model(ID='TextFile',nX=1,nY=1,par=list(alpha_1_k,lambda_1,alpha_2_k,lambda_2,beta_k),xtra=xtra)

  # Cooking
  mcmc_temp=RBaM::mcmcOptions(nCycles=nCyclesrec)
  cook_temp=RBaM::mcmcCooking(burn=burnrec,
                              nSlim=nSlimrec)
  # Error model
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Linear",
                                                par = list(RBaM::parameter(name="gamma1",
                                                                           init=1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))

  RBaM::BaM(mod=mod,
            data=data.object,
            workspace=temp.folder.Recession,
            mcmc=mcmc_temp,
            cook = cook_temp,
            dir.exe = file.path(find.package("RBaM"), "bin"),
            remnant = remnant_prior)

  residuals   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Summary.txt"),header=TRUE)

  return(list(residuals=residuals,
              summary_MCMC=summary.MCMC))
}

#' Estimate exponential recession model with four parameters
#'
#' @param data.object object, model object output by the RBaM function
#' @param equation_rec_model character, equation selected to modelling recessions
#' @param nCyclesrec integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burnrec real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlimrec integer, MCMC slim step
#' @param temp.folder.Recession directory, temporary directory to write computations
#'
#' @return list of following components:
#'  \enumerate{
#'    \item residuals: data frame, description below:
#'       \itemize{
#'          \item X1_obs: real value, duration of recession observed
#'          \item X1_true: real value, unbiased recession duration
#'          \item Y1_obs: real value, stage-recession observed
#'          \item Y1_unbiased: real value, unbiased stage-recession observed
#'          \item Y1_sim: real value, stage-recession simulated
#'          \item Y1_res: real value, residual between observed value and simulated value stage-recession
#'          \item Y1_stdres: real value, residual standar deviation at each stage-recession
#'      }
#'    \item summary_MCMC: data frame, statistical summary of all estimated parameters
#' }
#' @export
#' @importFrom RBaM parameter parameter_VAR xtraModelInfo mcmcOptions mcmcCooking remnantErrorModel BaM
#' @importFrom utils read.table
Estimate_Rec_BR1 <- function(data.object,
                                     equation_rec_model,
                                     temp.folder.Recession,
                                     nCyclesrec=100,
                                     burnrec=0.5,
                                     nSlimrec=max(nCyclesrec/10,1)){

  rec.data=data.object$data

  Ncurves=max(rec.data$indx)
  # Define prior information :
  # Parameters
  alpha_k=RBaM::parameter_VAR(name='alpha_k',
                              index='indx',
                              d=data.object,
                              # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                              init=rep(max(rec.data$hrec)-min(rec.data$hrec),Ncurves),
                              prior.dist=rep('FlatPrior+',Ncurves),
                              prior.par =rep(list(NULL), Ncurves))

  # Details : https://en.wikipedia.org/wiki/Exponential_decay
  lambda=RBaM::parameter(name='lambda',
                         init=log(2)/stats::median(rec.data$time_rec),
                         prior.dist='FlatPrior+')
  c=RBaM::parameter(name='c',
                    init=1,
                    prior.dist='FlatPrior+')

  beta_k=RBaM::parameter_VAR(name='beta_k',
                             index='indx',
                             d=data.object,
                             # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                             init=rep(min(rec.data$hrec),Ncurves), # first guesses
                             prior.dist=rep('FlatPrior',Ncurves),
                             prior.par =rep(list(NULL), Ncurves)) # prior distributions

  # use xtraModelInfo to pass the names of the inputs and the formulas
  xtra=RBaM::xtraModelInfo(object=list(inputs=c('t'),formulas=equation_rec_model))
  # model
  mod=RBaM::model(ID='TextFile',nX=1,nY=1,par=list(alpha_k,lambda,c,beta_k),xtra=xtra)

  # Cooking
  mcmc_temp=RBaM::mcmcOptions(nCycles=nCyclesrec)
  cook_temp=RBaM::mcmcCooking(burn=burnrec,
                              nSlim=nSlimrec)
  # Error model
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Linear",
                                                par = list(RBaM::parameter(name="gamma1",
                                                                           init=1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))

  RBaM::BaM(mod=mod,
            data=data.object,
            workspace=temp.folder.Recession,
            mcmc=mcmc_temp,
            cook = cook_temp,
            dir.exe = file.path(find.package("RBaM"), "bin"),
            remnant = remnant_prior)

  residuals   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Summary.txt"),header=TRUE)

  return(list(residuals=residuals,
              summary_MCMC=summary.MCMC))

}

