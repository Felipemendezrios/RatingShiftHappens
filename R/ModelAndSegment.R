recurvise.ModelAndSegmentation <- function(obs,
                                           time=1:length(obs),
                                           u=0*obs,
                                           nSmax=2,
                                           nMin= 1,
                                           nCycles=100,
                                           burn=0.5,
                                           nSlim=max(nCycles/10,1),
                                           temp.folder=file.path(tempdir(),'BaM')){
  # Initialization
  allRes=list() # store segmentation results for all nodes in a sequential list
  k=0 # Main counter used to control indices in allRes
  tree=data.frame() # store tree structure (parents - children relationship)
  p=1 # Auxiliary counter needed to keep track of children / parents indices
  level=0 # Recursion level. The tree is created level-by-level rather than branch-by-branch
  X=list(obs) # List of all nodes (each corresponding to a subseries of x) to be segmented at this level. Start with a unique node corresponding to the whole series
  TIME=list(time) # List of corresponding times
  U=list(u) # List of corresponding uncertainties
  indices=c(1) # Vector containing the indices of each node - same size as X
  parents=c(0) # Vector containing the indices of the parents of each node - same size as X
  continue=TRUE # Logical determining whether recursion should continue

  while(continue){
    level=level+1 # Increment recursion level
    nX=length(X) # Number of nodes at this level
    keepgoing=rep(NA,nX) # Should recursion continue for each node?
    newX=newTIME=newU=newIndices=newParents=c() # Will be used to update subseries, indices and parents at the end of each recursion level
    m=0 # Local counter used to control indices in the 4 vectors above => reset to 0 at each new level of the recursion
    for(j in 1:nX){ # Loop on each node
      k=k+1 # Increment main counter
      partial.segmentation=segmentation(obs=X[[j]],time=TIME[[j]],u=U[[j]],
                                        nSmax,nMin,nCycles,burn,nSlim,temp.folder) # Apply segmentation to subseries stored in node X[[j]]
      # Save results for this node
      allRes[[k]]=partial.segmentation
      # Save optimal number of segments
      nSopt=partial.segmentation$nS
      # Update recursion tree
      tree=rbind(tree,data.frame(indx=k,level=level,parent=parents[j],nS=nSopt))
      # This was the trickiest part: keeping track of indices and parents
      keepgoing[j]=nSopt>1 # if nS=1, segmentation will not continue for this node which is hence terminal
      if(keepgoing[j]){ # Save results for segmentation at next level
        for(i in 1:nSopt){ # Loop on each segment detected for the current node
          p=p+1 # Increment auxiliary counter
          m=m+1 # Increment local counter
          newX[[m]]=partial.segmentation$results[[nSopt]]$data.p$obs.p[[i]] # Save ith segment (on a total of nS)
          newTIME[[m]]=partial.segmentation$results[[nSopt]]$data.p$time.p[[i]] # Save corresponding times
          newU[[m]]=partial.segmentation$results[[nSopt]]$data.p$u.p[[i]] # Save corresponding uncertainty
          newParents[m]=indices[j] # At next level, the parent of this segment will be the index of current node
          newIndices[m]=p # At next level, the index of this segment will be p
        }
      }
    }
    # Check if recursion should continue at all, i.e. if at least one node is not terminal
    if(all(keepgoing==FALSE)) continue=FALSE
    # Update list of nodes to be further segmented at next level + parents and indices
    X=newX
    TIME=newTIME
    U=newU
    parents=newParents
    indices=newIndices
  }
  # Get terminal nodes
  terminal=which(tree$nS==1)

  # Get stable periods by adding information as period, segment,
  data <- c()
  for(i in 1:length(terminal)){
    data.stable.p=allRes[[terminal[i]]]$results[[1]]   #Save data from stable period

    node = data.frame(time=data.stable.p$data.p$time.p,
                      obs=data.stable.p$data.p$obs.p,
                      u=data.stable.p$data.p$u.p,
                      I95_lower=data.stable.p$data.p$obs.p+stats::qnorm(0.025)*data.stable.p$data.p$u.p,
                      I95_upper=data.stable.p$data.p$obs.p+stats::qnorm(0.975)*data.stable.p$data.p$u.p,
                      id = rep(i,length(data.stable.p$data.p$obs.p)))

    data = rbind(data,node)
  }

  data=data[order(data$time),]
  data$period <- rep(NA,nrow(data))

  period.counter=1
  for(i in 1:length(terminal)){
    data$period[which(unique(data$id)[i]==data$id)] <- period.counter
    period.counter=period.counter+1
  }

  data = data[,-which('id'==colnames(data))]

  # Get node with time shifts
  nodes.shift.time=which(tree$nS!=1)

  shift <- c()
  for(i in 1:length(nodes.shift.time)){
    nSopt.p = allRes[[nodes.shift.time[[i]]]]$nS
    results.p = allRes[[nodes.shift.time[[i]]]]$results

    shift.time.p=cbind(c(results.p[[nSopt.p]]$tau))

    for(j in 1:(nSopt.p-1)){

      shift.time.p.unc=data.frame(tau=as.numeric(shift.time.p[j,]),
                                  I95_lower=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                            probs=c(0.025)),
                                  I95_upper=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                            probs=c(0.975)))
      shift <- rbind(shift,
                     shift.time.p.unc)
    }
  }

  shift <- shift[order(shift$tau),]

  return(list(summary=list(data=data,
                           shift=shift),
              res=allRes,tree=tree))
}
