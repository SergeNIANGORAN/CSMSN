#' CENTRALIZED STATISTICAL MONITORING : HATAYAMA METHOD
#' Function to implemente the Hatayama CSM method which uses the Bayesian Finite
#' Mixture Models (FMM) in order to detect an atypical center in multicenter trial.
#' This method like the others of this package implementes the CSM for the
#'  continuous variable data
#'
#' @param dataHatayama
#'
#' @return Table
#'
#' @import data.table
#' @import jagsUI
#' @import dplyr
#' @import magrittr
#'
#' @export


CONSTRUCT_BAYES <-function(dataHatayama){

  dataHatayama <- dataHatayama %>% dplyr::rename(Site = center, Val=y.ij)
  Val=dataHatayama$Val
  Site=dataHatayama$Site

  # Calcul des param?tres d'effectifs et de moyennes sites
  Vect <- data.table::data.table(table(dataHatayama$Site))
  Vect$V1 <- as.numeric(Vect$V1)

  Vect2 <- Vect %>% dplyr::rename(Site = V1, Effectif=N)
  resMean <- aggregate(Val ~ Site, data = dataHatayama, FUN = mean)
  resMean$Val <- round(resMean$Val,1)
  # resMean$Site <- as.character(resMean$Site)

  Fus_NMoy <- merge(Vect2, resMean, by="Site")
  Fus_NMoy <- Fus_NMoy %>% dplyr::rename(MeanVal = Val)


  # Pour un mod?le de m?lange gaussien

  model <- paste("
model {
  # Likelihood:
  for(i in 1:N) {
    y[i] ~ dnorm(mu[i], phi[i])
    mu[i] <- muOfClust[ clust[i] ]
    phi[i] <- phiOfClust[ clust[i] ]
    clust[i] ~ dcat( pClust[1:Nclust] )
    sigma2[i] <- 1/phi[i]
  }
  # Prior:

  for ( clustIdx in 1: Nclust ) {
    muOfClust[clustIdx] ~ dnorm(mu_0, phi_0)
    phiOfClust[clustIdx] ~ dgamma(alpha, beta)
  }
  pClust[1:Nclust] ~ ddirch(onesRepNclust )
}
")
  writeLines(model,"modelHatayamaMix.txt")


  # The data specification:
  # Generate random data from known parameter values:

  y=Val
  N <- length(dataHatayama$Site)
  Nclust=3
  clust = rep(NA,N)
  clust[which.min(y)]=1 # smallest value assigned to cluster 1
  clust[which.max(y)]=2
  clust[which.min(y)]=3 # highest value assigned to cluster 2
  DT_Hatay <- list("y"=Val,
                   "N"=N, "mu_0"=0, "phi_0"=0.01,
                   "alpha"=0.001, "beta"=0.001,
                   Nclust = Nclust ,
                   clust = clust ,
                   onesRepNclust = rep(1,Nclust)
  )

  out_jags_Mix = jagsUI::jags(data = DT_Hatay,
                              parameters.to.save = c("muOfClust", "phiOfClust","pClust"),
                              model.file = "modelHatayamaMix.txt",
                              n.chains = 3,
                              parallel = FALSE,
                              n.adapt = 1000,
                              n.iter = 3000,
                              n.burnin = 1000,
                              n.thin = 1,
                              DIC = T)


  out_jags_Mix

  Out_jags_MixMean <- data.frame(out_jags_Mix$mean$pClust)
  MaxProb <-max(Out_jags_MixMean)


  Groupe=ifelse(MaxProb == Out_jags_MixMean[1,], 1, ifelse(MaxProb == Out_jags_MixMean[2,], 2,3))

  out_jags_Mix.sims.list <- data.frame(out_jags_Mix$sims.list)

  out_jags_Mix.sims.list$mu_Body_Dist <- out_jags_Mix.sims.list$muOfClust.Groupe

  if(Groupe == 1){
    out_jags_Mix.sims.list$mu_Body_Dist <- out_jags_Mix.sims.list$muOfClust.1
    out_jags_Mix.sims.list$SD_Body_Dist <- sqrt(1/out_jags_Mix.sims.list$phiOfClust.1)
  }
  if(Groupe == 2){
    out_jags_Mix.sims.list$mu_Body_Dist <- out_jags_Mix.sims.list$muOfClust.2
    out_jags_Mix.sims.list$SD_Body_Dist <- sqrt(1/out_jags_Mix.sims.list$phiOfClust.3)
  }
  if(Groupe == 3){
    out_jags_Mix.sims.list$mu_Body_Dist <- out_jags_Mix.sims.list$muOfClust.3
    out_jags_Mix.sims.list$SD_Body_Dist <- sqrt(1/out_jags_Mix.sims.list$phiOfClust.3)
  }

  Body_Distr_Data <- out_jags_Mix.sims.list[,c("mu_Body_Dist", "SD_Body_Dist")]
  NB_SITE <- nrow(table(dataHatayama$Site))

  BD_List_Final <- list()
  for(k in 1:NB_SITE){
    if (Fus_NMoy$Site[k]==k){ #Cette condition à ce niveau pourrait sauter
      repn_Site <- Fus_NMoy$Effectif[k]
    }

    Temp <- list()
    for (t in 1:repn_Site){
      Body_Distr_Data <- out_jags_Mix.sims.list[,c("mu_Body_Dist", "SD_Body_Dist")]
      Body_Distr_Data$Y_pred <- rnorm(n = nrow(Body_Distr_Data), Body_Distr_Data[,c("mu_Body_Dist")], Body_Distr_Data[,c("SD_Body_Dist")])
      Body_Distr_Data$rand <- runif(n = nrow(Body_Distr_Data))
      Temp[[t]] <- Body_Distr_Data
    }
    Temp2 <- do.call(rbind, Temp)
    Temp2 <- Temp2[order(Temp2$rand),]

    Temp2$repli<-rep(1:repn_Site,nrow(Temp2)/repn_Site)
    Temp2$repli_n<-rep(1:(nrow(Temp2)/repn_Site), each=repn_Site)

    MoyTmp <-data.frame(tapply(Temp2$Y_pred, Temp2$repli_n, mean))
    MoyTmp$repn <- 1:nrow(MoyTmp)
    MoyTmp$TailleSite <- max(Temp2$repli)
    MoyTmp$Site <- k
    MoyTmp <- MoyTmp %>% dplyr::rename(Y_predM = tapply.Temp2.Y_pred..Temp2.repli_n..mean.)
    BD_List_Final[[k]] <- MoyTmp
  }
  BD_List_Final_Temp <- do.call(rbind, BD_List_Final)


  Percent_CSM_Final <- list()
  for(k in 1:NB_SITE){
    PERC_BD_List_Final <- BD_List_Final_Temp[BD_List_Final_Temp$Site==k,]
    Site=k
    Seuil_05 <- quantile(PERC_BD_List_Final$Y_predM, probs = 0.05, na.rm=TRUE)
    Seuil_95 <- quantile(PERC_BD_List_Final$Y_predM, probs = 0.95, na.rm=TRUE)
    Percent_CSM <- cbind(Site,Seuil_05,Seuil_95)
    Percent_CSM_Final[[k]] <- Percent_CSM
  }

  Percent_CSM_Final1 <- do.call(rbind, Percent_CSM_Final)
  Percent_CSM_Final <- data.frame(Percent_CSM_Final1)
  Resultat_Bayes_CSM <- merge(Fus_NMoy,Percent_CSM_Final,by="Site")

  Resultat_Bayes_CSM$Flag = ifelse(Resultat_Bayes_CSM$MeanVal < Resultat_Bayes_CSM$Seuil_05 |
                                     Resultat_Bayes_CSM$MeanVal > Resultat_Bayes_CSM$Seuil_95,1,0)

  return(Resultat_Bayes_CSM)
}
