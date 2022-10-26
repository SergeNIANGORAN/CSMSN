#' CENTRALIZED STATISTICAL MONITORING : MASTER PROGRAM FOR CSM SIMULATION
#'
#' @param SiteAtypiq is the maximum number of atypical site in the simulations, that is the number of atypical site will varying from 1 to SiteAtypiq. His default value \code{6}
#' @param MaxEcart . His default value \code{6}
#' @param NombreSimul . His default value \code{6}
#' @param n.center indicate the number of investigative sites in the multicenter trials. His default value \code{6}
#' @param n.subject is the number of subject in a site (considering that all centers have the same subjects number)
#' @param mu is the mean of the non atypical centers. His default value \code{10}
#' @param sigma2c . His default value \code{1}
#' @param sigma2r . His default value \code{4} according the default value of the mean
#' @import data.table
#' @return table
#'
#' @import data.table
#' @import lme4
#' @import dplyr
#' @import sqldf
#' @import tcltk
#'
#' @export
#'
#' @examples
#'

MASTER_CSM_MOY_GLOBAL_SIMS <-function(NSiteAtypiq=2, MaxEcart=5, NombreSimul=2,
                                     n.subject = 50, n.center = 10, sigma2c = 1,
                                     sigma2r = 4, mu = 0){

  SORTIE_REPLIC_Delta_Taux <- data.frame(t(c(1:22)))
  colnames(SORTIE_REPLIC_Delta_Taux) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT",
                                         "Effectif", "MeanVal","Seuil_05","Seuil_95","Flag",
                                         "MoySite","MoyAutrSite","TestMoyenne",
                                         "MoyenneGlobaleD","VarianceGlobale","Distancei","Predict_i","DistDetect_i",
                                         "NbSimulation","DeltaMoyenne","Taux")

  pb   <- txtProgressBar(1, NSiteAtypiq+1, style=3)

  for (Taux in 1:NSiteAtypiq){
    for (DeltaVarie in 1:MaxEcart){
      # message(glue::glue("Scenario execute : {Taux} sites atypiques, soit un taux de {round(Taux/n.center*100,1)} et un ecart a la moyenne de {DeltaVarie} "))

      RESULT_REPLIC <- data.frame(t(c(1:20)))
      colnames(RESULT_REPLIC) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT",
                                  "MoySite","MoyAutrSite","TestMoyenne",
                                  "MoyenneGlobaleD","VarianceGlobale","Distancei","Predict_i","DistDetect_i",
                                  "Effectif", "MeanVal","Seuil_05","Seuil_95","Flag","NbSimulation")

      for (NbreSim in 1:NombreSimul) {
        message(glue::glue("Il a ete realise {NbreSim} simulations sur {NombreSimul} "))

        NbreSimulation = NbreSim
        data <- data.generebis( n.subject = n.subject, n.center = n.center,
                                sigma2c = sigma2c,
                                sigma2r = sigma2r,
                                mu = mu,
                                centreAtypic = Taux,
                                delta = DeltaVarie)

        RESULT_SIMUL_ZDI = CONSTRUCT_ZDI(data)
        RESULT_SIMUL_BAY = CONSTRUCT_BAYES(data)
        RESULT_SIMUL_TTEST = CONSTRUCT_TestMOY(data)
        RESULT_SIMUL_DIST = CONSTRUCT_DISTANCE(data)
        RESULT_SIMUL_ZTT <- merge(RESULT_SIMUL_ZDI, RESULT_SIMUL_TTEST, by="Site")
        RESULT_SIMUL_ZTT_DIST <- merge(RESULT_SIMUL_ZTT, RESULT_SIMUL_DIST, by="Site")
        RESULT_SIMUL <- merge(RESULT_SIMUL_ZTT_DIST, RESULT_SIMUL_BAY, by="Site")
        RESULT_SIMUL$NbSimulation <- NbreSim

        RESULT_REPLIC <- rbind(RESULT_REPLIC, RESULT_SIMUL)
      }

      RESULT_REPLIC_SORTIE <- RESULT_REPLIC[-1,]
      RESULT_REPLIC_SORTIE$DeltaMoyenne <- DeltaVarie
      RESULT_REPLIC_SORTIE$Taux <- Taux

      SORTIE_REPLIC_Delta_Taux <- rbind(SORTIE_REPLIC_Delta_Taux, RESULT_REPLIC_SORTIE)
    }
  }

  SORTIE_REPLIC_Delta_Taux_Final <- SORTIE_REPLIC_Delta_Taux[-1,]
  ASSIGNATION_MOY <- SORTIE_REPLIC_Delta_Taux_Final

  ASSIGNATION_MOY$TauPercent <- round(ASSIGNATION_MOY$Taux/n.center*100, 0)
  ASSIGNATION_MOY$Type_Assign = ifelse(ASSIGNATION_MOY$Site >= (n.center - ASSIGNATION_MOY$Taux + 1), 1, 2)
  ASSIGNATION_MOY$Detect = ifelse(ASSIGNATION_MOY$ProbT < 0.05, 1, 0)

  Req1 <- "SELECT COUNT(Distinct Site) AS NbreSiteAssig, SUM(Detect) AS NB_ATYPIQUE, MAX(NbSimulation) AS NbreReplic, DeltaMoyenne, Taux, Type_Assign FROM ASSIGNATION_MOY GROUP BY DeltaMoyenne, Taux, Type_Assign;"

  ESSAIPERFO <- sqldf::sqldf(Req1)
  ESSAIPERFOTRI <- ESSAIPERFO[order(ESSAIPERFO$Taux, ESSAIPERFO$DeltaMoyenne, ESSAIPERFO$NbreSiteAssig, ESSAIPERFO$Type_Assign),]
  ESSAIPERFOTRI$Cpt <- cumsum(!duplicated(ESSAIPERFOTRI[,c("Taux","DeltaMoyenne")]))

  SENS <- subset(ESSAIPERFOTRI, Taux == NbreSiteAssig & Type_Assign==1)
  SPEC <- subset(ESSAIPERFOTRI, !(Taux == NbreSiteAssig & Type_Assign==1))

  SENS_RESULT <- SENS
  SENS_RESULT$FN <- SENS_RESULT$NbreSiteAssig*SENS_RESULT$NbreReplic - SENS_RESULT$NB_ATYPIQUE
  SENS_RESULT$TP <- SENS_RESULT$NB_ATYPIQUE
  SENS_RESULT$SENSIBILITE <- round(SENS_RESULT$NB_ATYPIQUE*100/(SENS_RESULT$NbreSiteAssig*SENS_RESULT$NbreReplic),1)
  SENS_RESULT$SCENARIO <- SENS_RESULT$Cpt
  SENS_RESULT <- SENS_RESULT[, c(3:5,8:11)]

  SPEC_RESULT <- SPEC
  SPEC_RESULT$TN <- SPEC_RESULT$NbreSiteAssig*SPEC_RESULT$NbreReplic - SPEC_RESULT$NB_ATYPIQUE
  SPEC_RESULT$FP <- SPEC_RESULT$NB_ATYPIQUE
  SPEC_RESULT$SPECIFICITE <- round(SPEC_RESULT$TN*100/(SPEC_RESULT$NbreSiteAssig*SPEC_RESULT$NbreReplic),1)
  SPEC_RESULT$SCENARIO <- SPEC_RESULT$Cpt
  SPEC_RESULT <- SPEC_RESULT[, c(8:11)]

  Req_Bay <- " SELECT COUNT(Distinct Site) AS NbreSiteAssig, SUM(Flag) AS NB_ATYPIQUE_BAY, MAX(NbSimulation) AS NbreReplic, DeltaMoyenne, Taux, Type_Assign FROM ASSIGNATION_MOY GROUP BY DeltaMoyenne, Taux, Type_Assign;"

  ESSAIPERFO_BAY <- sqldf::sqldf(Req_Bay)
  ESSAIPERFOTRI_BAY <- ESSAIPERFO_BAY[order(ESSAIPERFO_BAY$Taux, ESSAIPERFO_BAY$DeltaMoyenne, ESSAIPERFO_BAY$NbreSiteAssig, ESSAIPERFO_BAY$Type_Assign),]
  ESSAIPERFOTRI_BAY$Cpt <- cumsum(!duplicated(ESSAIPERFOTRI_BAY[,c("Taux","DeltaMoyenne")]))

  SENS_BAY <- subset(ESSAIPERFOTRI_BAY, Taux == NbreSiteAssig & Type_Assign==1)
  SPEC_BAY <- subset(ESSAIPERFOTRI_BAY, !(Taux == NbreSiteAssig & Type_Assign==1))

  SENS_RESULT_BAY <- SENS_BAY
  SENS_RESULT_BAY$FN_BAY <- SENS_RESULT_BAY$NbreSiteAssig*SENS_RESULT_BAY$NbreReplic - SENS_RESULT_BAY$NB_ATYPIQUE_BAY
  SENS_RESULT_BAY$TP_BAY <- SENS_RESULT_BAY$NB_ATYPIQUE_BAY
  SENS_RESULT_BAY$SENSIBILITE_BAY <- round(SENS_RESULT_BAY$NB_ATYPIQUE_BAY*100/(SENS_RESULT_BAY$NbreSiteAssig*SENS_RESULT_BAY$NbreReplic),1)
  SENS_RESULT_BAY$SCENARIO <- SENS_RESULT_BAY$Cpt
  SENS_RESULT_BAY <- SENS_RESULT_BAY[, c(3:5,8:11)]

  SPEC_RESULT_BAY <- SPEC_BAY
  SPEC_RESULT_BAY$TN_BAY <- SPEC_RESULT_BAY$NbreSiteAssig*SPEC_RESULT_BAY$NbreReplic - SPEC_RESULT_BAY$NB_ATYPIQUE_BAY
  SPEC_RESULT_BAY$FP_BAY <- SPEC_RESULT_BAY$NB_ATYPIQUE_BAY
  SPEC_RESULT_BAY$SPECIFICITE_BAY <- round(SPEC_RESULT_BAY$TN_BAY*100/(SPEC_RESULT_BAY$NbreSiteAssig*SPEC_RESULT_BAY$NbreReplic),1)
  SPEC_RESULT_BAY$SCENARIO <- SPEC_RESULT_BAY$Cpt
  SPEC_RESULT_BAY <- SPEC_RESULT_BAY[, c(8:11)]

  ASSIGNATION_MOY$Detect_TT = ifelse(ASSIGNATION_MOY$TestMoyenne < 0.05, 1, 0)
  Req2 <- "SELECT COUNT(Distinct Site) AS NbreSiteAssig, SUM(Detect_TT) AS NB_ATYPIQ_TT, MAX(NbSimulation) AS NbreReplic, DeltaMoyenne, Taux, Type_Assign FROM ASSIGNATION_MOY GROUP BY DeltaMoyenne, Taux, Type_Assign;"

  ESSAIPERFO_TT <- sqldf::sqldf(Req2)
  ESSAIPERFOTRI_TT <- ESSAIPERFO_TT[order(ESSAIPERFO_TT$Taux, ESSAIPERFO_TT$DeltaMoyenne, ESSAIPERFO_TT$NbreSiteAssig, ESSAIPERFO_TT$Type_Assign),]
  ESSAIPERFOTRI_TT$Cpt <- cumsum(!duplicated(ESSAIPERFOTRI_TT[,c("Taux","DeltaMoyenne")]))

  SENS_TT <- subset(ESSAIPERFOTRI_TT, Taux == NbreSiteAssig & Type_Assign==1)
  SPEC_TT <- subset(ESSAIPERFOTRI_TT, !(Taux == NbreSiteAssig & Type_Assign==1))

  SENS_RESULT_TT <- SENS_TT
  SENS_RESULT_TT$FN_TT <- SENS_RESULT_TT$NbreSiteAssig*SENS_RESULT_TT$NbreReplic - SENS_RESULT_TT$NB_ATYPIQ_TT
  SENS_RESULT_TT$TP_TT <- SENS_RESULT_TT$NB_ATYPIQ_TT
  SENS_RESULT_TT$SENSIBILITE_TT <- round(SENS_RESULT_TT$NB_ATYPIQ_TT*100/(SENS_RESULT_TT$NbreSiteAssig*SENS_RESULT_TT$NbreReplic),1)
  SENS_RESULT_TT$SCENARIO <- SENS_RESULT_TT$Cpt
  SENS_RESULT_TT <- SENS_RESULT_TT[, c(3:5,8:11)]

  SPEC_RESULT_TT <- SPEC_TT
  SPEC_RESULT_TT$TN_TT <- SPEC_RESULT_TT$NbreSiteAssig*SPEC_RESULT_TT$NbreReplic - SPEC_RESULT_TT$NB_ATYPIQ_TT
  SPEC_RESULT_TT$FP_TT <- SPEC_RESULT_TT$NB_ATYPIQ_TT
  SPEC_RESULT_TT$SPECIFICITE_TT <- round(SPEC_RESULT_TT$TN_TT*100/(SPEC_RESULT_TT$NbreSiteAssig*SPEC_RESULT_TT$NbreReplic),1)
  SPEC_RESULT_TT$SCENARIO <- SPEC_RESULT_TT$Cpt
  SPEC_RESULT_TT <- SPEC_RESULT_TT[, c(8:11)]

  Req3 <- "SELECT COUNT(Distinct Site) AS NbreSiteAssig, SUM(DistDetect_i) AS NB_ATYPIQ_DIST, MAX(NbSimulation) AS NbreReplic, DeltaMoyenne, Taux, Type_Assign FROM ASSIGNATION_MOY GROUP BY DeltaMoyenne, Taux, Type_Assign;"

  ESSAIPERFO_DIST <- sqldf::sqldf(Req3)
  ESSAIPERFOTRI_DIST <- ESSAIPERFO_DIST[order(ESSAIPERFO_DIST$Taux, ESSAIPERFO_DIST$DeltaMoyenne, ESSAIPERFO_DIST$NbreSiteAssig, ESSAIPERFO_DIST$Type_Assign),]
  ESSAIPERFOTRI_DIST$Cpt <- cumsum(!duplicated(ESSAIPERFOTRI_DIST[,c("Taux","DeltaMoyenne")]))

  SENS_DIST <- subset(ESSAIPERFOTRI_DIST, Taux == NbreSiteAssig & Type_Assign==1)
  SPEC_DIST <- subset(ESSAIPERFOTRI_DIST, !(Taux == NbreSiteAssig & Type_Assign==1))

  SENS_RESULT_DIST <- SENS_DIST
  SENS_RESULT_DIST$FN_DIST <- SENS_RESULT_DIST$NbreSiteAssig*SENS_RESULT_DIST$NbreReplic - SENS_RESULT_DIST$NB_ATYPIQ_DIST
  SENS_RESULT_DIST$TP_DIST <- SENS_RESULT_DIST$NB_ATYPIQ_DIST
  SENS_RESULT_DIST$SENSIBILITE_DIST <- round(SENS_RESULT_DIST$NB_ATYPIQ_DIST*100/(SENS_RESULT_DIST$NbreSiteAssig*SENS_RESULT_DIST$NbreReplic),1)
  SENS_RESULT_DIST$SCENARIO <- SENS_RESULT_DIST$Cpt
  SENS_RESULT_DIST <- SENS_RESULT_DIST[, c(3:5,8:11)]

  SPEC_RESULT_DIST <- SPEC_DIST
  SPEC_RESULT_DIST$TN_DIST <- SPEC_RESULT_DIST$NbreSiteAssig*SPEC_RESULT_DIST$NbreReplic - SPEC_RESULT_DIST$NB_ATYPIQ_DIST
  SPEC_RESULT_DIST$FP_DIST <- SPEC_RESULT_DIST$NB_ATYPIQ_DIST
  SPEC_RESULT_DIST$SPECIFICITE_DIST <- round(SPEC_RESULT_DIST$TN_DIST*100/(SPEC_RESULT_DIST$NbreSiteAssig*SPEC_RESULT_DIST$NbreReplic),1)
  SPEC_RESULT_DIST$SCENARIO <- SPEC_RESULT_DIST$Cpt
  SPEC_RESULT_DIST <- SPEC_RESULT_DIST[, c(8:11)]

  RESULT_PERFORMANCE_LM <- merge(SENS_RESULT, SPEC_RESULT, by ="SCENARIO")
  RESULT_PERFORMANCE_TT <- merge(SENS_RESULT_TT, SPEC_RESULT_TT, by ="SCENARIO")
  RESULT_PERFORMANCE_DIST <- merge(SENS_RESULT_DIST, SPEC_RESULT_DIST, by ="SCENARIO")
  RESULT_PERFORMANCE_BAY <- merge(SENS_RESULT_BAY, SPEC_RESULT_BAY, by ="SCENARIO")

  RESULT_PERFORMANCE_ZTT <- merge(RESULT_PERFORMANCE_LM, RESULT_PERFORMANCE_TT[,-c(2:4)], by ="SCENARIO")
  RESULT_PERFORMANCE_DISTBAY <- merge(RESULT_PERFORMANCE_DIST, RESULT_PERFORMANCE_BAY[,-c(2:4)],  by ="SCENARIO")

  RESULT_PERFORMANCE <- merge(RESULT_PERFORMANCE_ZTT, RESULT_PERFORMANCE_DISTBAY[,-c(2:4)], by ="SCENARIO")
  RESULT_PERFORMANCE$TauPercent <- round(RESULT_PERFORMANCE$Taux/n.center*100, 0)
  return(list(RESULT_PERFORMANCE,ASSIGNATION_MOY))
}

