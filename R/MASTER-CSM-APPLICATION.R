#' CENTRALIZED STATISTICAL MONITORING : MASTER PROGRAM FOR CSM APPLICATION
#'
#' @param data is the table containing a multicenter clinical trial data on which the function program for each CSM methode will be performed.
#' @import data.table
#' @return table
#'
#' @import data.table
#' @import magrittr
#' @import stats
#' @import lme4
#' @import dplyr
#' @import sqldf
#' @import tcltk
#'
#' @export

## Construction du Programme Maître

APPLICATION_CSM_MOY_GLOBAL <-function(data=data){
  data <- data.table::data.table(data)
  EcartType_Tab <- data %>% dplyr::group_by(center) %>% dplyr::summarise(EcartypeSite=sd(y.ij_NCR, na.rm = TRUE),
                                                           MeanSite=mean(y.ij_NCR, na.rm = TRUE))
  data <- merge(data, EcartType_Tab, by=c("center"))
  data$y.ij <- (data$y.ij_NCR-mean(data$y.ij_NCR, na.rm=T))/data$EcartypeSite

  RESULT_REPLIC <- data.frame(t(c(1:21)))
  colnames(RESULT_REPLIC) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT","Detect_Zi",
                              "MoySite","MoyAutrSite","TestMoyenne","Detect_TT",
                              "MoyenneGlobaleD","VarianceGlobale","Distancei","Predict_i","DistDetect_i",
                              "Effectif", "MeanVal","Seuil_05","Seuil_95","Flag")

  RESULT_SIMUL_ZDI = CONSTRUCT_ZDI(data)
  RESULT_SIMUL_BAY = CONSTRUCT_BAYES(data)
  RESULT_SIMUL_TTEST = CONSTRUCT_TestMOY(data)
  RESULT_SIMUL_DIST = CONSTRUCT_DISTANCE(data)
  RESULT_SIMUL_ZTT <- merge(RESULT_SIMUL_ZDI, RESULT_SIMUL_TTEST, by="Site")
  RESULT_SIMUL_ZTT_DIST <- merge(RESULT_SIMUL_ZTT, RESULT_SIMUL_DIST, by="Site")
  RESULT_SIMUL <- merge(RESULT_SIMUL_ZTT_DIST, RESULT_SIMUL_BAY, by="Site")

  RESULT_REPLIC <- rbind(RESULT_REPLIC, RESULT_SIMUL)
  RESULT_REPLIC_SORTIE <- RESULT_REPLIC[-1,]

  data_restreint <- unique(data[,c("center","MeanSite","EcartypeSite")]) %>% dplyr::rename(Site=center)
  RESULT_REPLIC_SORTIE <- merge(data_restreint, RESULT_REPLIC_SORTIE, by=c("Site"))

  return(RESULT_REPLIC_SORTIE)
}


