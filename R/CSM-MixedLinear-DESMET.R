#' CENTRALIZED STATISTICAL MONITORING : DESMET METHOD
#'
#' @param data a multicenter clinical trial data containing one or several continuous variable
#'
#' @return table
#'
#' @import data.table
#' @import lme4
#' @import dplyr
#'
#'
#' @export
#'
#' @examples
#'



## Construction of the Desmet statistic test (ZSitei) - M?thode Desmet et al. (2014)

CONSTRUCT_ZDI <-function(data=data){

  ResultModel<-lme4::lmer(y.ij ~ 1 +(1|center),data=data)
  summary(ResultModel)

  ParVar <- data.frame(lme4::VarCorr(ResultModel))
  SigmaCcarre <- ParVar[1, c(4)]          # Center variance
  SigmaRcarre <- ParVar[2, c(4)]          # Residual variance
  MoyenneGlobale <- data.frame(lme4::fixef(ResultModel))[1,]

  Vecteur_centre <- unique(data$center)

  #  n.center <- length(unique(data$center))
  tab <- data.frame(t(c(1:7)))
  colnames(tab) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT","Detect_Zi")

  for (i in 1:length(Vecteur_centre)){
    Valeur_du_site <- Vecteur_centre[i]
    centre.uniq <- subset(data, center==Valeur_du_site)
    n.subject <- nrow(centre.uniq)

    temptab <- data.frame(t(c(1:7)))
    colnames(temptab) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT","Detect_Zi")

    temptab$Site <- Vecteur_centre[i]
    temptab$MoyenCentre <- mean(centre.uniq$y.ij, na.rm=TRUE)
    temptab$DiffMoyenne = round((temptab$MoyenCentre - MoyenneGlobale),4)
    temptab$Denominateur = round((SigmaCcarre + SigmaRcarre/n.subject),4)
    temptab$ZSitei = (temptab$MoyenCentre - MoyenneGlobale)/sqrt(SigmaCcarre + SigmaRcarre/n.subject)
    temptab$ProbT = round(2*pnorm(-abs(temptab$ZSitei)),4)
    temptab$Detect_Zi = ifelse(temptab$ProbT < 0.05, 1, 0)

    tab <- rbind(
      tab, temptab)

  }

  SortieResMoySite <- tab[-1,]
  return(SortieResMoySite)

}
