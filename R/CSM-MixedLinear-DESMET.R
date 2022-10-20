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
  SigmaCcarre <- ParVar[1, c(4)]          # Center Variance
  SigmaRcarre <- ParVar[2, c(4)]          # Residual Variance
  MoyenneGlobale <- data.frame(lme4::fixef(ResultModel))[1,]

  n.center <- length(unique(data$center))
  tab <- data.frame(t(c(1:6)))
  colnames(tab) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT")

  for (i in 1:n.center){
    centre.uniq <- subset(data, center==i)
    n.subject <- nrow(centre.uniq)

    temptab <- data.frame(t(c(1:6)))
    colnames(temptab) <-c("Site","MoyenCentre","DiffMoyenne","Denominateur","ZSitei","ProbT")

    temptab$Site <- i
    temptab$MoyenCentre <- mean(centre.uniq$y.ij)
    temptab$DiffMoyenne = temptab$MoyenCentre - MoyenneGlobale
    temptab$Denominateur = SigmaCcarre + SigmaRcarre/n.subject
    temptab$ZSitei = (temptab$MoyenCentre - MoyenneGlobale)/sqrt(SigmaCcarre + SigmaRcarre/n.subject)
    temptab$ProbT = 2*pnorm(-abs(temptab$ZSitei))

    tab <- rbind(
      tab, temptab)

  }

  SortieResMoySite <- tab[-1,]
  return(SortieResMoySite)

}
