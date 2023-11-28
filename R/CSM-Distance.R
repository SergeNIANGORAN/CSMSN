#' CENTRALIZED STATISTICAL MONITORING : DISTANCE METHOD
#'
#' @param data a multicenter clinical trial data containing one or several continuous variable
#'
#' @return table
#'
#' @importFrom data.table data.table
#' @import dplyr
#'
#' @export

### Développement du test de calcul de distance

CONSTRUCT_DISTANCE <-function(data=data){

  MoyenneGlobaleD <- mean(data$y.ij)        # Moyenne Globale
  VarianceGlobale <- var(data$y.ij)        # Variance Globale

  tab <- data.frame(t(c(1:6)))
  colnames(tab) <-c("Site","MoyenneGlobaleD","VarianceGlobale","Distancei","Predict_i","DistDetect_i")
  n.center <- length(unique(data$center))

  for (i in 1:n.center){
    centre.uniq <- subset(data, center==i)
    centre.uniq$EcartCarre.ij <- (centre.uniq$y.ij - MoyenneGlobaleD)^2

    temptab <- data.frame(t(c(1:6)))
    colnames(temptab) <-c("Site","MoyenneGlobaleD","VarianceGlobale","Distancei","Predict_i","DistDetect_i")

    temptab$Site <- i
    temptab$Distancei <- sum(centre.uniq$EcartCarre.ij)/VarianceGlobale
    temptab$Predict_i <- temptab$Distancei/(nrow(centre.uniq) - 1)

    temptab$DistDetect_i = ifelse(temptab$Predict_i > qf(0.95, df1=(nrow(centre.uniq) - 1), df2=(nrow(data) - 1), lower.tail=TRUE),1,0)
    temptab$MoyenneGlobaleD = MoyenneGlobaleD
    temptab$VarianceGlobale = VarianceGlobale
    tab <- rbind(
      tab, temptab)
  }

  SortieResDistSite <- tab[-1,]
  return(SortieResDistSite)
}
