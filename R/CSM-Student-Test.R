#' CENTRALIZED STATISTICAL MONITORING : STUDENT TEST METHOD
#'
#' @param data a multicenter clinical trial data containing one or several continuous variable
#'
#' @return table
#'
#' @importFrom data.table data.table
#' @import dplyr
#'
#'
#' @export
## Student t-test CSM method for mean comparisons

CONSTRUCT_TestMOY <-function(data=data){

  Vecteur_centre <- unique(data$center)
  tab <- data.frame(t(c(1:5)))
  colnames(tab) <- c("Site","MoySite","MoyAutrSite","TestMoyenne","Detect_TT")
  for (i in 1:length(Vecteur_centre)){
    tab.data <- data
    centre.uniq <- subset(tab.data, center== Vecteur_centre[i])
    centre.uniq.Diff <- subset(tab.data, !(center==Vecteur_centre[i]))

    tab.data$SiteRef <-ifelse(tab.data$center == Vecteur_centre[i], 1, 0)
    tempTestVar <- var.test(y.ij ~ SiteRef, data = tab.data)$p.value

    if(tempTestVar >= 0.05){
      tempTestMoyenne <- t.test(y.ij ~ SiteRef, data = tab.data, var.equal=T)$p.value
    }
    if(tempTestVar < 0.05){
      tempTestMoyenne <- t.test(y.ij ~ SiteRef, data = tab.data, var.equal=F)$p.value
    }

    temptab <- data.frame(t(c(1:5)))
    colnames(temptab) <-c("Site","MoySite","MoyAutrSite","TestMoyenne","Detect_TT")

    temptab$Site <- Vecteur_centre[i]
    temptab$MoySite <- mean(centre.uniq$y.ij, na.rm=TRUE)
    temptab$MoyAutrSite <- mean(centre.uniq.Diff$y.ij, na.rm=TRUE)
    temptab$TestMoyenne <- tempTestMoyenne
    temptab$Detect_TT = ifelse(temptab$TestMoyenne < 0.05, 1, 0)
    tab.data <- tab.data[,c(1:4)]
    tab <- rbind(
      tab, temptab)
  }

  SortieResMoyTestSite <- tab[-1,]
  return(SortieResMoyTestSite)

}

