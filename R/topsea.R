#' @title topSEA
#'
#' @description returns a permutation  of SEA-chart which rearranges
#' the feature-sets according to the selected argument into ascending or
#' descending order.
#'
#' @param object A SEA-chart object which is the output of \code{SEA} function
#'
#' @param by Variable name by which the ordering should happen. It should be a column of SEA-chart.
#' The default is TDP_bound.
#'
#' @param thresh A real number between 0 and 1. If specified the values of the variable defined in \code{by}
#' will be threshold accordingly.
#'
#' @param n Integer. Number of raws of the output chart
#'
#' @param descending Logical. If \code{TRUE} The output chart is organized in a descending order
#'
#' @param cover An optional threshold for coverage, which must be a real number between 0 and 1.
#' If specified, feature-sets with a coverage lower than or equal to this value are removed.
#'
#'
#'
#' @return Returns a subset of SEA_chart sorted according to the arguments
#'
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso \code{\link{SEA}}
#'
#' @references
#' Mitra Ebrahimpoor, Pietro Spitali, Kristina Hettne, Roula Tsonaka, Jelle Goeman,
#' Simultaneous Enrichment Analysis of all Possible Gene-sets: Unifying Self-Contained
#' and Competitive Methods, Briefings in Bioinformatics,bbz074
#'
#' @examples
#' #See the examples for \code{\link{SEA}}
#'
#' @export
#'
#' @importFrom hommel hommel tdp localtest
#'
#'
topSEA=topsea<- function(object, by, thresh=NULL, descending=TRUE, n=20, cover){

  #save the call to function
  cl <- match.call()

  #evaluate by argument
  if(missing(by)) byName <- "TDP.estimate"
  else byName <- deparse(cl$by)

  #check the arguments
    if(sum(colnames(object) %in% c("ID","Name","Size","Coverage","TDP.bound","TDP.estimate",
                                   "SC.adjP","Comp.adjP"))< 5)

       stop('Maybe the SEA-chart object is not specified correctly.')

    if(! byName %in% colnames(object))
          stop('The argument by should match the colnames of SEA-chart!')

    if( byName %in% c("ID","Name"))
          stop('The chart can not be reordered by ID or Name!')

  #remove low cover
      if(!missing(cover)){
        object<-with(object, object[round(Coverage,4) >= cover,])
          if(nrow(object)==0)
             stop('Nothing selected, modify Cover value!')}

  #remove thoes below threshold
      if(!missing(thresh)){
        object <- with(object, object[object[,byName]>=thresh,])
          if(nrow(object)==0)
             stop('Nothing selected, modify threshold!')}

      if(descending==TRUE)
          SEA_top<-object[order(object[,byName], decreasing = TRUE),]

      if(descending==FALSE)
          SEA_top<-object[order(object[,byName], decreasing = FALSE),]

  #select the top
      if(nrow(object)>n)
        SEA_top<-SEA_top[1:n,]


  return(SEA_top)
}
