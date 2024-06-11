#' @title topSEA
#'
#' @description returns a plotof SEA-chart which illustrates
#' proportion of discoveries per pathway.
#'
#' @param object A SEA-chart object which is the output of \code{SEA} function
#'
#' @param by the Variable which will we mapped.
#' It should be either the TDP estimate or TDP bound.The default is TDP bound.
#'
#' @param thresh A real number between 0 and 1. Which will be used as
#' a visual aid to distinguish significant pathways
#'
#' @param n Integer. Number of rows from SEA-chart object to be plotted.
#'
#'
#' @return Returns a plot of SEA_chart according to the selected arguments
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
#' @importFrom ggplot2 ggplot
#'
#'
plotSEA=plotsea<- function(object, by="TDP.estimate", threshold=0.005, n=20){

#save the call to function
   cl <- match.call()

#check the arguments
 if(sum(colnames(object) %in% c("ID","Name","Size","Coverage","TDP.bound","TDP.estimate",
                                  "SC.adjP","Comp.adjP"))< 5)

    stop('Maybe the SEA-chart object is not specified correctly.')

 if( sum(by %in% c("TDP.bound", "TDP.estimate"))<0)
    stop('By argument is not specified correctly.')


#make the plot
object<-object[1:n,]

val<-max(object[,by], na.rm = T)
if(threshold>val)
  stop(paste("Threshold is higher than",round(val,5),
             "which is th maximum value of", by))
if(nrow(object)<n)
  stop(paste("Data has",nrow(object),"rows","can not select",n))

val<-val+0.005

seabar<-ggplot(object, aes(x = object[,"Name"],
                           y = as.numeric(object[,by]))) +
     geom_bar(stat = "identity",color="grey80") +
     geom_text(aes(label=object[,"Size"], y = threshold+0.01), color="black") +
     coord_flip() +
     scale_y_continuous(limits = c(0, val),breaks = seq(0, val, 0.01),
                     labels = abs(seq(0, val, 0.01))) +
     theme_bw()+
     geom_hline(yintercept=threshold, linetype="dashed",
                color = "black", size=0.5)+
  xlab("Pathway Name \n")+
  ylab(paste(by))

  return(seabar)

}
