#' @title SEA
#'
#' @description returns SEA chart (a data.frame) including the test results and estimates for the specified
#' feature-sets from \code{pathlist}.
#'
#' @param pvalue Vector of p-values. It can be the name of the covariate representing the Vector of
#' all raw p-values in the \code{data} or a single vector but in the latter case it should match the
#' \code{featureIDs} vector
#'
#' @param featureIDs Vector of feature IDs. It can be the name of the covariate representing the IDs in the
#' \code{data} or a single vector but in the latter case it should match the \code{pvalue} vector
#'
#' @param data Optional data frame or matrix containing the variables in \code{pvalue} and \code{featureIDs}
#'
#' @param pathlist A list containing pathways defined by \code{featureIDs}. Checkout the vignette
#' for more details and available codes to create your own pathway
#'
#' @param select A vector. Number or names of pathways of interest from the \code{pathlist} of choice.
#' If missing, all pathways of the database will be included
#'
#' @param tdphat Logical. If \code{TRUE} the point estimate of the True Discoveries Proportion
#' within each pathway will be calculated
#'
#' @param selfcontained Logical. If \code{TRUE} the self-contained null hypothesis will be tested
#' for each pathway and the corresponding adj. p-value is returned
#'
#' @param competitive Logical. If \code{TRUE} the default competitive null hypothesis will be tested
#' for each pathway and the corresponding adj. p-value is returned, you can define a threshold with
#' \code{thresh} argument
#'
#' @param thresh A real number between 0 and 1. If specified, the competitive null hypothesis will be tested
#' against this threshold for each pathway and the corresponding adj. p-value is returned
#'
#' @param alpha The type I error allowed for TDP bound. The default is 0.05.
#'
#' @return A data.frame is returned including a list of pathways with corresponding TDP bound estimate,
#' and if specified, TDP point estimate and adjusted p-values
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso \code{\link{setTest}}, \code{\link{topSEA}},
#'
#' @references
#' Mitra Ebrahimpoor, Pietro Spitali, Kristina Hettne, Roula Tsonaka, Jelle Goeman,
#' Simultaneous Enrichment Analysis of all Possible Gene-sets: Unifying Self-Contained
#' and Competitive Methods, Briefings in Bioinformatics, , bbz074, https://doi.org/10.1093/bib/bbz074
#'
#' @examples
#'
#' \dontrun{
#' ##Generate a vector of pvalues
#' set.seed(159)
#'
#' m<- 100
#' pvalues <- runif(m,0,1)^5
#' featureIDs <- as.character(1:m)
#'
#' # perform a self-contained test for all features
#' settest(pvalues, featureIDs, testype = "selfcontained")
#'
#' # create 3 random pathway of size 60, 20 and 45
#' randpathlist=list(A=as.character(c(sample(1:m, 60))),
#'              B=as.character(c(sample(1:m, 20))),
#'              C=as.character(c(sample(1:m, 45))))
#'
#'
#' # get the seachart for the whole pathlist
#' set(pvalues, featureIDs, pathlist=randpathlist)
#'
#' # get the seachart for only first two pathways of the randpathlist
#' sea(pvalues, featureIDs, pathlist=randpathlist, select=1:2)
#' }
#' @export
#'
#' @importFrom hommel hommel
#'
#'

SEA=sea<- function(pvalue, featureIDs, data, pathlist, select,
               tdphat=TRUE, selfcontained=TRUE, competitive=TRUE, thresh=NULL, alpha= 0.05){

  #save the call to function
  cl<-match.call()

  #check data
  if(missing(data)){
    #evaluate pval and featureIDs
    if(length(pvalue)!=length(featureIDs)) stop('The arguments pvalue and featureIDs should match!')
    pval<-pvalue
    geneid<-featureIDs

  }else{

    if(is.matrix(data)) data <- data.frame(data)

    #evaluate pval and featureIDs, which may be one of the colnames of data
    pval <- eval(cl$pvalue, data, parent.frame())
    geneid <- eval(cl$featureIDs, data, parent.frame()) }

  if(any(pval<0)) stop('Negetive values passed as pvalue!')

  #the pathlist
   if (missing(pathlist))
    stop('Pathlist must be specified!')

  #creat the SEA_chart matching the number of pathways
    np<-length(pathlist)

    SEA_chart <- data.frame(ID=1:np,
                            Name=character(length = np),
                            Size=as.numeric(rep(NA,np)),
                            Coverage=as.numeric(rep(NA,np)),
                            TDP.bound=as.numeric(rep(NA,np)),
                            TDP.estimate=as.numeric(rep(NA,np)),
                            SC.adjP=as.numeric(rep(NA,np)),
                            Comp.adjP=as.numeric(rep(NA,np)),
                            Compc.adjP=as.numeric(rep(NA,np)), stringsAsFactors=FALSE)

    if(!missing(select)) {

      paths<-pathlist[select]
      SEA_chart<-SEA_chart[select,]
        } else paths<-pathlist

    if(sum(geneid %in% unlist(paths))==0) stop(paste("Pathway should match the data!"))

  #localtest functions per testtype

    obj<-hommel::hommel(pval)

    tdpfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                tdp(obj, ix =which(geneid  %in% x )),NA)}

    tdp50func<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                  tdp(obj, ix =which(geneid  %in% x ), alpha=0.5),NA)}

    selfadjPfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                     localtest(obj, ix =which(geneid  %in% x ),
                                            tdp=0),NA)}

    CompadjPfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                     localtest(obj, ix =which(geneid  %in% x ),
                                           tdp=tdp(obj)),NA)}

    thrCadjPfunc<-function(x){ifelse(length(which(geneid  %in% x ))>0,
                                     localtest(obj, ix =which(geneid  %in% x ),
                                             tdp=thresh),NA)}

  #fill in the relevent values

  SEA_chart$Name<-names(paths)

  SEA_chart$Size<-sapply(paths,FUN=function(x){length(x)},simplify = TRUE)

  SEA_chart$Coverage<-sapply(paths,FUN=function(x){round(length(which(geneid  %in% x ))/length(x),2)},
                            simplify = TRUE)


    SEA_chart$"TDP.bound"<-sapply(paths,tdpfunc,simplify = TRUE)

    if(tdphat==TRUE)
      SEA_chart$"TDP.estimate"<-sapply(paths,tdp50func,simplify = TRUE)
    else
      SEA_chart$"TDP.estimate"<-NA

    if(selfcontained==TRUE)
      SEA_chart$"SC.adjP"<-sapply(paths,selfadjPfunc,simplify = TRUE)

    else
      SEA_chart$"SC.adjP"<-NA

    if(competitive==TRUE)
      SEA_chart[,"Comp.adjP"]<-sapply(paths,CompadjPfunc,simplify = TRUE)
    else
      SEA_chart[,"Comp.adjP"]<-NA

    if(!missing(thresh)){
      SEA_chart[,"Compc.adjP"]<-sapply(paths,thrCadjPfunc,simplify = TRUE)
      colnames(SEA_chart)[colnames(SEA_chart)=="Compc.adjP"]<- paste0("Comp.",thresh,".adjP")}

    else
      SEA_chart$custom_adjP<-NA


  SEA_chart <- Filter(function(x) !(all(x=="")), SEA_chart)
  row.names(SEA_chart)<-NULL

  return(SEA_chart)

}
