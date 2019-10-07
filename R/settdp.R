#' @title setTDP
#'
#' @description Estimates the TDP of the specified set of features.
#'
#'
#' @param pvalue The vector of p-values. It can be the name of the covariate representing the Vector of
#' raw p-values in the \code{data} or a single vector but in the latter case it should match the
#' \code{featureIDs} vector
#'
#' @param featureIDs The vector of feature IDs. It can be the name of the covariate representing the IDs in the
#' \code{data} or a single vector but in the latter case it should match the \code{pvalue} vector
#'
#' @param data Optional data frame or matrix containing the variables in \code{pvalue} and \code{featureIDs}
#'
#' @param set The selection of features defining the feature-set based on the the \code{featureIDs}.
#' If missing, the set of all features is evaluated
#'
#' @param alpha The type I error allowed. The default is 0.05. NOTE: this shouls be consistent across the study
#'
#' @return A named vector including the lower bound and point estimate for the true discovery proportion (TDP)
#' of the specified test for the feature-set is returned.
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso \code{\link{setTest}}, \code{\link{SEA}}
#'
#'
#' @references
#' Mitra Ebrahimpoor, Pietro Spitali, Kristina Hettne, Roula Tsonaka, Jelle Goeman,
#' Simultaneous Enrichment Analysis of all Possible Gene-sets: Unifying Self-Contained
#' and Competitive Methods, Briefings in Bioinformatics, , bbz074, https://doi.org/10.1093/bib/bbz074
#'
#' @examples
#'
#' \dontrun{
#' set.seed(159)
#' #generate random p-values with pseudo IDs
#' m<- 100
#' pvalues <- runif(m,0,1)^5
#' featureIDs <- as.character(1:m)
#'
#' # perform a self-contained test for all features
#' settest(pvalues, featureIDs, testype = "selfcontained")
#'
#' # estimate the proportion of true discoveries among all m features
#' settdp(pvalues, featureIDs)
#'
#'# create a random pathway of size 60
#' randset=as.character(c(sample(1:m, 60)))
#'
#'
#' # estimate the proportion of true discoveries in a random set of size 50
#' settdp(pvalues, featureIDs, set=randset)
#'
#' }
#'
#' @export
#'
#' @importFrom hommel hommel tdp
#'
#'

setTDP = settdp <- function(pvalue, featureIDs, data, set, alpha=0.05){

  #save the call to function
  cltdp<-match.call()

  #check data

  if(missing(data)){
    allpval<-pvalue
    geneid<-featureIDs
  }else{
    if(is.matrix(data))  data <- data.frame(data)
    #evaluate pval and featureIDs, which may be one of the colnames of data
    allpval <- eval(cltdp$pvalue, data, parent.frame())
    geneid <- eval(cltdp$featureIDs, data, parent.frame())
  }

  if(length(allpval)!=length(geneid)) stop('The arguments pvalue and featureIDs should match!')


  #create a hommel object
  homobj<-hommel::hommel(allpval)

  if(missing(featureIDs) & !missing(set)) stop('You need featureIDs to define a set!')


  #set the initial values
  if (missing(set)) {
    n <- length(allpval)
    subpval <- allpval
    ix<-NA
  }
  else {
    ix<-which(geneid %in% set)
    subpval <- allpval[ix]
    n <- length(subpval)
  }


  if (any(is.na(subpval)))
    stop('The set of pvalues includes NA!')

  if (n == 0) {
    warning('No features in the set!')
    return(p=1, adjusted=1)
  }


  if(any(is.na(ix))){
    tdpbound <- hommel::tdp(homobj, alpha = alpha )
    tdpestimate <- hommel::tdp(homobj, alpha = 0.5)}

  else{
    tdpbound <- hommel::tdp(homobj, ix, alpha = alpha )
    tdpestimate <- hommel::tdp(homobj, ix, alpha = 0.5)}


    res<-list(TDP.bound=tdpbound, TDP.estimate=tdpestimate)

    return(res)
}
