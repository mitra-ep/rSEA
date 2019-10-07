#' @title setTest
#'
#' @description calculates the adjusted p-value for the local hypothesis as defined by \code{testtype}
#' and \code{testvalue}.
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
#' If missing, the set of all features is selected
#'
#' @param testype Character, type of the test: "selfcontained" or "competitive". Choosing the self-contained
#' option will automatically set the threshold to zero and the \code{testvalue} is ignored. Choosing the
#' competitive option without a \code{testvalue} will set the threshold to the overall estimated proportion
#' of true hypotheses
#'
#' @param testvalue Optional value to test against. Setting this value to c along with
#' \code{testype=="competitive"} will lead to testing the null hypothesis against a threshold c.
#' Note: this value needs to be a proportion
#'
#' @return The adjusted p-value of the specified test for the feature-set is returned.
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso \code{\link{setTDP}} \code{\link{SEA}}
#'
#' @references
#' Mitra Ebrahimpoor, Pietro Spitali, Kristina Hettne, Roula Tsonaka, Jelle Goeman,
#' Simultaneous Enrichment Analysis of all Possible Gene-sets: Unifying Self-Contained
#' and Competitive Methods, Briefings in Bioinformatics, , bbz074, https://doi.org/10.1093/bib/bbz074
#'
#' @examples
#'
#' \dontrun{
#' #Generate a vector of pvalues
#' set.seed(159)
#'
#' m<- 100
#' pvalues <- runif(m,0,1)^5
#' featureIDs <- as.character(1:m)
#'
#' # perform a self-contained test for all features
#' settest(pvalues, featureIDs, testype = "selfcontained")
#'
#' # create a random pathway of size 60
#' randset=as.character(c(sample(1:m, 60)))
#'
#' # perform a competitive test for the random pathway
#' settest(pvalues, featureIDs, set=randset, testype = "competitive")
#'
#' # perform a unified null hypothesis test against 0.2 for a set of size 50
#' settest(pvalues, featureIDs, set=randset, testype = "competitive", testvalue = 0.2 )
#'
#' }
#' @export
#'
#' @importFrom hommel hommel tdp localtest
#'
#'

setTest=settest <- function(pvalue, featureIDs, data, set, testype, testvalue){

  #save the call to function
  cltest<-match.call()

  #check data

  if(missing(data)){
    pval<-pvalue
    geneid<-featureIDs
  }else{
    if(is.matrix(data))  data <- data.frame(data)
    #evaluate pval and featureIDs, which may be one of the colnames of data
     pval <- eval(cltest$pvalue, data, parent.frame())
     geneid <- eval(cltest$featureIDs, data, parent.frame())
     }

  if(length(pval)!=length(geneid)) stop('The arguments pvalue and featureIDs should match!')


  #create a hommel object
  homobj<-hommel::hommel(pval)

      #set the initial values
      if (missing(set)) {
        n <- length(pval)
        ix<-NA
      }
      else {
        ix<-which(geneid %in% set)
        pval <- pval[ix]
        n <- length(pval)
      }


      #evaluate the test type and threshold value
      total_tdp<- ceiling(hommel::tdp(homobj, alpha=0.05)*n)/n

      if (any(is.na(pval)))
        stop('The set includes NA!')

      if (n == 0) {
        warning('No features in the set')
        return(p=1, adjusted=1)
      }

      if(testype!="selfcontained" & testype!="competitive")
        stop("Test type is not correctly specifed")


      if (testype=="selfcontained")
        thr<-0
      if (testype=="competitive" & !missing(testvalue))
        thr<-testvalue
      if (testype=="competitive" & missing(testvalue))
        thr<-total_tdp

      #Do the test
      if (any(is.na(ix)))
      adjustedp <- localtest(homobj, tdp=thr)
          else
      adjustedp <- localtest(homobj, ix, tdp=thr)

    #return the p-value
    res <- adjustedp
    names(res) <- c("adj.P-value")

   return(res)
 }
