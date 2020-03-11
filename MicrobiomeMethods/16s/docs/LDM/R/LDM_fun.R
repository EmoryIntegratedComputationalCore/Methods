gower = function (d, square=TRUE, center=TRUE) 
{
    
    a <- as.matrix(d)
    
    #--------------------------------------------------------------
    # squaring the matrix 
    #--------------------------------------------------------------
    
    if (square) a <- a^2
    
    #--------------------------------------------------------------
    # centering rows and columns of the (squared) matrix of d
    #--------------------------------------------------------------
    
    if (center) {

        a =  sweep(a, 1, rowMeans(a) )
        a = -sweep(a, 2, colMeans(a) )/2

    }
    
    return(a)
    
}# gower


fdr.Sandev = function(p.otu) {
    
    m = length(p.otu)
    
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))
    
    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)   
    qval.orig = qval.sort[mat]
    results = qval.orig
    return(results)
    
} # fdr.Sandev




calculate.dist <- function(otu.table, tree=NULL, dist.method="bray", 
                           binary=FALSE, rarefy=0, scale.otu.table=TRUE) {
    
    rowsum = rowSums(otu.table)
    if (min(rowsum)==0) {
        warning("There exists sample(s) with zero reads at every OTU!")
    }
    
    rowsum[which(rowsum==0)] = 1
    if (scale.otu.table) freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowsum ) )
    else freq.table <- otu.table
    
    if (grepl(dist.method, "wt-unifrac")) {
        dist <- GUniFrac(otu.table, tree, alpha=c(1))$unifrac[,,"d_1"]
    }  else if (grepl(dist.method, "hellinger")) {
        dist <- 0.5*dist( x=sqrt(freq.table), method='euclidean')
    } else {
        dist <- vegdist(x=freq.table, method=dist.method)
    }
    
    # if (dist.method==tolower("bray")) {
    #     dist <- 0.5*as.matrix( dist( x=freq.table, method='manhattan') )
    # } else if (dist.method==tolower("Hellinger")) {
    #     dist <- 0.5*as.matrix( dist( x=sqrt(freq.table), method='euclidean') )
    # } else if (dist.method==tolower("Jaccard")) {
    #     otu.rff <- Rarefy(otu.table)$otu.tab.rff
    #     dist <- vegdist(otu.rff, method="jaccard", binary=TRUE)
    # } else if (dist.method=="unwt-unifrac") {
    # dist <- GUniFrac(freq.table, tree, alpha=c(1))$unifrac[,,"d_UW"]
    
    dist <- as.matrix(dist)
    
    return(dist)
    
} # calculate.dist


#' Adjust data by covariates
#' 
#' This function produces adjusted distance matrix and otu table (if provided) 
#' after removing the effects of covariates (e.g., confounders). 
#' Observations with any missing data are removed.
#' 
#' @param formula a symbolic description of the covariate model with the form \code{ ~ model}, 
#' where \code{model} is specified in the same way as for \code{lm} or \code{glm}. For example, 
#' \code{~ a + b} specifies a model with the main effects of covariates \code{a} and \code{b}, and 
#' \code{~ a*b}, equivalently \code{~ a + b + a:b}, specifies a model with the main effects of 
#' \code{a} and \code{b} as well as their interaction.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the covariates. 
#' If not found in \code{data}, the covariates are taken from environment(formula), 
#' typically the environment from which \code{adjust.data.by.covariates} is called. 
#' The default is NULL.
#' @param otu.table the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' If provided, the adjusted (and column-centered) otu.table at both the frequency scale 
#' and arcsin-root-transformed frequency scale are outputted. If provided, 
#' it is also used for calculating the distance matrix unless the distance matrix is directly 
#' imported through \code{dist}.
#' The default is NULL.
#' @param tree a phylogenetic tree. Only used for calculating a
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of
#'   requested distance does not require a phylogenetic tree, or if the distance
#'   matrix is directly imported through \code{dist}. The default is NULL.
#' @param dist.method method for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package
#'  (i.e., "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", 
#'  "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis")
#'   as well as "hellinger" and "wt-unifrac". 
#'   The default is "bray". 
#'   For more details, see the \code{dist.method} argument in the \code{ldm} function.
#' @param dist a distance matrix. Can be either an object of class "dist" or "matrix".
#'   The elements of the distance matrix will be squared and then the matrix will be centered if the default choices 
#'   \code{square.dist=TRUE} and \code{center.dist=TRUE} are used. If \code{dist=NULL}, the distance matrix is 
#'   calculated from the \code{otu.table}, using the values of \code{dist.method} (and \code{tree} if required). 
#'   The default is NULL.
#' @param square.dist a logical variable indicating whether to square the 
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale the rows of the OTU table 
#'   for the frequency scale.  For count data, this corresponds to dividing by the library size to give 
#'   relative frequencies. The default is TRUE. 
#' @param center.otu.table a logical variable indicating whether to center the 
#'   columns of the OTU table. The OTU table should be centered if the distance 
#'   matrix has been centered. Applies to both frequency and transformed scales.  The default is TRUE.
#' @param freq.scale.only a logical variable indicating whether to provide adjusted frequency-scale data matrix only (not adjusted data on the arcsin-root transformed frequency scale). The default is FALSE.
#' @return a list consisting of 
#'   \item{adj.dist}{the (squared/centered) distance matrix
#'   after adjustment of covariates.}
#'   \item{x.freq}{the (column-centered) frequency-scale data matrix after adjustment of covariates.} 
#'   \item{x.tran}{the (column-centered) arcsin-root-transformed 
#'   data matrix after adjustment of covariates.} 

#' @keywords microbiome PCA ordination distance
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @examples
#' adj.data <- adjust.data.by.covariates(formula= ~ Sex + AntibioticUse, data=throat.meta,
#'                                       otu.table=throat.otu.tab, dist.method="bray")
#'
#' #-------------------------------------------------
#' # Use the adjusted distance matrix for ordination
#' #-------------------------------------------------
#'
#' PCs <- eigen(adj.data$adj.dist, symmetric=TRUE)
#' 
#' color = rep("blue", length(throat.meta$SmokingStatus))
#' w = which(throat.meta$SmokingStatus=="Smoker")
#' color[w] = "red"
#' 
#' plot(PCs$vectors[,1], PCs$vectors[,2], xlab="PC1", ylab="PC2", 
#'      col=color, main="Smokers vs. non-smokers")
#' legend(x="topleft", legend=c("smokers","non-smokers"), pch=c(21,21), 
#'        col=c("red","blue"), lty=0)



adjust.data.by.covariates = function(formula=NULL, data=NULL, 
                                     otu.table=NULL, tree=NULL, dist.method="bray", dist=NULL, 
                                     square.dist=TRUE, center.dist=TRUE, 
                                     scale.otu.table=TRUE, center.otu.table=TRUE,
                                     freq.scale.only=FALSE) {
    
    #------------------------
    # covariates (e.g., confounders)
    #------------------------
    
    options(na.action=na.pass)
    m1 = model.matrix(object=formula, data=data)
    
    if (!is.null(dist)) {
        dist = as.matrix(dist)
        if (dim(m1)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between covariates and dist' )
    }
    if (!is.null(otu.table)) {
        if (dim(m1)[1] != dim(otu.table)[1]) 
            otu.table <- t(otu.table)
        if (dim(m1)[1] != dim(otu.table)[1]) stop( 'numbers of observations mismatch between covariates and otu.table' )
    }
    
    if (any(is.na(m1))) {
        w = apply(is.na(m), 1, any)
        warning(paste(sum(w), 'observation(s) with any missing data are removed', sep=" "))
        
        w = !w
        if (!is.null(dist)) dist = dist[w, w]
        if (!is.null(otu.table)) otu.table = otu.table[w,]
        m1 = m1[w,]
    }
    
    center.m1 = TRUE
    if (center.m1) m1 = scale( m1, center=TRUE, scale=FALSE )
    

    #------------------------
    # dist matrix
    #------------------------
    
    if (is.null(dist) & is.null(otu.table)) {
        stop( 'must specify one of dist and otu.table' )
    }
    
    if (!is.null(otu.table) & !is.null(dist)) {
        if (dim(otu.table)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between otu.table and dist' )
    }
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.method=dist.method, otu.table=otu.table, tree=tree, scale.otu.table=scale.otu.table)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)

   
    #---------------------
    # calculate d.resid
    #---------------------
    
    tol.d=10^-8
    svd.m1 = svd(m1)
    use = (svd.m1$d>tol.d)
    
    hat.matrix = svd.m1$u[, use] %*% t( svd.m1$u[, use] )
    hat.matrix.bar = diag(dim(d.gower)[1]) - hat.matrix
    
    d.resid = hat.matrix.bar %*% d.gower
    d.resid = d.resid %*% hat.matrix.bar
    
    #---------------------
    # calculate adj.otu.table
    #---------------------
    
    x.freq = NULL
    x.tran = NULL
    
    if (!is.null(otu.table)) {
        
        rowsum = rowSums(otu.table)
        if (min(rowsum)==0) {
            warning("There exists sample(s) with zero reads at every OTU!")
        }
        
        # freq
        if (scale.otu.table) {
            rowsum[which(rowsum==0)] = 1
            freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowsum ) )
        } else {
            freq.table <- otu.table
        }
        x.freq <- scale( freq.table, center=center.otu.table, scale=FALSE )
        
        x.tran <- NULL
        if (!freq.scale.only) {
            # arcsin
            theta <- asin(sqrt(freq.table))
            x.tran <- scale( theta, center=center.otu.table, scale=FALSE)
        }
        
        # check discordant centering
        max.center.gower <- max(abs(rowMeans(d.gower)))
        max.center.xfreq  <- max(abs(rowMeans(x.freq)))
        if ( (max.center.gower - 10^-6) * (max.center.xfreq - 10^-6) < 0) {
            stop( 'discordant centering of the OTU table and distance matrix' )
        }
        
        # b.model
        b.model = svd.m1$u[, use]
        
        # adj.otu.table
        x1.tilda.freq = t(b.model) %*% x.freq
        x1.tilda.freq = b.model %*% x1.tilda.freq
        x.freq = x.freq - x1.tilda.freq
        
        if (!freq.scale.only) {
            x1.tilda.tran = t(b.model) %*% x.tran
            x1.tilda.tran = b.model %*% x1.tilda.tran
            x.tran = x.tran - x1.tilda.tran
        }
    }
    
    res <- list(x.freq=x.freq,
                x.tran=x.tran,
                adj.dist=d.resid)
    
    return(res)
    
} # adjust.data.by.covariates




#' Testing hypotheses using a linear decomposition model (LDM)
#' 
#' This function allows you to simultaneously test the global association with the overall  
#' microbiome composition and individual OTU associations to give coherent 
#' results. It is capable of handling complex design features such as 
#' confounders, interactions, and clustered data.
#' 
#' The formula has the form 
#' 
#' \code{otu.table ~ (first set of covariates) + (second set of covariates)
#' ... + (last set of covariates)} 
#' 
#' or 
#' 
#' \code{otu.table | confounders ~ (first set of covariates) + (second set of covariates)
#' ... + (last set of covariates)} 
#' 
#' where \code{otu.table} is
#' the OTU table with rows for samples and columns for OTUs and each set of 
#' covariates are enclosed in parentheses. The covariates in each submodel (set of covariates) are tested jointly,
#' after projecting off terms in submodels that appear earlier in the model.
#' 
#' For example, given OTU table \code{y} and a data frame \code{metadata} that contains 4 covariates, 
#' \code{a}, \code{b}, \code{c} and \code{d},  
#' some valid formulas would be:
#' 
#' \code{y ~ a + b + c + d} ### no confounders, 4 submodels (i.e., sets of covariates)
#' 
#' \code{y ~ (a+b) + (c+d)} ### no confounders, 2 submodels each having 
#' 2 covariates; 
#' 
#' \code{y | b ~ (a+c) + d} ### \code{b} is a confounder, submodel 1 is 
#' \code{(a+c)}, and submodel 2 is \code{d}
#' 
#' \code{y | b+c ~ a*d}     ### there are 2 confounders \code{b} 
#' and \code{c}; there is 1 submodel consisting of the three terms \code{a}, \code{d}, and \code{a:d} (interaction). 
#' This example is equivalent to \code{y | b+c ~ (a+d+a:d)}.
#' 
#' \code{y | as.factor(b) ~ (a+d) + a:d}  ### now confounder 
#' \code{b} will be treated as a factor variable, submodel 1 will have the main 
#' effects \code{a} and \code{d}, and submodel 2 will have only the interaction 
#' between \code{a} and \code{d}
#' 
#' \code{y | as.factor(b) ~ (a) + (d) + (a:d)} ### there are 3 submodels \code{a}, \code{d}, and \code{a:d}.
#' Putting paratheses around a single variable is allowed but not necessary.
#'
#' Submodels that combine character and numeric values are allowed; character-valued variables are coerced into factor 
#' variables.  Confounders are distinguished from other covariates as test statistics are not calculated for confounders
#' (which are included for scientific reasons, not by virtue of significance test results); 
#' consequently they also do not contribute to stopping criteria.  If tests of confounders are desired, confounders should
#' put on the right hand side of the formula as the first submodel.
#' 
#' LDM uses two sequential stopping criteria. For the global test, LDM uses the 
#' stopping rule of Besag and Clifford (1991), which stops permutation when a 
#' pre-specified minimum number (default=20) of rejections (i.e., the permutation 
#' statistic exceeded the observed test statistic) has been reached. For the 
#' OTU-specific tests, LDM uses the stopping rule of Sandve et al. (2011), 
#' which stops permutation when every OTU test has either reached the pre-specified 
#' number (default=20) of rejections or yielded a q-value that is below the 
#' nominal FDR level (default=0.1). As a convention, we call a test "stopped"
#' if the corresponding stopping criterion has been satisfied. Although all tests 
#' are always terminated if a pre-specified maximum number (see description of \code{n.perm.max} in Arguments list) of 
#' permutations have been generated, some tests may not have "stopped."  This typically occurs when
#' the relevant p-value is small or near the cutoff for inclusion in a list of significant findings; 
#' for global tests meeting the stopping criterion is not critical, but 
#' caution is advised when interpreting OTU-level tests that have not stopped as additional OTUs may be found 
#' with a larger number of permutations.
#' 
#' 
#' @param formula a symbolic description of the 
#'   model to be fitted. The details of model specification are given under 
#'   "Details".
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the covariates of interest and 
#' confounding covariates. 
#' If not found in \code{data}, the covariates are taken from environment(formula), 
#' typically the environment from which \code{ldm} is called. The default is .GlobalEnv.
#' @param tree a phylogenetic tree. Only used for calculating a 
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of 
#'   the requested distance does not involve a phylogenetic tree, or if the 
#'   distance matrix is directly imported through \code{dist}.
#' @param dist.method method for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package
#'  (i.e., "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", 
#'  "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
#'  "chao", "cao", "mahalanobis") as well as "hellinger" and "wt-unifrac". 
#'  Note that the option \code{binary} is set to \code{FALSE} 
#'  (i.e., not allowing calculation of presence-absence distances) 
#'  when we internally call \code{vegdist} for its supported distances. 
#'  The Hellinger distance measure (\code{dist.method="hellinger"}) takes the form
#'  \code{0.5*E}, where E is the Euclidean distance between the square-root-transformed 
#'  frequency data. The weighted UniFrac distance (\code{dist.method="wt-unifrac"}) 
#'  is calculated by interally calling \code{GUniFrac} in the \code{GUniFrac} package.
#'   Not used when anything other than \code{dist=NULL} is specified for \code{dist}.
#'   The default is "bray".
#' @param dist a distance matrix. Can be an object of class either "dist" or "matrix".
#'   The elements of the distance matrix will be squared and then the matrix will be centered if the default choices 
#'   \code{square.dist=TRUE} and \code{center.dist=TRUE} are used. If \code{dist=NULL}, the distance matrix is 
#'   calculated from the \code{otu.table}, using the values of \code{dist.method} (and \code{tree} if required). 
#'   The default is NULL.
#' @param cluster.id character or factor variable that identifies clusters. The default value
#'   cluster.id=NULL if the observations are not clustered (i.e., are independent).
#' @param perm.within.type a character string that takes values "free", "none", "series", or "grid".  
#'   The default is "free" (for random permutations).
#' @param perm.within.nrow a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See documentation for permute package for additional details
#' @param perm.within.ncol a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See documentation for permute package for additional details
#' @param perm.between.type a character string that takes values "free", "none", or "series".  
#'   The default is "none".
#' @param strata a character or factor variable that defines strata (groups), within which to constrain permutations. 
#'   The default is NULL.
#' @param how a permutation control list, for users who want to specify their own call to the \code{how} function from the \code{permute} package.  
#'   The default is NULL.
#' @param n.perm.max the maximum number of permutations.  The default is NULL, in which case a maximum of
#'   5000 permutations are used for the global test and a maximum of \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal}) 
#'   are used for the OTU test, where \code{n.otu} is the number of OTUs.  If a numeric value for \code{n.otu} is specified, 
#'   this value is used for both global and OTU-level tests.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. The 
#'   default is 20.
#' @param seed an integer seed for the random number generator in the 
#'   permutation procedure. The default is NULL; with the default value, an integer seed will be 
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case 
#'   the user wants to reproduce the permutations.
#' @param test.global a logical value indicating whether to perform the global 
#'   test. The default is TRUE.
#' @param test.otu a logical value indicating whether to perform the 
#'   OTU-specific tests. The default is TRUE.
#' @param fdr.nominal the nominal FDR value. The default is 0.1.
#' @param square.dist a logical variable indicating whether to square the 
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale the rows of the OTU table for the freq scale.  For count 
#'   data, this corresponds to dividing by the library size to give relative frequencies.  Does not affect the tran scale.  
#'   The default is TRUE. 
#' @param center.otu.table a logical variable indicating whether to center the 
#'   columns of the OTU table. The OTU table should be centered if the distance 
#'   matrix has been centered. Applies to both the frequency and transformed scales.  The default is TRUE.
#' @param freq.scale.only a logical variable indicating whether to perform analysis of the frequency-scale data only (not the arcsin-root transformed frequency data and 
#' the omnibus test). The default is FALSE.
#' @return a list consisting of 
#'   \item{b}{the matrix B as defined in Hu and Satten (2018)} 
#'   \item{dist}{the (squared/centered) distance matrix} 
#'   \item{x.freq}{the frequency-scale data matrix, scaled and centered if so specified} 
#'   \item{d.freq}{a vector of the nonnegative diagonal elements of \code{D} that satisfies
#'   \code{b^T x.freq = D v^T}}
#'   \item{v.freq}{the v matrix with unit columns that satisfies
#'   \code{b^T x.freq = D v^T}}
#'   \item{x.tran}{the (column-centered) arcsin-root-transformed 
#'   data matrix} 
#'   \item{d.tran}{a vector of the nonnegative diagonal elements of \code{D} that satisfies
#'   \code{b^T x.tran = D v^T}}
#'   \item{v.tran}{the v matrix with unit columns that satisfies
#'   \code{b^T x.tran = D v^T}}
#'   \item{low}{a vector of lower indices for confounders (if there is any) and submodels}
#'   \item{up}{a vector of upper indices for confounders (if there is any) and submodels}
#'   \item{VE.global.freq.confounders}{Variance explained (VE) by confounders, based on the frequency-scale data}
#'   \item{VE.global.freq.submodels}{VE by each submodel, based on the frequency-scale data}
#'   \item{VE.global.freq.residuals}{VE by each component in the residual distance, based on the frequency-scale data}
#'   \item{VE.otu.freq.confounders}{Contribution of each OTU to VE by confounders, based on the frequency-scale data}
#'   \item{VE.otu.freq.submodel}{Contribution of each OTU to VE by each submodel, based on the frequency-scale data}
#'   \item{VE.global.tran.confounders}{Variance explained (VE) by confounders, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.global.tran.submodels}{VE by each submodel, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.global.tran.residuals}{VE by each component in the residual distance, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.otu.tran.confounders}{Contribution of each OTU to VE by confounders, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.otu.tran.submodels}{Contribution of each OTU to VE by each submodel, based on 
#'   the arcsin-root-transformed frequency data}
#'   \item{VE.df.confounders}{Degree of freedom (i.e., number of components) associated with the VE for confounders}
#'   \item{VE.df.submodels}{Degree of freedom (i.e., number of components) associated with the VE for each submodel}
#'   \item{F.global.freq}{F statistics for testing each submodel, based on
#'   the frequency-scale data} 
#'   \item{F.global.tran}{F statistics for testing each submodel, based on 
#'   the arcsin-root-transformed frequency data} 
#'   \item{F.otu.freq}{F statistics for testing each OTU for each submodel, based on the frequency-scale data} 
#'   \item{F.otu.tran}{F statistics for testing each OTU for each submodel, based on the arcsin-root-transformed data} 
#'   \item{p.global.freq}{p-values for the global test of each set of covariates
#'   based on the frequency-scale data} 
#'   \item{p.global.tran}{p-values for the global test of each set of covariates
#'   based on the arcsin-root-transformed frequency data} 
#'   \item{p.global.omni}{p-values for the global test of each set of covariates 
#'   based on the omnibus statistics, which are the minima of the p-values obtained 
#'   from the frequency scale and the arcsin-root-transformed frequency data 
#'   as the final test statistics, and use the corresponding minima from the 
#'   permuted data to simulate the null distributions} 
#'   \item{p.otu.freq}{p-values for the OTU-specific tests based on the 
#'   frequency scale data} 
#'   \item{p.otu.tran}{p-values for the OTU-specific tests based on the 
#'   arcsin-root-transformed frequency data} 
#'   \item{p.otu.omni}{p-values for the OTU-specific tests based on the 
#'   omnibus statistics} 
#'   \item{q.otu.freq}{q-values (i.e., FDR-adjusted p-values) 
#'   for the OTU-specific tests based on the frequency scale data} 
#'   \item{q.otu.tran}{q-values for the OTU-specific tests based on 
#'   the arcsin-root-transformed frequency data} 
#'   \item{q.otu.omni}{q-values for the OTU-specific tests based on the 
#'   omnibus statistics} 
#'   \item{n.perm.completed}{number of permutations completed} 
#'   \item{global.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all global tests} 
#'   \item{otu.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all OTU-specific tests}
#'   \item{seed}{a single-value integer seed that is user supplied or internally generated}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @references Hu YJ, Satten GA. (2019) Testing hypotheses about microbiome 
#'   using the linear decomposition model. 
#'   bioRXiv:doi.org/10.1101/229831.
#' @examples
#' 
#'#-----------------------------------------------
#'# fit only
#'#-----------------------------------------------
#'fit <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
#'           data=throat.meta, dist.method="bray", n.perm.max=0)
#'
#'#-----------------------------------------------
#'# test the global hypothese only
#'#-----------------------------------------------
#'res1.ldm <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
#'                data=throat.meta, dist.method="bray", 
#'                test.global=TRUE, test.otu=FALSE, seed=123)
#'                     
#'#----------------------------------------------------
#'# test both the global hypothese and individual OTUs
#'#----------------------------------------------------
#'res2.ldm <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
#'                data=throat.meta, dist.method="bray", 
#'                test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1, seed=123)
#' 
#'#----------------------------------------------------
#'# clustered data
#'#----------------------------------------------------
#'res4.ldm <- ldm(formula=sim.otu.tab | X ~ Y, data=sim.meta, dist.method="bray", 
#'                cluster.id=ID, perm.between.type="free", perm.within.type="none",
#'                test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1, seed=123)


ldm = function( formula, data=.GlobalEnv, tree=NULL, dist.method="bray", dist=NULL, 
                     cluster.id=NULL, strata=NULL, how=NULL,
                     perm.within.type="free", perm.between.type="none",
                     perm.within.ncol=0, perm.within.nrow=0,
                     n.perm.max=NULL, 
                     n.rej.stop=20, seed=NULL, 
                     test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
                     square.dist=TRUE, center.dist=TRUE, 
                     scale.otu.table=TRUE, center.otu.table=TRUE,
                     freq.scale.only=FALSE) {  

    
    #------------------------
    # form.call
    #------------------------
    
    object=formula
    #
    #   extract cluster.id from dataframe
    #
    cl=match.call()
    mf=match.call(expand.dots=FALSE)
    m=match( x='cluster.id', table=names(mf) )
    mf.string=as.character( mf[c(1L,m)] )
    cluster.name=mf.string[2]
    if (cluster.name=='NULL') {
        cluster.id=NULL
    } else {   
        loc.dollar=tail( gregexpr('\\$', cluster.name)[[1]] , n=1 )
        if (loc.dollar<0)  {
            cluster.id=getElement(data,cluster.name)
            if( is.null(cluster.id) ) cluster.id=get(cluster.name)
        } else {   
            df.name=get( substr(cluster.name, start=1, stop=loc.dollar-1) )
            var.name=substr(cluster.name, start=loc.dollar+1, stop=nchar(cluster.name))            
            cluster.id= getElement(df.name,var.name) 
        }
    }
    #        
    #   extract model from formula    
    #    
    obj=toString(object)
    obj=gsub('\\s','',obj)
    prefix=' ~ + 0 + '
    loc.comma=gregexpr(',',obj)[[1]]
    start.terms=loc.comma[2]
    terms=substr(obj,start=start.terms+1, stop=nchar(obj))
    #
    #   find n.obs and full set of rownames
    #   
    if (class(data)=='data.frame') {
        row.names=rownames(data)
        n.obs=length(row.names)
    } else {   
        df=model.frame( paste('~',terms) , na.action=na.pass )
        row.names=rownames(df)
        n.obs=length(row.names)
    }
    #
    #   check for missing values in cluster.id
    #        
    
    if (is.null(cluster.id)) {
        use.rows=row.names
    } else {   
        use=!is.na(cluster.id)
        use.rows=row.names[use]
    }
    #
    #   check for and extract confounders
    #
    model=list()
    j=1
    loc.bar=regexpr('\\|',obj)[1]
    loc.minus=regexpr('-',obj)[1]
    loc.delim=max( loc.bar, loc.minus)
    if (loc.delim>0) {
        end.confound=loc.comma[2]
        c=substr(obj,start=loc.delim+1, stop=end.confound-1)
        conf=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        model[[j]]=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
#       use.rows=intersect( use.rows, rownames(conf) )
        use.rows=rownames(model[[1]]) 
        j=j+1
    } else {
        conf=NULL
    }     
    #
    #   extract model terms
    #
#   j=1
    continue=TRUE
    while (continue) {
        if (substr(terms,1,1)=='(') {
            stop=regexpr(')\\+',terms)[1]
        } else {
            stop=regexpr('\\+',terms)[1] - 1
        }          
        
        if (stop<=0) stop=nchar(terms) 
        m=substr(terms, start=1, stop=stop)
        model[[j]]=model.matrix( as.formula( paste(prefix,m) ) , data=data)
        use.rows=intersect( use.rows, rownames(model[[j]]) )
        #        if (j==1) {
        #            use.rows=rownames(model[[1]])
        #            }
        #        else {
        #            use.rows=intersect( use.rows, rownames(model[[j]]) )
        #            }         
        if (stop+2<=nchar(terms)) {
            terms=substr(terms, start=stop+2, stop=nchar(terms))
            j=j+1
        } else {
            continue=FALSE
        }             
    }   
    n.model=j    
    #
    #  extract OTU table
    #      
    if (is.null(conf)) loc.delim=loc.comma[2]
    otu.name=substr(obj, start=loc.comma[1]+1, stop=loc.delim-1)
    #   loc.dollar=regexpr('\\$', otu.name)[1]
    loc.dollar=tail( gregexpr('\\$', otu.name)[[1]] , n=1 )
    if (loc.dollar<0)  {
        if (class(data)=='data.frame') {
            otu.table=getElement(data, otu.name)
            if (is.null(otu.table)) otu.table= get(otu.name) 
            otu.table=as.matrix(otu.table)
        } else {
            otu.table=as.matrix( get(otu.name) )
        }
    } else {
        df.name=get( substr(otu.name, start=1, stop=loc.dollar-1) )
        var.name=substr(otu.name, start=loc.dollar+1, stop=nchar(otu.name))
        otu.table=as.matrix( getElement(df.name,var.name) )
    }        
    #    if (is.null(otu.table)) otu.table=as.matrix( getElement(.GlobalEnv,otu.name) )
    if ( nrow(otu.table) != n.obs ) {
        if (ncol(otu.table)==n.obs ) {
            otu.table=t(otu.table)
        } else {   
            print('warning: OTU table and covariates have different number of observations')
            return
        }
    }   
    #
    #   remove rows having NA 
    #    
    for (j in 1:n.model) {
        keep =  rownames( model[[j]] ) %in% use.rows
        model[[j]]=model[[j]][keep,,drop=FALSE]
    }
    if (!is.null(conf)) {
        keep =  rownames(conf) %in% use.rows 
        conf=conf[keep,,drop=FALSE]
    }
    keep=row.names %in% use.rows    
    otu.table=otu.table[keep,,drop=FALSE]    
    if (!is.null(cluster.id)) cluster.id=cluster.id[keep]
    
    # transpose
    if (dim(model[[1]])[1] != dim(otu.table)[1]) 
        otu.table <- t(otu.table)
    if (dim(model[[1]])[1] != dim(otu.table)[1]) stop( 'numbers of observations mismatch between covariates and the OTU table' )
    
    # OTU names
    if (is.null(colnames(otu.table))) { 
        colnames(otu.table) = 1:ncol(otu.table)
    }
    
    # remove zero OTUs
    w = which(colSums(otu.table)>0)
    otu.table = otu.table[,w]
    
    
    #------------------------
    # setup model
    #------------------------
    
    n.obs = dim(model[[1]])[1]
    adjust.for.confounders = !is.null(conf)
    
    n.var = length(model)
    n.var1 = n.var - as.numeric(adjust.for.confounders)
    
    center.vars=TRUE
    
    index = rep(0, n.var)
    
    for (i in 1:n.var) {
        m.i = model[[i]]
        if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
        
        if (i==1) {
            m = m.i
            index[i] = dim(m.i)[2] 
        } else {
            m = cbind(m, m.i)   
            index[i] = index[i-1] + dim(m.i)[2]    
        }
        
    }
    
    
    #------------------------
    # setup permutation
    #------------------------
    
    if (class(how)=='how') {
        CTRL=how                   # user-provided how list
    } else {
        if (is.null(cluster.id)) {
            if (is.null(perm.within.type) & is.null(perm.between.type)) {
                # default when no unclustered data has no type specified is 'free'
                perm.within.type='free'    
            }
            if (is.null(strata)) {
                # setup for unclustered permutation
                CTRL = how( within=Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))  
            } else {
                # setup for unclustered, stratified permutation
                strata=as.factor(strata)
                CTRL = how( blocks=strata, within=Within(type=perm.within.type, 
                                                         nrow=perm.within.nrow, 
                                                         ncol=perm.within.ncol))  
            }    
        } else {        
            cluster.id=as.factor(cluster.id)
            if (is.null(strata)) {            
                #  clustered but unstratified data
                CTRL = how( plots=Plots(cluster.id, type=perm.between.type ), 
                            within=Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            } else {
                #   clustered and stratified data
                strata=as.factor(strata)             
                CTRL = how( blocks=strata, 
                            plots=Plots(cluster.id, type=perm.between.type ), 
                            within=Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            }
        }
    }    
    
    
    #------------------------
    # dist matrix
    #------------------------
    
    if (!is.null(dist)) {
        dist <- as.matrix(dist)
        if (dim(otu.table)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between the OTU table and dist' )
    }
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.method=dist.method, otu.table=otu.table, tree=tree, scale.otu.table=scale.otu.table)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)
    
    
    #------------------------
    # data matrix X
    #------------------------
    
    rowsum = rowSums(otu.table)
    if (min(rowsum)==0) {
        warning("There exists sample(s) with zero reads at every OTU!")
    }
    
    
    # freq
    if (scale.otu.table) {
        rowsum[which(rowsum==0)] = 1
        freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowsum ) )
    } else {
        freq.table <- otu.table
    }
    x.freq <- scale( freq.table, center=center.otu.table, scale=FALSE )
    
    x.tran <- NULL
    if (!freq.scale.only) {
        # arcsin
        theta <- asin(sqrt(freq.table))
        x.tran <- scale( theta, center=center.otu.table, scale=FALSE)
    }
    
    
    max.center.gower <- max(abs(rowMeans(d.gower)))
    max.center.xfreq  <- max(abs(rowMeans(x.freq)))
    if ( (max.center.gower - 10^-6) * (max.center.xfreq - 10^-6) < 0) {
        stop( 'discordant centering of the OTU table and distance matrix' )
    }
    
    #---------------------
    # model fitting
    #---------------------
  
    fit.ldm = calculate.b.and.resid( d.gower=d.gower, x.freq=x.freq, x.tran=x.tran, 
                                     index=index, m=m, adjust.for.confounders=adjust.for.confounders)  

    #---------------------
    # observed statistic
    #---------------------
    
    ldm.obs.freq = ldm.stat(b=fit.ldm$b, low=fit.ldm$low, up=fit.ldm$up, resid=fit.ldm$resid.freq, ss.tot=fit.ldm$ss.tot.freq, adjust.for.confounders=adjust.for.confounders)
    ldm.obs.tran = NULL
    if (!freq.scale.only) ldm.obs.tran = ldm.stat(b=fit.ldm$b, low=fit.ldm$low, up=fit.ldm$up, resid=fit.ldm$resid.tran, ss.tot=fit.ldm$ss.tot.tran, adjust.for.confounders=adjust.for.confounders)
    

    p.global.freq = NULL
    p.global.tran = NULL
    p.global.omni = NULL
    
    p.otu.freq = NULL
    p.otu.tran = NULL
    p.otu.omni = NULL
    
    q.otu.freq = NULL
    q.otu.tran = NULL
    q.otu.omni = NULL
    
    
    n.perm.completed = NULL
    n.global.perm.completed = NULL
    n.otu.perm.completed = NULL
    
    global.tests.stopped = NULL
    otu.tests.stopped    = NULL
    
    if (is.null(n.perm.max)) {
        n.global.perm.max = ifelse(test.global, 5000, NA)
        n.otu.perm.max = ifelse(test.otu, ncol(x.freq) * n.rej.stop * (1/fdr.nominal), NA)
        n.perm.max = max(n.global.perm.max, n.otu.perm.max, na.rm=TRUE)
    } else {
        n.global.perm.max = ifelse(test.global, n.perm.max, NA)
        n.otu.perm.max = ifelse(test.otu, n.perm.max, NA)
    }
        
    if (n.perm.max > 0) {
        
        
        n.pvalue.cal.min = 1000
        n.pvalue.cal.step = 100
        n.perm.block = 1000
        
        n.stable.max=10 # n.pvalue.cal.step*n.stable.max=1000
        
        
        if (is.null(seed)) {
            seed = sample(1:10^6, 1)
        }
        
        n.otu = dim(x.freq)[2]
        
        otu.smallp = 1:n.otu
        n.otu.smallp = n.otu
        
        if (test.global) 
        {
            global.tests.stopped = FALSE
            
            global.freq = ldm.obs.freq$ve.global
            global.freq.perm = matrix(NA, n.var1, n.global.perm.max)
            n.global.freq = 0
            n.global.tran = n.global.perm.max # to satisfy the "if" condition
            if (!freq.scale.only) {
                global.tran = ldm.obs.tran$ve.global
                global.tran.perm = matrix(NA, n.var1, n.global.perm.max)
                n.global.tran = 0
            }
    
        }    
        
      
        if (test.otu)
        {
            otu.tests.stopped  = FALSE
            otu.freq = ldm.obs.freq$ve.otu
            
            n.perm.tmp = ifelse(test.global, max(n.global.perm.max, n.pvalue.cal.min), n.pvalue.cal.min)
            
            otu.freq.perm = array(NA, c(n.var1, n.otu, n.perm.tmp))
            
            n.otu.freq = 0
            
            Aset.freq = matrix(TRUE, n.var1, n.otu)

            p.otu.freq = matrix(NA, n.var1, n.otu)
            p.otu.tran = NULL
            p.otu.omni = NULL
            
            if (!freq.scale.only) {
                otu.tran = ldm.obs.tran$ve.otu
                otu.tran.perm = array(NA, c(n.var1, n.otu, n.perm.tmp))
                n.otu.tran = 0
                
                Aset.tran = matrix(TRUE, n.var1, n.otu)
                Aset.omni = matrix(TRUE, n.var1, n.otu)
                
                p.otu.tran = matrix(NA, n.var1, n.otu)
                p.otu.omni = matrix(NA, n.var1, n.otu)
            }
            
        }
      
        tol.eq = 10^-16
        
        n.perm.completed = 0
        n.global.perm.completed = n.global.perm.max
        n.otu.perm.completed = n.otu.perm.max
        
        n.stable.freq = 0
        if (!freq.scale.only) {
            n.stable.tran = 0
            n.stable.omni = 0
        }
      
        set.seed(seed)
        
        for (i.sim in 1:n.perm.max) {
            
            i.sim.r = i.sim%%n.perm.block
            if (i.sim.r==0) i.sim.r = n.perm.block
            
            if (i.sim.r==1) {
                cat("permutations:", i.sim, "\n")
                perm = shuffleSet(n.obs, n.perm.block, CTRL)
            }
            
            
            n.perm.completed = n.perm.completed + 1
            inv.n.perm.completed = 1/n.perm.completed
            inv.n.perm.completed.1 = 1/(n.perm.completed+1)
            
            #------------------------------------------------
            # reduce OTUs to speed-up computation
            #------------------------------------------------
            
            if (i.sim > 100000) {
                n.pvalue.cal.min = 10000
                n.pvalue.cal.step = 1000
            }
            
            if (test.otu & i.sim > n.pvalue.cal.min & i.sim %% n.pvalue.cal.min==1){ # 1001, 2001, 3001, ...
                
                global.done = TRUE
                if (test.global) if (!global.tests.stopped & i.sim <= n.global.perm.max) global.done = FALSE
                
                if (global.done) {
                
                    if (!freq.scale.only) {
                        w.otu.smallp = which(apply(rbind(Aset.freq, Aset.tran, Aset.omni), 2, any))
                    } else {
                        w.otu.smallp = which(apply(Aset.freq, 2, any))
                    }
                    
                    cat("number of OTUs do not meet early stopping criterion:", length(w.otu.smallp), "\n")
                    
                    if (length(w.otu.smallp) != n.otu.smallp) {
                    
                        otu.smallp = otu.smallp[w.otu.smallp]
                        n.otu.smallp = length(otu.smallp)
                        
                        n.otu.freq = n.otu.freq[,w.otu.smallp,drop=FALSE]
                        Aset.freq = Aset.freq[,w.otu.smallp,drop=FALSE]
                        if (!freq.scale.only) {
                            n.otu.tran = n.otu.tran[,w.otu.smallp,drop=FALSE]
                            Aset.tran = Aset.tran[,w.otu.smallp,drop=FALSE]
                            Aset.omni = Aset.omni[,w.otu.smallp,drop=FALSE]
                        }
                    }
                    
                    otu.freq.perm.smallp = otu.freq.perm[,w.otu.smallp,1:(i.sim-1),drop=FALSE]
                    otu.freq.perm = array(NA, c(n.var1, n.otu.smallp, i.sim-1+n.pvalue.cal.min))
                    otu.freq.perm[,,1:(i.sim-1)] = otu.freq.perm.smallp
                    rm(otu.freq.perm.smallp)
                    
                    if (!freq.scale.only) {
                        otu.tran.perm.smallp = otu.tran.perm[,w.otu.smallp,1:(i.sim-1),drop=FALSE]
                        otu.tran.perm = array(NA, c(n.var1, n.otu.smallp, i.sim-1+n.pvalue.cal.min))
                        otu.tran.perm[,,1:(i.sim-1)] = otu.tran.perm.smallp
                        rm(otu.tran.perm.smallp)
                    }
                    
                }
            }


            # perform permutations                   
        
            b.perm = fit.ldm$b[perm[i.sim.r,], ]
            
            ldm.perm.freq = ldm.stat(b=b.perm, low=fit.ldm$low, up=fit.ldm$up, 
                                   resid=fit.ldm$resid.freq[,otu.smallp,,drop=FALSE], ss.tot=fit.ldm$ss.tot.freq[,otu.smallp,drop=FALSE], adjust.for.confounders=adjust.for.confounders)
            if (!freq.scale.only) ldm.perm.tran = ldm.stat(b=b.perm, low=fit.ldm$low, up=fit.ldm$up, 
                                   resid=fit.ldm$resid.tran[,otu.smallp,,drop=FALSE], ss.tot=fit.ldm$ss.tot.tran[,otu.smallp,drop=FALSE], adjust.for.confounders=adjust.for.confounders)

            if (test.global & i.sim <= n.global.perm.max) {
                if (!global.tests.stopped) {
                    global.freq.perm[,i.sim] = ldm.perm.freq$ve.global
                    n.global.freq = n.global.freq + (ldm.perm.freq$ve.global > global.freq + tol.eq) + (ldm.perm.freq$ve.global > global.freq - tol.eq)
                    if (!freq.scale.only) {
                        global.tran.perm[,i.sim] = ldm.perm.tran$ve.global
                        n.global.tran = n.global.tran + (ldm.perm.tran$ve.global > global.tran + tol.eq) + (ldm.perm.tran$ve.global > global.tran - tol.eq)
                    }
                    
                }
            }
        
            if (test.otu) {
                if (!otu.tests.stopped) {
  
                    otu.freq.perm[,,i.sim] = ldm.perm.freq$ve.otu
                    n.otu.freq = n.otu.freq + (ldm.perm.freq$ve.otu>otu.freq[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.freq$ve.otu>otu.freq[,otu.smallp,drop=FALSE]-tol.eq)
                    if (!freq.scale.only) {
                        otu.tran.perm[,,i.sim] = ldm.perm.tran$ve.otu
                        n.otu.tran = n.otu.tran + (ldm.perm.tran$ve.otu>otu.tran[,otu.smallp,drop=FALSE]+tol.eq) + (ldm.perm.tran$ve.otu>otu.tran[,otu.smallp,drop=FALSE]-tol.eq)
                    }
                }
            }
        
        
            if ((i.sim %% n.pvalue.cal.min == 0) 
                | (i.sim > n.pvalue.cal.min & i.sim %% n.pvalue.cal.step == 0) 
                | (i.sim==n.global.perm.max) 
                | (i.sim==n.perm.max)) {
          
                if (test.global & i.sim <= n.global.perm.max) {
                    if (!global.tests.stopped) {
                        
                        ################
                        # test global
                        ################
    
                        if ((all(n.global.freq >= n.rej.stop*2) & all(n.global.tran >= n.rej.stop*2)) 
                            | (i.sim==n.global.perm.max) 
                            | (i.sim==n.perm.max)) {
                            
                            p.global.freq = ifelse((n.global.freq >= n.rej.stop*2), 0.5*n.global.freq*inv.n.perm.completed, (0.5*n.global.freq+1)*inv.n.perm.completed.1)

                            if (!freq.scale.only) {
                                p.global.tran = ifelse((n.global.tran >= n.rej.stop*2), 0.5*n.global.tran*inv.n.perm.completed, (0.5*n.global.tran+1)*inv.n.perm.completed.1)
                            
                                p.global.freq.tmp <- 0.5*n.global.freq # omit (1/n.perm.completed) for computation efficiency
                                p.global.tran.tmp <- 0.5*n.global.tran 
                  
                                pmin.global.omni <- pmin(p.global.freq.tmp, p.global.tran.tmp)
                                pnull.global.freq <- (n.perm.completed + 0.5 - apply(global.freq.perm[,1:n.perm.completed,drop=FALSE], 1, rank))
                                pnull.global.tran <- (n.perm.completed + 0.5 - apply(global.tran.perm[,1:n.perm.completed,drop=FALSE], 1, rank))
                                pnullmin.global.omni <- pmin(pnull.global.freq, pnull.global.tran)
                  
                                n.global.omni <- apply((t(pnullmin.global.omni) < pmin.global.omni + tol.eq), 1, sum) + 0.5*apply((abs(t(pnullmin.global.omni) - pmin.global.omni) < tol.eq), 1, sum)  
                                p.global.omni = ifelse((n.global.omni >= n.rej.stop), n.global.omni*inv.n.perm.completed, (n.global.omni+1)*inv.n.perm.completed.1)

                                if (all(n.global.freq >= n.rej.stop*2) 
                                    & all(n.global.tran >= n.rej.stop*2)
                                    & all(n.global.omni >= n.rej.stop)) {
                                    global.tests.stopped = TRUE
                                    cat("global test stopped at permutation", i.sim, "\n")
                                    n.global.perm.completed = i.sim
                                }
                            } else {
                                if (all(n.global.freq >= n.rej.stop*2)) {
                                    global.tests.stopped = TRUE
                                    cat("global test stopped at permutation", i.sim, "\n")
                                    n.global.perm.completed = i.sim
                                }
                            }
                            
                        } 
                    }
                } 
          
                if (test.otu) {
                    if (!otu.tests.stopped) {
                        
                        ################
                        # test otu
                        ################
                
                        if (any(Aset.freq)) {
                            AtoB.freq <- Aset.freq & (n.otu.freq >= n.rej.stop*2)
                            Aset.freq <- Aset.freq & !AtoB.freq
                            p.otu.freq[,otu.smallp][AtoB.freq] <- 0.5*n.otu.freq[AtoB.freq]*inv.n.perm.completed
                            p.otu.freq[,otu.smallp][Aset.freq] <- (0.5*n.otu.freq[Aset.freq]+1)*inv.n.perm.completed.1
                        
                            q.otu.freq <- t(apply(p.otu.freq, 1, fdr.Sandev))
                        
                            Aset.freq.meet.criteria <- apply(((q.otu.freq[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.freq) | (!Aset.freq), 1, all)
                            n.stable.freq <- ifelse(Aset.freq.meet.criteria, n.stable.freq + Aset.freq.meet.criteria, 0)
                            Aset.freq.rm.row <- (n.stable.freq >= n.stable.max)
                            Aset.freq[Aset.freq.rm.row,] = FALSE
                        }
                
                        if (!freq.scale.only) {
                            if (any(Aset.tran)) {
                                AtoB.tran <- Aset.tran & (n.otu.tran >= n.rej.stop*2)
                                Aset.tran <- Aset.tran & !AtoB.tran
                                p.otu.tran[,otu.smallp][AtoB.tran] <- 0.5*n.otu.tran[AtoB.tran]*inv.n.perm.completed
                                p.otu.tran[,otu.smallp][Aset.tran] <- (0.5*n.otu.tran[Aset.tran]+1)*inv.n.perm.completed.1
                          
                                q.otu.tran <- t(apply(p.otu.tran, 1, fdr.Sandev))
                          
                                Aset.tran.meet.criteria <- apply(((q.otu.tran[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.tran) | (!Aset.tran), 1, all)
                                n.stable.tran <- ifelse(Aset.tran.meet.criteria, n.stable.tran + Aset.tran.meet.criteria, 0)
                                Aset.tran.rm.row <- (n.stable.tran >= n.stable.max)
                                Aset.tran[Aset.tran.rm.row,] = FALSE
                            }
                        
                            if ((!any(Aset.freq) & !any(Aset.tran)) 
                                | (i.sim %% n.pvalue.cal.min==0)
                                | (i.sim==n.perm.max)) {
                                
                                p.otu.freq.tmp <- 0.5*n.otu.freq # *edited*: omit (1/n.perm.completed) to save time
                                p.otu.tran.tmp <- 0.5*n.otu.tran
                                pmin.otu.omni <- pmin(p.otu.freq.tmp, p.otu.tran.tmp)
                        
                                pnull.otu.freq <- n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                                pnull.otu.tran <- n.perm.completed + 0.5 - apply(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank)
                                pnullmin.otu.omni <- pmin(pnull.otu.freq, pnull.otu.tran)
                                
                                if (length(dim(pnullmin.otu.omni))==3) {
                                    pnullmin.otu.omni <- aperm(pnullmin.otu.omni, c(2,3,1))
                                } else {
                                    pnullmin.otu.omni <- array(pnullmin.otu.omni, c(dim(pnullmin.otu.omni), 1))
                                }
                                
                                n.otu.omni <- apply( (pnullmin.otu.omni < c(pmin.otu.omni) - tol.eq), c(1,2), sum) + 0.5 * apply( (abs(pnullmin.otu.omni - c(pmin.otu.omni)) < tol.eq), c(1,2), sum)
                                
                                if (any(Aset.omni)) {
                                    
                                    AtoB.omni <- Aset.omni & (n.otu.omni >= n.rej.stop)
                                    Aset.omni <- Aset.omni & !AtoB.omni
                                    
                                    p.otu.omni[,otu.smallp][AtoB.omni] <- n.otu.omni[AtoB.omni]*inv.n.perm.completed
                                    p.otu.omni[,otu.smallp][Aset.omni] <- (n.otu.omni[Aset.omni]+1)*inv.n.perm.completed.1
                                    
                                    q.otu.omni <- t(apply(p.otu.omni, 1, fdr.Sandev))
                                    
                                    Aset.omni.meet.criteria <- apply(((q.otu.omni[,otu.smallp,drop=FALSE] < fdr.nominal) & Aset.omni) | (!Aset.omni), 1, all)
                                    n.stable.omni <- ifelse(Aset.omni.meet.criteria, n.stable.omni + Aset.omni.meet.criteria, 0)
                                    Aset.omni.rm.row <- (n.stable.omni >= n.stable.max)
                                    Aset.omni[Aset.omni.rm.row,] = FALSE
                                }
    
                                if (!any(Aset.freq) 
                                    & !any(Aset.tran) 
                                    & !any(Aset.omni)) {
                                    otu.tests.stopped = TRUE 
                                    cat("otu test stopped at permutation", i.sim, "\n")
                                    n.otu.perm.completed = i.sim
                                }
                            }
                        } else {
                            if (!any(Aset.freq)) {
                                otu.tests.stopped = TRUE 
                                cat("otu test stopped at permutation", i.sim, "\n")
                                n.otu.perm.completed = i.sim
                            }
                            q.otu.tran <- NULL
                            q.otu.omni <- NULL
                        }
                    }
                }# if test.otu
          
                if (test.global + test.otu == 2) {
                    if (global.tests.stopped + otu.tests.stopped == 2) break
                    if (otu.tests.stopped == TRUE & i.sim >= n.global.perm.max) break
                } 
          
                if ( test.global + test.otu == 1) {
                    if (!is.null(global.tests.stopped)) {if (global.tests.stopped | i.sim >= n.global.perm.max) break}
                    if (!is.null(otu.tests.stopped)) {if (otu.tests.stopped) break}
                }
          
            }# check if stop early 
        
        }# permutation
        
        n.perm.completed = max(n.global.perm.completed, n.otu.perm.completed, na.rm=TRUE)
        
    }# if (n.perm.max > 0)

    otu.names <- colnames(otu.table)
    colnames(ldm.obs.freq$ve.otu) <- otu.names
    if (!freq.scale.only) colnames(ldm.obs.tran$ve.otu) <- otu.names
    if (!is.null(p.otu.freq)) colnames(p.otu.freq) <- otu.names
    if (!is.null(p.otu.tran)) colnames(p.otu.tran) <- otu.names
    if (!is.null(p.otu.omni)) colnames(p.otu.omni) <- otu.names
    if (!is.null(q.otu.freq)) colnames(q.otu.freq) <- otu.names
    if (!is.null(q.otu.tran)) colnames(q.otu.tran) <- otu.names
    if (!is.null(q.otu.omni)) colnames(q.otu.omni) <- otu.names
    
    # confounder
    
    VE.df.confounders = NULL
    VE.global.freq.confounders = NULL
    VE.otu.freq.confounders = NULL
    VE.global.tran.confounders = NULL
    VE.otu.tran.confounders = NULL
    
    if (adjust.for.confounders) {
        i.conf <- fit.ldm$low[1]:fit.ldm$up[1]
        VE.df.confounders <- length(i.conf)
        
        VE.global.freq.confounders <- sum((fit.ldm$d.freq[i.conf])^2)
        wt <- fit.ldm$d.freq[i.conf] * t(fit.ldm$v.freq[, i.conf]) 
        if (is.vector(wt)) VE.otu.freq.confounders = wt^2
        else               VE.otu.freq.confounders = colSums(wt^2)
        
        if (!freq.scale.only) {
            VE.global.tran.confounders <- sum((fit.ldm$d.tran[i.conf])^2)
            wt <- fit.ldm$d.tran[i.conf] * t(fit.ldm$v.tran[, i.conf]) 
            if (is.vector(wt)) VE.otu.tran.confounders = wt^2
            else               VE.otu.tran.confounders = colSums(wt^2)
        }
    }
    
    # submodels
    
    VE.df.submodels = NULL
    VE.global.freq.submodels = NULL
    VE.otu.freq.submodels = NULL
    VE.global.tran.submodels = NULL
    VE.otu.tran.submodels = NULL
    
    for (k in 1:n.var1) {
        k1 = k + as.numeric(adjust.for.confounders)
        i.m <- fit.ldm$low[k1]:fit.ldm$up[k1]
        VE.df.submodels <- c(VE.df.submodels, length(i.m))
        
        VE.global.freq.submodels <- c(VE.global.freq.submodels, sum((fit.ldm$d.freq[i.m])^2))
        wt <- fit.ldm$d.freq[i.m] * t(fit.ldm$v.freq[, i.m]) 
        if (is.vector(wt)) VE.otu.freq.submodels = rbind(VE.otu.freq.submodels, wt^2)
        else               VE.otu.freq.submodels = rbind(VE.otu.freq.submodels, colSums(wt^2))

        if (!freq.scale.only) {
            VE.global.tran.submodels <- c(VE.global.tran.submodels, sum((fit.ldm$d.tran[i.m])^2))
            wt <- fit.ldm$d.tran[i.m] * t(fit.ldm$v.tran[, i.m]) 
            if (is.vector(wt)) VE.otu.tran.submodels = rbind(VE.otu.tran.submodels, wt^2)
            else               VE.otu.tran.submodels = rbind(VE.otu.tran.submodels, colSums(wt^2))
        }
    }
    
    # residuals
    
    i.all <- fit.ldm$low[1]:fit.ldm$up[n.var]
    
    VE.global.freq.residuals <- fit.ldm$d.freq[-i.all]^2
    VE.global.tran.residuals <- NULL
    if (!freq.scale.only) VE.global.tran.residuals <- fit.ldm$d.tran[-i.all]^2
    
    
    res = list( b=fit.ldm$b,
                dist=d.gower,
                x.freq=x.freq,
                v.freq=fit.ldm$v.freq,
                d.freq=fit.ldm$d.freq,
                x.tran=x.tran,
                v.tran=fit.ldm$v.tran,
                d.tran=fit.ldm$d.tran,
                low=fit.ldm$low,
                up=fit.ldm$up,
                
                VE.global.freq.confounders=VE.global.freq.confounders,
                VE.global.freq.submodels=VE.global.freq.submodels,
                VE.global.freq.residuals=VE.global.freq.residuals,
                VE.otu.freq.confounders=VE.otu.freq.confounders,
                VE.otu.freq.submodels=VE.otu.freq.submodels,
                VE.global.tran.confounders=VE.global.tran.confounders,
                VE.global.tran.submodels=VE.global.tran.submodels,
                VE.global.tran.residuals=VE.global.tran.residuals,
                VE.otu.tran.confounders=VE.otu.tran.confounders,
                VE.otu.tran.submodels=VE.otu.tran.submodels,
                VE.df.confounders=VE.df.confounders,
                VE.df.submodels=VE.df.submodels,
                
                F.global.freq=ldm.obs.freq$ve.global, 
                F.global.tran=ldm.obs.tran$ve.global, 
                F.otu.freq=ldm.obs.freq$ve.otu, 
                F.otu.tran=ldm.obs.tran$ve.otu,
                p.global.freq=p.global.freq, 
                p.global.tran=p.global.tran, 
                p.global.omni=p.global.omni,
                p.otu.freq=p.otu.freq,
                p.otu.tran=p.otu.tran,
                p.otu.omni=p.otu.omni,
                q.otu.freq=q.otu.freq,
                q.otu.tran=q.otu.tran,
                q.otu.omni=q.otu.omni,
                n.perm.completed=n.perm.completed,
                global.tests.stopped=global.tests.stopped,
                otu.tests.stopped=otu.tests.stopped,
                seed=seed)
      
    return(res)
    
} # ldm  



calculate.b.and.resid = function( d.gower, x.freq, x.tran, index, m, adjust.for.confounders) {
    
    n.var = length(index)
    n.otu = ncol(x.freq)
    n.sam = nrow(d.gower)
    ndf.nominal = rep(0, n.var+1)
    
    tol.d = 10^-8
    
    #--------------------------------------------------------------------------
    # construct directions matrix b 
    # from each set of covariates in the list vars
    #--------------------------------------------------------------------------
    
    d.resid = d.gower
    
    for (i in 1:n.var) 
    {
        var = m[,1:index[i]]
        
        svd.var = svd(var)   
        use = (svd.var$d>tol.d)    
        
        hat.matrix = svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        #---------------------
        # calculate direction
        #---------------------
        
        n.dim = dim( hat.matrix)[1]
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) # es: eigen system in Mathematica
        
        use = ( abs(es.model$values)>tol.d )
        ndf.model = sum( use )
        
        b.model = es.model$vectors[, use]
        e.model = es.model$values[use]
        
        hat.matrix.bar = diag(n.dim)  - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # end of calculating direction
        #-----------------------------    
        
        if (i==1) {
            b = b.model
            e = e.model
        } else {   
            b = cbind(b, b.model)
            e = c(e, e.model )
        }
        
        ndf.nominal[i] = ndf.model
        
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE)
    use = which( abs(es.resid$values)>tol.d )
    
    ndf.nominal[n.var+1] = length(use)
    b = cbind(b, es.resid$vectors[, use])
    e = c(e, es.resid$values[use])
    
    
    #---------------------------------------------------------
    # fit LDM to b
    #---------------------------------------------------------
    
    b = matrix(b, nrow=nrow(d.gower)) 
    
    wt.freq = t(b) %*% x.freq    
    d.freq = sqrt(apply(wt.freq^2, 1, sum))
    v.freq = t((1/d.freq)*wt.freq)
    d.tran = NULL
    v.tran = NULL
    if (!is.null(x.tran)) {
        wt.tran = t(b) %*% x.tran    
        d.tran = sqrt(apply(wt.tran^2, 1, sum))
        v.tran = t((1/d.tran)*wt.tran)
    }
    
    #-------------------------------------------------
    # low, up
    #-------------------------------------------------
    
    low = rep(NA, n.var)
    up = rep(NA, n.var)
    
    up.prev = 0
    
    for (k in 1:n.var)
    {
        low[k] = up.prev + 1
        up[k] = up.prev + ndf.nominal[k]
        up.prev = up[k]
    }
    
    #-------------------------------------------------
    # calculate resid, ss.tot
    #-------------------------------------------------
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    ss.tot.freq = matrix( rep(NA, n.otu*n.var1), nrow=n.var1)
    resid.freq = array( NA, dim=c( dim(x.freq), n.var1 ) ) 
    ss.tot.tran = NULL
    resid.tran = NULL
    if (!is.null(x.tran)) {
        ss.tot.tran = matrix( rep(NA, n.otu*n.var1), nrow=n.var1)
        resid.tran = array( NA, dim=c( dim(x.tran), n.var1 ) ) 
    }
    
    for (k in 1:n.var1) {
        
        k1 = ifelse(adjust.for.confounders, k+1, k)
        use = setdiff( 1:up[n.var], low[k1]:up[k1] )
        
        resid.freq[,,k] = x.freq - b[,use,drop=FALSE] %*% wt.freq[use,,drop=FALSE]
        ss.tot.freq[k,] = colSums( resid.freq[,,k]^2 )
        if (!is.null(x.tran)) {
            resid.tran[,,k] = x.tran - b[,use,drop=FALSE] %*% wt.tran[use,,drop=FALSE]
            ss.tot.tran[k,] = colSums( resid.tran[,,k]^2 )
        } 
    }

    res = list( b=b,
                d.freq=d.freq,
                v.freq=v.freq,
                d.tran=d.tran,
                v.tran=v.tran,
                low=low,
                up=up,
                resid.freq=resid.freq,
                resid.tran=resid.tran,
                ss.tot.freq=ss.tot.freq,
                ss.tot.tran=ss.tot.tran,
                ndf=ndf.nominal)
    
    return(res)
    
} # calculate.b.and.resid


ldm.stat = function(b, low, up, resid, ss.tot, adjust.for.confounders) {
    
    #---------------------------------------------
    #  calculate FL statistics for each model
    #---------------------------------------------
    
    n.var = length(low)
    n.otu = dim(resid[,,1,drop=FALSE])[2]
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    ve.otu = matrix(rep(NA, n.otu*n.var1), nrow=n.var1 )
    ve.global = rep(NA, n.var1)

    for (k in 1:n.var1) {
        
        k1 = k + as.numeric(adjust.for.confounders)
        use = low[k1]:up[k1]
        
        wt = t( b[, use] ) %*% resid[,,k]
        ve.otu.k = colSums( wt^2 )
        
        if (n.var==1) {
            sigma.k = ss.tot[k,] - ve.otu.k
        } else {
            use = 1:up[n.var]
            wt.cum = t( b[, use] ) %*% resid[,,k]
            sigma.k = ss.tot[k,] - colSums( wt.cum^2 )
        }
        
        ve.otu[k,] = ve.otu.k/sigma.k
        ve.global[k] = sum( ve.otu.k )/sum( sigma.k )
        
    }
    
    out = list( ve.otu=ve.otu, 
                ve.global=ve.global)   
    
    return(out)
    
} # ldm.stat


#' PERMANOVA test of association
#' 
#' This function performs the PERMANOVA test that can allow adjustment of
#' confounders and control of clustered data. As in \code{ldm}, 
#' \code{permanovaFL} allows multiple sets of covariates to be tested, 
#' in the way that the sets are entered sequentially and the variance 
#' explained by each set is that part that remains after the previous 
#' sets have been fit. 
#' 
#' @param formula a symbolic description of the model to be fitted in the form
#'   of \code{data.matrix ~ sets of covariates} or \code{data.matrix |
#'   confounders ~ sets of covariates}. The details of model specification are
#'   given in "Details" of \code{ldm}. Additionally, in \code{permanovaFL}, the \code{data.matrix}
#'   can be either an OTU table or a distance matrix. If it is an OTU table,
#'   the distance matrix will be calculated internally using the OTU table, \code{tree} (if required), and 
#'   \code{dist.method}. If \code{data.matrix} is a distance
#'   matrix (having class \code{dist} or \code{matrix}), it can be squared and//or centered by
#'   specifying \code{square.dist} and \code{center.dist} (described below).  Distance matrices are distinguished
#'   from OTU tables by checking for symmetry of \code{as.matrix(data.matrix)}.
#' @param data an optional data frame, list or environment (or object coercible 
#' to a dataframe) containing the covariates of interest and confounding covariates. 
#' If not found in \code{data}, the covariates are taken from environment(formula), 
#' typically the environment from which \code{permanovaFL} is called. The default is .GlobalEnv.
#' @param tree a phylogenetic tree. Only used for calculating a 
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of 
#'   the requested distance does not involve a phylogenetic tree, or if a 
#'   distance matrix is directly imported through \code{formula}.
#' @param dist.method method for calculating the distance measure, partial
#' match to all methods supported by \code{vegdist} in the \code{vegan} package
#'  (i.e., "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", 
#'  "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis")
#'   as well as "hellinger" and "wt-unifrac". 
#'   Not used if a distance matrix is specified in \code{formula}.
#'   The default is "bray". 
#'   For more details, see the \code{dist.method} argument in the \code{ldm} function.
#' @param cluster.id cluster identifiers. The default is value of NULL should be used if the observations are 
#' not in clusters (i.e., independent).
#' @param perm.within.type a character string that takes values "free", "none", "series", or "grid".  
#'   The default is "free" (for random permutations).
#' @param perm.within.nrow a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See the documentation for the R package \code{permute} for further details.
#' @param perm.within.ncol a positive integer, only used if perm.within.type="grid". 
#'   The default is 0.  See the documentation for the R package \code{permute} for further details.
#' @param perm.between.type a character string that takes values "free", "none", or "series".  
#'   The default is "none".
#' @param strata a factor variable (or, character variable converted into a factor) to define strata (groups), within which to constrain permutations. 
#'   The default is NULL.
#' @param how a permutation control list, for users who want to specify their permutation control list using the \code{how} function 
#'   from the \code{permute} R package.  The default is NULL.
#' @param n.perm.max the maximum number of permutations.
#'   The default is 5000.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. The 
#'   default is 20.
#' @param seed an integer seed for the random number generator in the 
#'   permutation procedure. The default is NULL; with the default value, an integer seed will be 
#'   generated internally and randomly. In either case, the integer seed will be 
#'   stored in the output object.
#' @param square.dist a logical variable indicating whether to square the
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale the OTU table.
#'   For count data, this corresponds to dividing by the library size to give
#'   relative frequencies.  The default is TRUE.
#' @return  a list consisting of 
#' \item{F.statistics}{F statistics for the global test of each set of covariates}
#' \item{p.permanova}{p-values for the global test
#'   of each set of covariates} 
#' \item{n.perm.completed}{number of permutations completed}
#' \item{permanova.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all global tests}
#' \item{seed}{a single-value integer seed that is user supplied or internally generated}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @examples
#' res.permanova <- permanovaFL(formula=throat.otu.tab | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
#'                             data=throat.meta, dist.method="bray", seed=123)


permanovaFL = function(formula, data=.GlobalEnv, tree=NULL, dist.method="bray",
                       cluster.id=NULL, strata=NULL, how=NULL,
                       perm.within.type="free", perm.between.type="none",
                       perm.within.ncol=0, perm.within.nrow=0,
                       n.perm.max=5000, n.rej.stop=20, seed=NULL,
                       square.dist=TRUE, center.dist=TRUE, scale.otu.table=TRUE) {  
    
    # dm = form.call( object=formula, data=data, cluster.id=cluster.id )
    
    #------------------------
    # form.call
    #------------------------
    
    object=formula
    #
    #   extract cluster.id from dataframe
    #
    cl=match.call()
    mf=match.call(expand.dots=FALSE)
    m=match( x='cluster.id', table=names(mf) )
    mf.string=as.character( mf[c(1L,m)] )
    cluster.name=mf.string[2]
    if (cluster.name=='NULL') {
        cluster.id=NULL
    } else {   
        loc.dollar=tail( gregexpr('\\$', cluster.name)[[1]] , n=1 )
        if (loc.dollar<0)  {
            cluster.id=getElement(data,cluster.name)
            if( is.null(cluster.id) ) cluster.id=get(cluster.name)
        } else {   
            df.name=get( substr(cluster.name, start=1, stop=loc.dollar-1) )
            var.name=substr(cluster.name, start=loc.dollar+1, stop=nchar(cluster.name))            
            cluster.id= getElement(df.name,var.name) 
        }
    }
    #        
    #   extract model from formula    
    #    
    obj=toString(object)
    obj=gsub('\\s','',obj)
    prefix=' ~ + 0 + '
    loc.comma=gregexpr(',',obj)[[1]]
    start.terms=loc.comma[2]
    terms=substr(obj,start=start.terms+1, stop=nchar(obj))
    #
    #   find n.obs and full set of rownames
    #   
    if (class(data)=='data.frame') {
        row.names=rownames(data)
        n.obs=length(row.names)
    } else {   
        df=model.frame( paste('~',terms) , na.action=na.pass )
        row.names=rownames(df)
        n.obs=length(row.names)
    }
    #
    #   check for missing values in cluster.id
    #        
    
    if (is.null(cluster.id)) {
        use.rows=row.names
    } else {   
        use=!is.na(cluster.id)
        use.rows=row.names[use]
    }
    #
    #   check for and extract confounders
    #
    model=list()
    j=1
    loc.bar=regexpr('\\|',obj)[1]
    loc.minus=regexpr('-',obj)[1]
    loc.delim=max( loc.bar, loc.minus)
    if (loc.delim>0) {
        end.confound=loc.comma[2]
        c=substr(obj,start=loc.delim+1, stop=end.confound-1)
        conf=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        model[[j]]=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
        #       use.rows=intersect( use.rows, rownames(conf) )
        use.rows=rownames(model[[1]]) 
        j=j+1
    } else {
        conf=NULL
    }     
    #
    #   extract model terms
    #
    #   j=1
    continue=TRUE
    while (continue) {
        if (substr(terms,1,1)=='(') {
            stop=regexpr(')\\+',terms)[1]
        } else {
            stop=regexpr('\\+',terms)[1] - 1
        }          
        
        if (stop<=0) stop=nchar(terms) 
        m=substr(terms, start=1, stop=stop)
        model[[j]]=model.matrix( as.formula( paste(prefix,m) ) , data=data)
        use.rows=intersect( use.rows, rownames(model[[j]]) )
        #        if (j==1) {
        #            use.rows=rownames(model[[1]])
        #            }
        #        else {
        #            use.rows=intersect( use.rows, rownames(model[[j]]) )
        #            }         
        if (stop+2<=nchar(terms)) {
            terms=substr(terms, start=stop+2, stop=nchar(terms))
            j=j+1
        } else {
            continue=FALSE
        }             
    }   
    n.model=j    
    #
    #  extract OTU table
    #      
    if (is.null(conf)) loc.delim=loc.comma[2]
    otu.name=substr(obj, start=loc.comma[1]+1, stop=loc.delim-1)
    #   loc.dollar=regexpr('\\$', otu.name)[1]
    loc.dollar=tail( gregexpr('\\$', otu.name)[[1]] , n=1 )
    if (loc.dollar<0)  {
        if (class(data)=='data.frame') {
            otu.table=getElement(data, otu.name)
            if (is.null(otu.table)) otu.table= get(otu.name) 
            otu.table=as.matrix(otu.table)
        } else {
            otu.table=as.matrix( get(otu.name) )
        }
    } else {
        df.name=get( substr(otu.name, start=1, stop=loc.dollar-1) )
        var.name=substr(otu.name, start=loc.dollar+1, stop=nchar(otu.name))
        otu.table=as.matrix( getElement(df.name,var.name) )
    }        
    #    if (is.null(otu.table)) otu.table=as.matrix( getElement(.GlobalEnv,otu.name) )
    if ( nrow(otu.table) != n.obs ) {
        if (ncol(otu.table)==n.obs ) {
            otu.table=t(otu.table)
        } else {   
            print('warning: OTU table and covariates have different number of observations')
            return
        }
    }   
    #
    #   remove rows having NA 
    #    
    for (j in 1:n.model) {
        keep =  rownames( model[[j]] ) %in% use.rows
        model[[j]]=model[[j]][keep,,drop=FALSE]
    }
    if (!is.null(conf)) {
        keep =  rownames(conf) %in% use.rows 
        conf=conf[keep,,drop=FALSE]
    }
    keep=row.names %in% use.rows    
    otu.table=otu.table[keep,,drop=FALSE]    
    if (!is.null(cluster.id)) cluster.id=cluster.id[keep]
    
    otu.or.dist <- as.matrix(otu.table)
    
    otu.table <- NULL
    dist <- NULL
    if (isSymmetric(otu.or.dist)) {
        dist <- otu.or.dist
        if (dim(model[[1]])[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between covariates and the distance matrix' )
    } else {
        otu.table <- otu.or.dist
        if (dim(model[[1]])[1] != dim(otu.table)[1]) 
            otu.table <- t(otu.table)
        if (dim(model[[1]])[1] != dim(otu.table)[1]) stop( 'numbers of observations mismatch between covariates and the OTU table' )
    }
    
    #------------------------
    # dist matrix
    #------------------------
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.method=dist.method, otu.table=otu.table, tree=tree, scale.otu.table=scale.otu.table)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)
    
    #------------------------
    # setup model
    #------------------------
    
    n.var = length(model)
    n.obs = dim(model[[1]])[1]
    adjust.for.confounders = !is.null(conf)
    
    center.vars=TRUE
    
    index = rep(0, n.var)
    
    for (i in 1:n.var) {
        m.i = model[[i]]
        if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
        
        if (i==1) {
            m = m.i
            index[i] = dim(m.i)[2] 
        } else {
            m = cbind(m, m.i)   
            index[i] = index[i-1] + dim(m.i)[2]    
        }
        
    }    
    
    #------------------------
    # setup permutation
    #------------------------
    
    if (class(how)=='how') {
        CTRL=how                   # user-provided how list
    }
    else {
        if (is.null(cluster.id)) {
            if (is.null(perm.within.type) & is.null(perm.between.type)) {
                # default when no unclustered data has no type specified is 'free'
                perm.within.type='free'    
            }
            if (is.null(strata)) {
                # setup for unclustered permutation
                CTRL = how( within=Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))  
            }
            else {
                # setup for unclustered, stratified permutation
                strata=as.factor(strata)
                CTRL = how( blocks=strata, within=Within(type=perm.within.type, 
                                                         nrow=perm.within.nrow, 
                                                         ncol=perm.within.ncol))  
            }    
        }
        else {        
            cluster.id=as.factor(cluster.id)
            if (is.null(strata)) {            
                #  clustered but unstratified data
                CTRL = how( plots=Plots(cluster.id, type=perm.between.type ), 
                            within=Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            }
            else {
                #   clustered and stratified data
                strata=as.factor(strata)             
                CTRL = how( blocks=strata, 
                            plots=Plots(cluster.id, type=perm.between.type ), 
                            within=Within(type=perm.within.type, 
                                          nrow=perm.within.nrow, 
                                          ncol=perm.within.ncol))
            }
        }
    }  
    
    #---------------------
    # observed statistic
    #---------------------
    
    fit.res = fit.permanova( d.gower=d.gower, index=index, m=m, adjust.for.confounders=adjust.for.confounders) 
    
    permanova.obs = permanova.stat(b=fit.res$b, low=fit.res$low, up=fit.res$up, resid.dist=fit.res$resid.dist, ndf=fit.res$ndf, adjust.for.confounders=adjust.for.confounders)
    
    p.permanova = NULL
    n.perm.completed = NULL
    permanova.stopped = NULL
    
    #---------------------
    # permutation
    #---------------------
    
    if (n.perm.max > 0) {
        
        tol.eq = 10^-16
        n.perm.block = 1000
        
        n.permanova = 0
        n.perm.completed = 0
        
        if (is.null(seed)) {
            seed = sample(1:2^15, 1)
        }
        set.seed(seed)
        
        for (i.sim in 1:n.perm.max) {
            
            i.sim.r = i.sim%%n.perm.block
            if (i.sim.r==0) i.sim.r = n.perm.block
            
            if (i.sim.r==1) {
                cat("permutations:", i.sim, "\n")
                perm = shuffleSet(n.obs, n.perm.block, CTRL)
            }
            
            n.perm.completed = n.perm.completed + 1
            
            # perform permutations                   
            
            b.perm = fit.res$b[perm[i.sim.r,], ]   
            
            permanova.perm = permanova.stat(b=b.perm, low=fit.res$low, up=fit.res$up, resid.dist=fit.res$resid.dist, ndf=fit.res$ndf, adjust.for.confounders=adjust.for.confounders)
            
            n.permanova <- n.permanova + (permanova.perm$permanova > permanova.obs$permanova + tol.eq) + (permanova.perm$permanova > permanova.obs$permanova - tol.eq)
            
            if (all(n.permanova >= n.rej.stop*2)) break
            
        }# permutation
        
        permanova.stopped=all(n.permanova >= n.rej.stop*2)
        
        p.permanova <- ifelse((n.permanova >= n.rej.stop*2), 0.5*n.permanova*(1/n.perm.completed), (0.5*n.permanova+1)*(1/(n.perm.completed+1)))
        
    }# if (n.perm.max > 0)
    
    res = list( F.statistics=permanova.obs$permanova,
                p.permanova=p.permanova, 
                n.perm.completed=n.perm.completed, 
                permanova.stopped=permanova.stopped,
                seed=seed)
    return(res)
    
}# permanovaFL


fit.permanova = function( d.gower, index, m, adjust.for.confounders) {
    
    n.var = length(index)
    n.otu = ncol(d.gower)
    n.sam = nrow(d.gower)
    ndf.nominal = rep(0, n.var+1)
    
    tol.d = 10^-8
    
    #--------------------------------------------------------------------------
    # construct directions matrix b 
    # from each set of covariates in the list vars
    #--------------------------------------------------------------------------
    
    d.resid = d.gower
    
    for (i in 1:n.var) 
    {
        
        var = m[,1:index[i]]
        
        svd.var = svd(var)   
        use = (svd.var$d>tol.d)    
        
        hat.matrix = svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        #---------------------
        # calculate direction
        #---------------------
        
        n.dim = dim( hat.matrix)[1]
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) # es: eigen system in Mathematica
        
        use = ( abs(es.model$values)>tol.d )
        ndf.model = sum( use )
        
        b.model = es.model$vectors[, use]
        e.model = es.model$values[use]
        
        hat.matrix.bar = diag(n.dim)  - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # end of calculating direction
        #-----------------------------    
        
        if (i==1) {
            b = b.model
            e = e.model
        } else {   
            b = cbind(b, b.model)
            e = c(e, e.model )
        }
        
        ndf.nominal[i] = ndf.model
        
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE)
    use = which( abs(es.resid$values)>tol.d )
    
    ndf.nominal[n.var+1] = length(use)
    b = cbind(b, es.resid$vectors[, use])
    e = c(e, es.resid$values[use])
    
    #-------------------------------------------------
    # low, up
    #-------------------------------------------------
    
    low = rep(NA, n.var)
    up = rep(NA, n.var)
    
    up.prev = 0
    
    for (k in 1:n.var)
    {
        low[k] = up.prev + 1
        up[k] = up.prev + ndf.nominal[k]
        up.prev = up[k]
    }
    
    #---------------------
    # permanova: resid.dist
    #---------------------
    
    if (n.var==1) {
        resid.dist = array( NA, dim=c( dim(d.gower), n.var ) ) 
        resid.dist[,,1] = d.gower
    }
    else {
        n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
        
        resid.dist = array( NA, dim=c( dim(d.gower), n.var1 ) ) 
        
        for (k in 1:n.var1) {
            k1 = ifelse(adjust.for.confounders, k+1, k)
            use = setdiff( 1:up[n.var], low[k1]:up[k1] )
            
            hat.matrix = b[,use,drop=FALSE] %*% t( b[,use,drop=FALSE] )
            hat.matrix.bar = diag(n.sam) - hat.matrix
            resid.dist[,,k] = hat.matrix.bar %*% d.gower # the other hat.matrix.bar is not needed
        }
    }
    
    
    res = list( b=b,
                low=low,
                up=up,
                resid.dist=resid.dist,
                ndf = ndf.nominal)
    return(res)
    
} # fit.permanova


permanova.stat = function(b, low, up, resid.dist, ndf, adjust.for.confounders) {
    
    #---------------------------------------------
    #  calculate FL statistics for each model
    #---------------------------------------------
    
    n.var = length(low)
    n.sam = nrow(b)
    
    n.var1 = ifelse(adjust.for.confounders, n.var-1, n.var)
    
    permanova = rep(NA, n.var1)
    
    if (n.var==1) {
        
        use = 1:up[1]    
        H = b[, use] %*% t( b[, use] )
        I_H = diag(n.sam) - H
        
        permanova[1] = sum( diag(H %*% resid.dist[,,1]) )/sum( diag(I_H %*% resid.dist[,,1]) ) # the other H and I_H is not needed
        
    }
    else {
        
        use = 1:up[n.var]
        Hcum = b[, use] %*% t( b[, use] )
        I_Hcum = diag(n.sam) - Hcum
        
        for (k in 1:n.var1) {
            
            k1 = k + as.numeric(adjust.for.confounders)
            
            use = low[k1]:up[k1]
            Hk = b[, use] %*% t( b[, use] )
            
            permanova[k] = sum( diag(Hk %*% resid.dist[,,k]) )/sum( diag(I_Hcum %*% resid.dist[,,k]) ) # the other Hk and I_Hcum is not needed
        }
    }   
    
    var = ifelse(adjust.for.confounders, 2:n.var, 1:n.var)
    out = list( permanova=permanova * ndf[n.var+1] / ndf[var]) 
    
    return(out)
    
} # permanova.stat


