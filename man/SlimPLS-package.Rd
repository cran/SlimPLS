\name{SlimPLS-package}
\alias{SlimPLS-package}
\docType{package}
\title{
slimPLS - Partial Least Squares based multivariate feature selection.
}
\description{
This package introduces a multivariate feature selection algorithm developed in the Computational
Genomics lab of Tel-Aviv University primarily for dimension reduction of gene expression data for the
purpose of classification. The package defines methods for reading an expression data from a file, 
and extracting features using SlimPLS from such a matrix. For the purposes of testing it also provides
a wrapper to a number of classification packages, allowing for learning and predicting of samples'
classes learned after the SlimPLS dimension reduction. For comparison, it also provides an interface
for a number of univariate feature selection methods based on scoring.
}
\details{
The usual flow of using SlimPLS starts at reading an expression matrix from a standard
tab delimited file into an \code{ExpMat} object, by calling:



\code{my_exp_mat <- readExpMat("myFile.txt")}



Next, a set of feature is selected from this expression matrix by calling either 



\code{my_feature_set <- selectFeatures(expMat,...)}



 or

 
 
\code{my_feature_set <- selectFeaturesSlimPLS(expMat,...)}.
Then the features and the matrix are used to learn a classification model using 



\code{my_model <- trainClassifier(...)}. 



Last, the model learned can be used to predict classes of another (or the same) expression matrix, using



\code{getClassification(...)}.


The method is used for binary classification (i.e. two classes only), and is supervised (i.e. one needs to
supply a labeling of samples as belonging to one of two classes)  
}
\author{
The method was developed by Michael Gutkin, Prof. Ron Shamir and Prof. Gideon Dror, and programming was
done by Michael Gutkin. It is maintained by Ofer Lavi.
}
\references{
1. Gutkin, M., MSc thesis, April 2008, Feature selection methods for classification of gene expression profiles,
Tel Aviv Univetsity, Tel Aviv. 

2. Gutkin, M., Shamir, R., Dror, G., 2009 SlimPLS: A Method for Feature Selection in Gene Expression-Based Disease 
Classification. PLoS ONE 4(7): e6416. doi:10.1371/journal.pone.0006416 
}

\keyword{ package }
\seealso{
}
\examples{

}
