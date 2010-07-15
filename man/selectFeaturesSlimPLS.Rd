\name{selectFeaturesSlimPLS}
\alias{selectFeaturesSlimPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selects and constructs a set of features using the SlimPLS algorithm methods from a given expression matrix}
\description{
  The function constructs a small number of features, each of them is a linear combination of basic features
  from the given expression matrix with respect to a given labeling of the samples. The new features extracted
  can later be used for reducing the dimension of an expression matrix in order to learn a classification model more efficiently.
  While learning, either the extracted features (linear combinations, aka components) which are output of this function,
  or merely the list of features that are members of any component can be used.
}
\usage{
selectFeaturesSlimPLS(exp_mat, num_class_a=0, num_class_b=0, class_a="",
		      class_b="", num_features, component_size, p_value_threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{exp_mat}{ \code{exp_mat} is an expression matrix of type \code{expMat}, usually created by reading an
  expression matrix from a file using \code{readExpMat}. The matrix is supposed to have all samples belonging to class
  1 grouped together in the first columns, following all samples from class 2 grouped together as well. }
  \item{num_class_a}{ \code{num_class_a} is the number of samples belonging to class 1 }
  \item{num_class_b}{ \code{num_class_b} is the number of samples belonging to class 2 }
  \item{class_a}{ \code{class_a} is the class label of class a. May be used if the expression matrix has labels in its second row.}
  \item{class_b}{ \code{class_b} is the class label of class b. May be used if the expression matrix has labels in its second row.}
  \item{num_features}{ The maximal number of basic features that will be selected }
  \item{component_size}{ The number of features in each constructed component. }
  \item{p_value_threshold}{ An optional p-value threshold on the features selected. If 0, this parameter is ignored. If not,
  only features below the given p-value will be selected. }
}
\details{
 The function selects up to \code{num_features} basic features, grouped together as a linear combination of them 
 in components of the given \code{component_size}. If \code{p_value_threshold} is not 0, it selects only features
 with a p-value lower than the given threshold. See references for method details.
 The selection is supervised, based on two classes, so the user must supply a class parameter for every sample. This
 is done by providing a sorted matrix, where all the samples from class a precede all the samples from class b.
 In addition to the sorted matrix, the user must provide the number of samples from class a and b.
 Alternatively the user may provide an unsorted matrix, with labels in its second row, denoting the class of each 
 sample. In this case the user must also provide two labels to the selection function - one for class a and one
 for class b. 
 }
\value{
  The function returns an object of type \code{featureSet}, which holds the selected features, and can later be used for
  learning a classification model by \code{trainClassifier}. The object holds both the features, and their weighted grouping into
  components. Each component is a linear combination of basic features.
  }
\references{
1. Gutkin, M, MSc thesis, April 2008, Feature selection methods for classification of gene expression profiles,
Tel Aviv Univetsity, Tel Aviv. 

2. Gutkin, M., Shamir, R., Dror, G., 2009 SlimPLS: A Method for Feature Selection in Gene Expression-Based Disease 
Classification. PLoS ONE 4(7): e6416. doi:10.1371/journal.pone.0006416 
}
\author{ Michael Gutkin }
\note{ }
\seealso{ \code{\link{selectFeatures}}, \code{\link{trainClassifier}},\code{\link{getClassification}},\code{\link{readExpMat} }}
\examples{
# reads an expression matrix with no class labels into exp_mat
\dontrun{
exp_mat <- readExpMat("golub_leukemia_data_training.csv", FALSE)}
	
# selects a set of features into the features variable. The matrix we read is
# sorted by classes, and 
# it is known to have 41 samples of class a and 20 samples of class b. The selection is
# done using the SlimPLS method. Up to two components with 25 features in each component
# will be selected.

\dontrun{
features <- selectFeaturesSlimPLS(exp_mat, num_class_a=41, num_class_b=20, class_a="",
				  class_b="", num_features=50, component_size=25,
				  p_value_threshold=0)}
 	
# reads an expression matrix with class labels into exp_mat2
\dontrun{
exp_mat2 <- readExpMat("golub_leukemia_data_with_classes_training.csv", FALSE)}

# selects a set of features into the features2 variable. The matrix we read has a class
# label for each sample in its second row. Labels are either "AML" or "ALL". 
# Selection is done using the SlimPLS method. Up to two components with 25 features in
# each component will be selected.
\dontrun{
features2 <- selectFeaturesSlimPLS(exp_mat, class_a="AML", class_b="ALL",
				   num_features=50, component_size=25,
				   p_value_threshold=0)}
	
# the found list of features:
\dontrun{
features2@features}
	
# the found components (weighted matrix for each feature and each component, of size
# 2X50 - 25 non-zero feature weights for each feature:
\dontrun{
features2@w_mat}
	
}

