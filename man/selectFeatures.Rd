\name{selectFeatures}
\alias{selectFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Selects the best features using various scoring methods from a given expression matrix}
\description{
  The function selects the top scoring features according to the given scoring method in a supervised manner.
  Every scoring method computes a score for each feature in the expression matrix, with respect to a given
  labeling of the samples. The selected features can later be used for reducing the dimension of an expression matrix
  in order to learn a classification model more efficiently. }

\usage{
selectFeatures(select_method, exp_mat, num_class_a=0, num_class_b=0, class_a="",
	       class_b="", num_features)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{select_method}{ \code{select_method} can be one of the following strings, denoting the scoring function for a feature:
  \code{"CORREL"} - Pearson correlation score of the expression levels and the labeling vector.
  \code{"TTEST"} - Two samples Student's t-test score for the set of expression levels from class 1 and the set of expression levels from class 2. 
  \code{"FC"} - Fold change score between the average of values of samples from class 1 and samples from class 2
  \code{"GOLUB"} - Golub criterion score.
  \code{"MI"} - Mutual information between the expression levels and the labeling vector. Can be used only if bioDist is installed.}
  \item{exp_mat}{ \code{exp_mat} is an expression matrix of type \code{expMat}, usually created by reading an
  expression matrix from a file using \code{readExpMat}. The matrix is supposed to have all samples belonging to class
  1 grouped together in the first columns, following all samples from class 2 grouped together as well. }
  \item{num_class_a}{ \code{num_class_a} is the number of samples belonging to class 1}
  \item{num_class_b}{ \code{num_class_b} is the number of samples belonging to class 2 }
  \item{class_a}{ \code{class_a} is the class label of class a. May be used if the expression matrix has labels in its second row.}
  \item{class_b}{ \code{class_b} is the class label of class b. May be used if the expression matrix has labels in its second row.}
  \item{num_features}{ \code{num_features} is the total number of features to be selected. }
}
\details{
 The function selects \code{num_features} features that achieve the best scoring using the selected scoring method.
 The selection is supervised, based on two classes, so the user must supply a class parameter for every sample. This
 is done by providing a sorted matrix, where all the samples from class a precede all the samples from class b.
 In addition to the sorted matrix, the user must provide the number of samples from class a and b.
 Alternatively the user may provide an unsorted matrix, with labels in its second row, denoting the class of each 
 sample. In this case the user must also provide two labels to the selection function - one for class a and one
 for class b. 
 For the mutual information selection filter, bioDist needs to be installed first.
 See http://www.bioconductor.org/packages/release/bioc/html/bioDist.html
}
\value{
  returns a \code{featureSet} object, which holds the selected features, and can later be used for
  learning a classification model by \code{trainClassifier}
}
\references{ 1. Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP,
Coller H, Loh ML, Downing JR, Caligiuri MA et al: Molecular classification
of cancer: class discovery and class prediction by gene expression
monitoring. Science 1999, 286(5439):531-537.

2. Hamming RW: Coding and Information Theory: Prentice-Hall Inc.; 1980.

3. Everitt BS, Hothorn T: A Handbook of Statistical Analyses Using R:
Chapman & Hall/CRC Taylor & Francis Group; 2006.}
\author{ Michael gutkin, Ofer Lavi }
\note{ }
\seealso{ \code{\link{selectFeatures}}, \code{\link{selectFeaturesSlimPLS}}, \code{\link{getClassification}}, \code{\link{readExpMat}} }
\examples{
# reads an expression matrix with no class labels into exp_mat
	
\dontrun{
exp_mat <- readExpMat("golub_leukemia_data_training.csv", FALSE)}
	
# selects a set of features into the features variable. The matrix we read is sorted by
# classes, and  it is known to have 41 samples of class a and 20 samples of class b.
# The selection is done using best 100 two-samples t-test scores.
	
\dontrun{
features <- selectFeatures("TTEST", exp_mat, num_class_a=41, num_class_b=20,
			    num_features=100)}

# reads an expression matrix with class labels into exp_mat2

\dontrun{
exp_mat2 <- readExpMat("golub_leukemia_data_with_classes_training.csv", TRUE)}

# selects a set of features into the features2 variable. The matrix we read has a class
# label for each sample in its second row. Labels are either "AML" or "ALL". 
# Selection is done using the fold change score, selecting the 200 features with the 
# largest fold change between the "AML" samples and the "ALL" samples.
\dontrun{
features2 <- selectFeatures("FC", exp_mat, num_class_a=41, num_class_b=20,
			     num_features=200)}
	
}


