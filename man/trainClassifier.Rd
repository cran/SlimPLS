\name{trainClassifier}
\alias{trainClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates a classification model using a labeled expression matrix and a
 set of features}
\description{
   This function creates a binary classification model using various supervised learning algorithms and
   the set of features selected or extracted by \code{selectFeatures} \code{selectFeaturesSlimPLS} which 
   output a \code{featureSet} object. The model can be later used to classify new instances using 
   \code{get_prdeiction} function. The function wraps svm and random forest classifiers supplied by external packages
   \code{e1701} and \code{randomForest}.} 
\usage{
trainClassifier(learn_method, exp_mat, num_class_a=0, num_class_b=0, class_a="",
		class_b="", feature_set, use_components_as_features = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{learn_method}{ \code{learn_method} can receive one of three strings denoting differen classification algorithms:
  \code{"SVM_LINEAR"} - a linear kerenl svm, \code{"SVM_RADIAL"} - a radial kerenl svm, 
  \code{"RANDOM_FOREST"} - a random forest based collection of decision trees. The svm variants 
  chooses the cost parameter by leave one out cross validation. The random forest size bulit is at most 1500 trees.
  }
  \item{exp_mat}{ \code{exp_mat} is an expression matrix of type \code{expMat}, usually created by reading an
  expression matrix from a file using \code{readExpMat}. The matrix is supposed to have all samples belonging to class
  1 grouped together in the first columns, following all samples from class 2 grouped together as well. }
  \item{num_class_a}{ \code{num_class_a} is the number of samples belonging to class 1}
  \item{num_class_b}{ \code{num_class_b} is the number of samples belonging to class 2 }
  \item{class_a}{ \code{class_a} is the class label of class a. May be used if the expression matrix has labels in its second row.}
  \item{class_b}{ \code{class_b} is the class label of class b. May be used if the expression matrix has labels in its second row.}
  \item{feature_set}{ \code{feature_set} is an object of type \code{featureSet}, created by
  invoking either \code{selectFeatures} or \code{selectFeaturesSlimPLS}. }
  \item{use_components_as_features}{ The \code{use_components_as_features} parameter is only applicable for SlimPLS features.
  If \code{FALSE} features are used as is. If \code{TRUE} new features are constructed from the basic features. New features
  are a set of linear combinations of the basic features.  }
}
\details{
  The selection method by \code{feature_set} is created is transparent to the \code{trainClassifier} function.
  The \code{feature_set} already holds both information about the set of selected features, and for SlimPLS based 
  features it also holds their weights, and information required for constructing new features which are linear 
  combinations (aka compnents) of the selected features.
  Learning is based on external packages (\code{e1071} and \code{randomForest}.
  The function simply reduces the dimension of the given expression matrix using the selected features supplied,
  and invokes a learning algorithm on the reduced matrix.
  The learning is supervised, based on two classes, so the user must supply a class parameter for every sample. This
  is done by providing a sorted matrix, where all the samples from class a precede all the samples from class b.
  In addition to the sorted matrix, the user must provide the number of samples from class a and b.
  Alternatively the user may provide an unsorted matrix, with labels in its second row, denoting the class of each 
  sample. In this case the user must also provide two labels to the selection function - one for class a and one
  for class b. 
  }
\value{
	The function returns a \code{classificationModel} object, which holds the actual model
	created by the external packages, and the feature set used for building the model. The object
	can later be used to get prediction from new instances.
}
\references{1. Boser, B.E., Guyon, I.M., Vapnik, V.N., A training algorithm for optimal margin classifiers. 
In COLT '92: Proceedings of the fifth annual workshop on Computational Learning Theory, pages 144-152, New York,
NY, USA, 1992. ACM.

2. Liaw, A. and Wiener M., Classification and Regression by randomForest, R News 2(3), 18-22, 2002.
}
\author{ Michael Gutkin, Ofer Lavi }
\seealso{ \code{\link{selectFeatures}}, \code{\link{selectFeaturesSlimPLS}}, \code{\link{getClassification}}, \code{\link{readExpMat}}
}
\examples{

# reads an expression matrix with class labels into exp_mat2
\dontrun{
exp_mat2 <- readExpMat("golub_leukemia_data_with_classes_training.csv", TRUE)}

# selects a set of features into the features2 variable. The matrix we read has a class
# label for each sample in its second row. Labels are either "AML" or "ALL".
# Selection is done using the SlimPLS method. Up to two components with 25 features in
# each component will be selected.
\dontrun{
features2 <- selectFeaturesSlimPLS(exp_mat, class_a="AML", class_b="ALL",
				   num_features=50, component_size=25,
				   p_value_threshold=0)}
# train a SVM classifier with linear kernel on the expression matrix using the
# individual feature that are part of the components selected earlier, and are now in
# features2.
\dontrun{
model_t <- trainClassifier("SVM_LINEAR", exp_mat2, 0,0, "AML", "ALL", features2,
			   FALSE )}
	
}
