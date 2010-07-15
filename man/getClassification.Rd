\name{getClassification}
\alias{getClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predicts classes of a samples in a given expression matrix }
\description{
  Predicts classes of a samples in a given expression matrix, using a given classification model. 
  The function runs a classification prediction algorithm on the given expression matrix using a model learned
  earlier by \code{trainClassifier}. The model holds both a set of features and a specific model
  suited for a specific algorithm. 
}
\usage{
getClassification(exp_mat, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp_mat}{ \code{exp_mat} is an expression matrix of type \code{expMat}, usually created by reading an
  expression matrix from a file using \code{readExpMat}. The matrix is supposed to have samples in its columns, and features
  in its rows.
  }
  \item{model}{ \code{model} is a \code{classificationModel} object, output by \code{trainClassifier}. It holds
  both the features selected earlier for dimension reduction of the expression matrix, and a specific algorithm
  classification model learned earlier using the exact set of features }
}
\details{
  The function first takes the set of features and reduces the matrix dimension according to the given set
  of features that reside in the classification model supplied. It then runs the specific prediction algorithm
  of the given model. The function is indifferent to the model and feature set supplied.
} 
\value{
  Returns a named vector, with samples as names and value 1 or 2 as class prediction for 
  each sample.
}
\references{  }
\author{  }
\note{ 
}
\seealso{\code{\link{selectFeatures}}, \code{\link{selectFeaturesSlimPLS}}, \code{\link{trainClassifier}}, \code{\link{readExpMat} }}
\examples{
# reads an expression matrix with class labels into exp_mat_train
\dontrun{
exp_mat_train <- readExpMat("golub_leukemia_data_with_classes_training.csv", TRUE)}
# selects a set of features into the features2 variable. The matrix we read has a class
# label for each sample in its second row. Labels are either "AML" or "ALL".
# Selection is done using the SlimPLS method. Up to two components with 25 features in
# each component will be selected.
\dontrun{
features2 <- selectFeaturesSlimPLS(exp_mat_train, class_a="AML", class_b="ALL",
				   num_features=50, component_size=25,
				   p_value_threshold=0)}

# train a SVM classifier with linear kernel on the expression matrix using the
# individual feature that are part of the components selected earlier, and are now in
# features2.
\dontrun{
model_t <- trainClassifier("SVM_LINEAR", exp_mat_train, 0,0, "AML", "ALL",
			   features2, FALSE )}

# reads an unlabeled expression matrix exp_mat_test
\dontrun{
exp_mat_test <- readExpMat("golub_leukemia_data_test.csv", FALSE)}

# classify the samples in the test matrix using the model we build:
\dontrun{
getClassification(exp_mat_test, model_t)}
	
# the results should be:
#
# ALL-T-cell-10 ALL-T-cell-11 ALL-T-cell-14   ALLB-cell71   ALLB-cell72  ALL-B-cell-8
#             2             2             2             2             2             2          
# AML60         AML61         AML62         AML63         AML64 
#     1             1             1             1 	      1
	
}
