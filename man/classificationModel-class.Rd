\name{classificationModel-class}
\docType{class}
\alias{classificationModel-class}

\title{Class "classificationModel" }
\description{The class wraps a model learned using a specific classification algorithm,
that can later be used to predict the class of a sample. It holds the specific algorithm data
as well as the set of features used while learning the model. The set of features reduces the
dimension of the problem by selecting or extracting a small number of features from a much 
larger number of basic features.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("classificationModel", ...)}, but are 
usually created by calling \code{"trainClassifier"}, with a given feature set, data and algorithm.
The object can later be used by calling \code{"getClassification"} for predicting classes of
unseen samples. 
}
\section{Slots}{
	 \describe{
    \item{\code{model_type}:}{Object of class \code{"character"} - a string the represents the specific classification algorithm.

	
    Can be one of \code{"SVM_LINEAR"}, \code{"SVM_RADIAL"} or \code{"RandomForest"}}
    \item{\code{model}:}{Object of class \code{"svmOrRF"} - the model itself }
    \item{\code{features}:}{Object of class \code{"featureSet"} - the set of features used for learning the model }
    \item{\code{pls_use_components}:}{Object of class \code{"logical"} - whether the model was learned by using 
    linear combinations (components) of basic features or not. }
  }
}
\section{Methods}{
No methods defined with class "classificationModel" in the signature.
}
\references{  }
\author{ Michael Gutkin, Ofer Lavi }
\note{ }
\seealso{ \code{\link{selectFeatures}}, \code{\link{trainClassifier}}, \code{\link{getClassification}}, \code{\link{readExpMat}}, \code{\linkS4class{featureSet}}, \code{\linkS4class{classificationModel}}
}
\examples{
}
\keyword{classes}
