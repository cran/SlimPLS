\name{featureSet-class}
\docType{class}
\alias{featureSet-class}

\title{Class "featureSet" }
\description{Holds a set of features to be used for classification of expression
data. The class is usually the output of a feature selection process on the data. The features
are constructed of basic features which are usually genes or probe ids in an expression matrix.
The set of features can be either a simple vector of names denoting the probes or genes, or a 
set of weighted linear combinations of such simple features. Such a linear combination is known
as a component.}
\section{Objects from the Class}{
Objects can be created by calls of the form 



\code{new("featureSet", ...)}



but are usually created by calling either 



\code{"selectFeatures"} 


or



\code{"selectFeaturesSlimPLS"}.



While the former creates a set of basic features selected from the basic features in the data by various methods,
the latter creates a set of component features, where each component is a linear combination of
basic features from the data. 
}
\section{Slots}{
	 \describe{
    \item{\code{features}:}{Object of class \code{"vector"} - A list of strings representing the features }
    \item{\code{t_mat_learn}:}{Object of class \code{"matrix"} - for internal use  }
    \item{\code{w_mat}:}{Object of class \code{"matrix"} - a weight matrix that represent component features, where
    each column represents a component, and each row represents a basic feature. The values within the matrix are
    the weight contribution of the basic feature in the given row, to the component in the given column.}
    \item{\code{p_mat}:}{Object of class \code{"matrix"} - for internal use }
    \item{\code{q_mat}:}{Object of class \code{"matrix"} - for internal use  }
    \item{\code{is_pls}:}{Object of class \code{"logical"} - a flag whether to use the simple features vector, or the
    components }
  }
}
\section{Methods}{
No methods defined with class "featureSet" in the signature.
}
\references{  }
\author{ Michael Gutkin, Ofer Lavi }
\note{  }

\seealso{ \code{\link{selectFeatures}}, \code{\link{trainClassifier}}, \code{\link{getClassification}}, \code{\link{readExpMat}}
, \code{\linkS4class{featureSet}}, \code{\linkS4class{classificationModel}}
}
\examples{
}
\keyword{classes}
