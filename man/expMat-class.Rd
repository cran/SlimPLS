\name{expMat-class}
\docType{class}
\alias{expMat-class}

\title{Class "expMat" }
\description{The class represent a general expression matrix, which is an output of a series
of micro-array based experiments.
Each column represent the result of a single experiment. At
each row of this column, the matrix holds the expression value of a specific feature. An experiment
is also known as a condition, or a sample, and a feature might be a probe on the micro-array, a
gene represented by the probe, or any other normalized feature. All values are numeric. The
matrix also holds the set of features, and set of samples, and possibly a set of labels, one for
each sample.}
\section{Objects from the Class}{
Objects can be created by calls of the form 



\code{new("expMat", ...)}



where one can initialize the set of samples, features, the data itself, and possibly a set of class
labels for the samples.
However, the more convenient way of creating such an object is by calling 



\code{readExpMat("file.txt")}



which reads an expression matrix from a file and returns an \code{expMat} object. 
}
\section{Slots}{
	 \describe{
    \item{\code{data}:}{Object of class \code{"matrix"} - the data itself, with samples in rows and 
    features in columns }
    \item{\code{features}:}{Object of class \code{"vector"} - the set of features. Strings that
    usually represent porbes, genes .}
    \item{\code{samples}:}{Object of class \code{"vector"} - the set of samples. String that represent
    a single experiment. Could be a replicate number in a time series experiment, a patient in a disease-control
    experiment, or any other name that represent the specific micro-array result}
    \item{\code{classes}:}{Object of class \code{"vector"} - An optional set of labels for the samples.
    Could be a time in a time series, or a patient condition (sick/healthy) in a disease-control
    experiment" }
  }
}
\section{Methods}{
No methods defined with class "expMat" in the signature.
}
\references{  }
\author{ Michael Gutkin, Ofer Lavi  }
\note{ }


\seealso{ \code{\link{selectFeatures}}, \code{\link{selectFeaturesSlimPLS}}, \code{\link{getClassification}}, \code{\link{readExpMat}}, code{\linkS4class{featureSet}}, \code{\linkS4class{classificationModel}}
}
\examples{
}
\keyword{classes}
