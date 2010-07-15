\name{readExpMat}
\alias{readExpMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Reads an expression matrix file to be used for feature selection and classification }
\description{
Read an expression matrix from a file into a expMat object.
Feature labels are expected to be in first column
Sample labels are expected to be in first row
Optionally, the matrix can also hold classes for the different samples. If exist, they
are expected to be in the second row
File should be tab delimited.

}
\usage{
readExpMat(file_name, with_classes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{ \code{file_name} is a tab-delimited matrix of expression values with 
  headers for the rows (features/probes/genes) and for the columns (samples). }
  \item{with_classes}{ If \code{with_classes} is \code{TRUE}, the second row should include
  a class label for each sample }
}
\details{
  
}
\value{
  Returns an object of type \code{expMat}
}
\references{  }
\author{ Michael Gutkin, Ofer Lavi }

\seealso{ \code{\linkS4class{expMat}}}
\examples{
# reads an expression matrix with no class labels into exp_mat
\dontrun{
exp_mat <- readExpMat("golub_leukemia_data_training.csv", FALSE)}
# reads an expression matrix with class labels into exp_mat2
\dontrun{
exp_mat2 <- readExpMat("golub_leukemia_data_with_classes_training.csv", TRUE)}
}

