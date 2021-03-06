% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/external_resources_input.R
\name{getURL}
\alias{getURL}
\title{getURL}
\usage{
getURL(
  URL,
  N.TRIES = 2L,
  read_method = c("delim", "lines", "html"),
  skip = 0,
  delim = "\\t",
  col_names = TRUE
)
}
\arguments{
\item{URL}{Indicated the url that will be read}

\item{N.TRIES}{Integer, how often should the function try to read the URL?}

\item{read_method}{Method to be used for reading of URL content 
("delim" -> \code{readr::read_delim}, 
"lines" -> \code{readr::read_lines}, "html" -> \code{xml2::read_html})}

\item{skip}{integer indicating how many lines to skip when reading URL}

\item{delim}{pattern used for delim 
(passed to \code{delim} of read functions)}

\item{col_names}{list of colnames to use}
}
\value{
returns a the content of the URL. The format of the return object
depends on the read_method that was used.
}
\description{
getURL
}
\keyword{internal}
