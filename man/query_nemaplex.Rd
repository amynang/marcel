\name{query_nemaplex}
\alias{query_nemaplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Query NEMAPLEX to retrieve information about nematode taxa
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
query_nemaplex(taxa, complete = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxa}{
%%     ~~Describe \code{x} here~~
A vector with the taxa for which the database will be queried.
}
\item{complete}{
%%     ~~Describe \code{x} here~~
Logical. Default returns a subset of all the information provided by nemaplex. This does not make the query faster, it just makes the output less messy.
}
}
\details{
If the input contains taxa not present in Nemaplex, these will appear on the console.
query_nemaplex() currently uses Firefox, so the browser needs to be installed in 
your computer. It also requires JAVA. For 64bit(32bit) R you need 64bit(32bit) JAVA.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

A dataframe with Ecophysiological Parameters for the genera and families in 'taxa'.



}
\references{
Ferris, Howard. "Nemaplex". UC-Davis <http://nemaplex.ucdavis.edu/>
}
\author{
Angelos Amyntas
}
\note{
%%  ~~further notes~~
  My experience is that the first time you will call the function, RSelenium will 
  do some setup procedure. It is a good idea to start with a small number of taxa
  (1 or 2 are fine) to make sure everything runs ok before you input the complete 
  list of taxa you want to retrieve information for.
  
  If any error occurs that breaks the process, you will need to call
  system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)
  to kill java manually, before you try calling the function again.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# the input can contain genera and/or families

query_nemaplex("Actinolaimidae", 
               "Alaimidae", 
               "Acrobeles", 
               "Acrobeloides")

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
