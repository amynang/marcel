\name{query_nemaplex}
\alias{query_nemaplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Query NEMAPLEX to retrive information about nematode taxa
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
query_nemaplex(taxa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxa}{
%%     ~~Describe \code{x} here~~
A vector with the taxa for which the database will be queried.
}
}
\details{
If the input contains taxa not present in Nemaplex, these will appear on the console
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list, containing:
   \item{gen}{
A dataframe with Ecophysiological Parameters for the genera in 'taxa'
}
   \item{fam}{
A dataframe with Ecophysiological Parameters for the families in 'taxa'
}


}
\references{
Ferris, Howard. "Nemaplex". UC-Davis <http://nemaplex.ucdavis.edu/>
}
\author{
Angelos Amyntas
}
\note{
%%  ~~further notes~~
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
