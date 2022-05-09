\name{Community Indices}
\alias{all.indices}
\alias{Enrichment}
\alias{Structure}
\alias{Channel}
\alias{Maturity}
\alias{Plant.Parasite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
all.indices(x, nemaplex)

Enrichment(x, nemaplex)
Structure(x, nemaplex)
Channel(x, nemaplex)
Maturity(x, nemaplex)
Plant.Parasite(x, nemaplex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A dataframe with locations/samples as rows and nematode taxa (families or genera) as columns. Values can be raw counts, densities or frequencies of the taxa in each sample
}
  \item{nemaplex}{
  The output of query_nemaplex(), containing information about the c-p scores and the feeding habits of the taxa in x
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~

b = (Ba2 + Fu2) * W2, where W2 = 0.8,

e = (Ba1 * W1) + (Fu2 * W2), where W1 = 3.2 and W2 = 0.8

s = (Ban*Wn + Can*Wn + Fun*Wn + Omn*Wn)

where n=3-5, W3.= 1.8, W4.= 3.2, W5.= 5.0

The Enrichment index is calculated as EI = 100 * e / (e + b)

The Structure index is calculated as SI = 100 * s / (s + b)

The Channel index is calculated as CI = 100 * (0.8*Fu2 / (3.2*Ba1 + 0.8*Fu2))
}
\value{
Either a single column dataframe (in the case of individual index functions), or a dataframe with samples as rows and all the indices calculated by maRcel as columns
}
\references{
Ferris, H., Bongers, T. and Goede, R.G.M. de (2001) A framework for soil food web diagnostics:
extension of the nematode faunal analysis concept. Applied Soil Ecology 18, 13-29.

Ferris, H., & Bongers, A. M. T. (2009). Indices Developed Specifically for analysis of Nematode Assemblages. In M. J. Wilson, & T. Kakouli-Duarte (Eds.), Nematodes as Environmental Indicators (pp. 124-145). CAB International. 
}
\author{
Angelos Amyntas
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as 
\code{
  \link{query_nemaplex}
  }
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
