\name{extrapolate}
\alias{extrapolate}
\alias{extrapolate.robust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extrapolate nematode community composition from a subsample of identified individuals
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extrapolate(df, counts)
extrapolate.robust(df, counts, n = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
A sites by taxa matrix or dataframe, containing the numbers (or frequencies) of 
identified nematodes. Either way, it will be standardised by dividing with row totals.
}
  \item{counts}{
A vector with the number of counted nematodes in each site.
}
  \item{n}{
number of iterations
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
\itemize{
\item extrapolate() simply multiplies the frequency of each taxon in the 
subsample with the total number of nematodes in each site. This implicitly 
assumes that the composition of the identified subsample is perfectly 
representative of the total.
\item extrapolate.robust() is a bootstraping approach (I think). It treats the 
frequencies of taxa in the identified subsample as probabilities. It then draws 
randomly, with replacement, from the possible taxa, a number equal to the counted 
nematodes, where the probability of each taxon being drawn is its frequency in 
the subsample. The process is iterated n times yielding n versions of each site.
}

}
\value{
\item{extrapolate()}{ 
returns a matrix of the same dimensions as the df.
}

\item{extrapolate.robust()}{
returns a list of n such matrices.
}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angelos Amyntas
}
\note{
%%  ~~further notes~~
What you should NOT do with the output of .robust (or any downstream product 
such as the Enrichment Index), is treat the iterations as levels of a random 
effect in a regression; this would artificially inflate the number of observations, 
likely making any difference, however small and bogus, highly significant! 

What you should do is open to discussion; one simple idea is to fit n individual 
regressions and examine the distribution of the lower and upper bounds of their 
confidence intervals.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# community matrix with frequencies
example1 = matrix(c(.1,0,.2,.7,0,0,.4,.6,.8,.1,.1,0),
                  byrow = T, nrow = 3, ncol = 4,
                  dimnames = list(c("site1","site2","site3"),
                                  c("sp1","sp2","sp3","sp4")))

# community matrix with raw counts
example2 = matrix(c(10,0,20,70,0,0,40,60,80,10,10,0),
                  byrow = T, nrow = 3, ncol = 4,
                  dimnames = list(c("site1","site2","site3"),
                                  c("sp1","sp2","sp3","sp4")))

# counted nematodes
counts = c(132,312,213)

test1 = extrapolate.robust(example1,counts)
test2 = extrapolate.robust(example2,counts)
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
