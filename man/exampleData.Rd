\name{exampleData}
\alias{exampleData}
\docType{data}
\title{ example data to test the functionality of the package}
\description{Example data generated with the generate.Guilds function using the parameters: theta = 200, alpha X = 0.005, alpha Y = 0.001}
\usage{data(exampleData)}
\format{The data is stored in a list with:\cr
     \describe{ 
		\item{guildX}{- Species abundance distribution for guild X}
		\item{guildY}{- Species abundance distribution for guild Y}
     }
}

\examples{
data(exampleData);
parametervals <- c(200,0.005,0.001) #theta = 200, alpha X = 0.005, alpha Y = 0.001
LL = logLikelihood.Guilds(parametervals,model="D1",
			  exampleData$guildX,exampleData$guildY,verbose=TRUE);
}
\keyword{datasets}
