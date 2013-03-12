### Definition ###

setClass("Rico", representation(id="integer"))

### Construction ###

## DEBUG
rico.num <- function(x) {
  new("Rico", id = .Call("createRico_Number", x))
}

rico.tri <- function(x, y, z) {
  if(missing(x)) return(new("Rico", id = .Call("create_Triangular")))
  if(missing(y)) return(new("Rico", id = .Call("create_TriangularM", x)))
  if(missing(z)) return(new("Rico", id = .Call("create_TriangularLR", x, y)))
  return(               new("Rico", id = .Call("create_TriangularLMR", x, y, z)))
}

rico.norm <- function(x, y) {
  if(missing(x)) return(new("Rico", id = .Call("create_Normal")))
  if(missing(y)) return(new("Rico", id = .Call("create_Normal_Mean", x)))
  return(               new("Rico", id = .Call("create_Normal_Mean_Stdev", x, y)))
}

rico.chisq <- function(x) {
  new("Rico", id = .Call("create_ChiSquared", x))
}

rico.t <- function(x) {
  new("Rico", id = .Call("create_Students_t", x))
}

rico.f <- function(x, y) {
  new("Rico", id = .Call("create_Fisher_F", x, y))
}

rico.cauchy <- function(x, y) {
  if(missing(x)) return(new("Rico", id = .Call("create_Cauchy")))
  if(missing(y)) return(new("Rico", id = .Call("create_Cauchy_Mean", x)))
  return(               new("Rico", id = .Call("create_Cauchy_Mean_Stdev", x, y)))
}

rico.lnorm <- function(x, y) {
  if(missing(x)) return(new("Rico", id = .Call("create_LogNormal")))
  if(missing(y)) return(new("Rico", id = .Call("create_LogNormal_Mean", x)))
  return(               new("Rico", id = .Call("create_LogNormal_Mean_Stdev", x, y)))
}

rico.beta <- function(x, y) {
  new("Rico", id = .Call("create_Beta", x, y))
}

rico.exp <- function(x) {
  new("Rico", id = .Call("create_Exponential", x))
}

rico.logis <- function(x, y) {
  new("Rico", id = .Call("create_Logistic", x, y))
}

rico.weibull <- function(x, y) {
  if(missing(y)) return(new("Rico", id = .Call("create_Weibull_Shape", x)))
  return(               new("Rico", id = .Call("create_Weibull_Shape_Scale", x, y)))
}

rico.unif <- function(x, y) {
  if(missing(x)) return(new("Rico", id = .Call("create_Uniform_Std")))
  return(               new("Rico", id = .Call("create_Uniform_Left_Right", x, y)))
}

### Display ###

setMethod("show", "Rico", function(object) {
  .Call("show_Rico", object@id)
})

setMethod("plot", signature(x="Rico", y="ANY"),
          function(x, y, ...) {
            r = .Call("getSuggestedPlotRange", x@id)
            ##cat("TRF: Rico.plot[", r[1], ", ", r[2], "]\n", sep="")
            left  <- r[1]
            right <- r[2]
            N     <- 1000    ## default number of data points
            Xo    <- 1:N
            X     <- (Xo-1)/(N-1)*(right-left) + left   ## includes endpoints
            Y     <- .Call("getPlotPoints", x@id, X)
            ##Y     <- X ** 2
            plot(X, Y, type="l", ...)
          })

### Arithmetic Operations ###

setMethod("Arith", signature(e1 = "Rico", e2 = "Rico"),
          function(e1, e2) {
            ##cat("TRF: called Arith[", .Generic, "](Rico,Rico)\n", sep="")
            .rico.arith(e1, e2, .Generic)
          })

.rico.negate <- function(x) {new("Rico", id = .Call("createRico_Negate", x@id))}

setMethod("Arith", signature(e1 = "Rico"),
          function(e1, e2) {
            ##cat("TRF: called Arith[", .Generic, "](Rico,ANY)\n", sep="")
            if(missing(e2)) {
              switch(.Generic,
                     "+" = e1,
                     "-" = .rico.negate(e1),
                     )
            } else {
              .rico.arith(e1, as(e2,"Rico"), .Generic)
            }
          })

setMethod("Arith", signature(e2 = "Rico"),
          function(e1, e2) {
            ##cat("TRF: called Arith[", .Generic, "],(ANY,Rico)\n", sep="")
            .rico.arith(as(e1,"Rico"),e2,.Generic)
          })

setMethod("sqrt", signature(x = "Rico"),
          function(x) {
            ##cat("TRF: called log(Rico)\n", sep="")
            new("Rico", id = .Call("createRico_Sqrt", x@id))
          })

setMethod("log", signature(x = "Rico"),
          function(x) {
            ##cat("TRF: called log(Rico)\n", sep="")
            new("Rico", id = .Call("createRico_Log", x@id))
          })

setMethod("exp", signature(x = "Rico"),
          function(x) {
            ##cat("TRF: called exp(Rico)\n", sep="")
            new("Rico", id = .Call("createRico_Exp", x@id))
          })

setGeneric("D", function(x, y) standardGeneric("D"))
setMethod("D", signature(x = "Rico", y = "missing"),
          function(x,y) {
            new("Rico", id = .Call("differentiateRico_X", x@id))
          })
setMethod("D", signature(x = "Rico", y = "Rico"),
          function(x,y) {
            new("Rico", id = .Call("differentiateRico_WRT", x@id, y@id))
          })

setGeneric("Simplify", function(x) standardGeneric("Simplify"))
setMethod("Simplify", signature(x = "Rico"),
          function(x) {
            new("Rico", id = .Call("simplifyRico", x@id))
          })

## The helper functions do most of decision making.

.rico.number <- function(i) {new("Rico", id = .Call("createRico_Number", i))}

setAs("numeric", "Rico",
      function(from) {
        #cat("TRF: setAs[numeric,Rico] called.\n");           ## DEBUG
        .rico.number(from);
      })

.rico.add      <- function(x,y) {new("Rico", id = .Call("createRico_Add",      x@id, y@id))}
.rico.subtract <- function(x,y) {new("Rico", id = .Call("createRico_Subtract", x@id, y@id))}
.rico.multiply <- function(x,y) {new("Rico", id = .Call("createRico_Multiply", x@id, y@id))}
.rico.divide   <- function(x,y) {new("Rico", id = .Call("createRico_Divide",   x@id, y@id))}
.rico.power    <- function(x,y) {new("Rico", id = .Call("createRico_Power",    x@id, y@id))}

.rico.arith <- function(e1, e2, .Generic) {
  ##cat("TRF: called .rico.arith(",.Generic,")\n",sep="")
  #print(e1); print(e2)
  switch(.Generic, 
         "+"   = .rico.add(e1, e2),
         "-"   = .rico.subtract(e1, e2),
         "*"   = .rico.multiply(e1, e2),
         "/"   = .rico.divide(e1, e2),
         "^"   = .rico.power(e1, e2),
         stop(paste(.Generic, "not allowed on Rico objects."))    ## default
         #"%/%" = ??,
         #"%%"  = ??,
         )
}

## Comparison operations

.rico.compare <- function(x,y) {.Call("compareRico", x@id, y@id)}

.CompareRico <- function(e1, e2) {
  switch(.Generic,
         "==" =   .rico.compare(e1, e2),
	 "!=" = ! .rico.compare(e1, e2),
	 stop("unsupported comparison of Rico Objects"))
}

setMethod("Compare", signature(e1 = "Rico", e2 = "Rico"), .CompareRico)
#setMethod("Compare", signature(e2 = "polynomial"), .ComparePolynomial)

### Testing ###########################################

rico.test <- function() {.C("RicoTestRun"); invisible(0)}

rico.showArgs         <- function(...) invisible(.External("showArgs", ...))
rico.showArgs1        <- function(...) invisible(.Call("showArgs1", list(...)))
rico.testGetList      <- function() .Call("testGetList")
rico.testGetArray     <- function() .Call("testGetArray")
rico.testGetArrayList <- function() .Call("testGetArrayList")
