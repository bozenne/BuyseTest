### normexp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  6 2020 (14:06) 
## Version: 
## Last-Updated: maj  6 2020 (20:18) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## cumulative distribution fuction for Z = X + \rho Y
## where X follows a standard normal distribution
##   and Y an exponential distribution with rate parameter \lambda
## denoting \Phi the cumulative distribution function of the standard normal distribution:
## F_Z(z) = \Prob[X + \rho Y < z]
##        = \int f(x,y) \Ind[x + \rho y < z] dx dy
##        = \int_{x \in [-\inf;z]}  f(x) \int_{y \in [0;(z-x)/\rho] f(y) dx dy
##        = \int_{x \in [-\inf;z]}  f(x) (1-\exp(-(z-x)*\lambda/\rho)) dx
##        = \Phi(z) - \exp(-z\lambda/rho)/\sqrt{2\pi} \int_{x \in [-\inf;z]} \exp(x^2/2)\exp(x*\lambda/\rho) dx
##        = \Phi(z) - \exp(-z\lambda/rho+\lambda^2/(2\rho^2))/\sqrt{2\pi} \int_{x \in [-\inf;z]} \exp(-(x-\lambda/rho)^2/2) dx
##        = \Phi(z) - \exp(-z\lambda/rho+\lambda^2/(2\rho^2)) \Phi(z-\lambda/\rho)

## examples
## c(qnormexp(0.5, rate = 2, rho = 1.5), quantile(rnorm(1e5) + 1.5 * rexp(1e5, rate = 2), 0.5))
## c(qnormexp(0.95, rate = 1/10, rho = 1.5), quantile(rnorm(1e5) + 1.5 * rexp(1e5, rate = 1/10), 0.95))
pnormexp <- function(q, rate, rho){
    if(abs(rho)<1e-12){
        out <- stats::pnorm(q)
    }else{
        out <- stats::pnorm(q) - exp(-(rate/rho)*q+(rate/rho)^2/2)*stats::pnorm(q, mean = rate/rho)
    }
    return(out)
}
qnormexp <- function(p, rate, rho){
    if(abs(rho)<1e-12){
        out <- stats::qnorm(p)
    }else{
        out <- sapply(p, function(iP){
            stats::uniroot(function(x){pnormexp(x, rate = rate, rho = rho) - iP},
                           lower = stats::qnorm(iP),
                           upper = (stats::qnorm(iP)+3) + (stats::qexp(iP, rate = rate) + 5/rate))$root
        })
        ## iP <- tail(p,1)
        ## pnormexp(stats::qnorm(iP), rate = rate, rho = rho) - iP
        ## pnormexp((stats::qnorm(iP)+3) + (stats::qexp(iP, rate = rate) + 5/rate), rate = rate, rho = rho) - iP
        ## hist(rnorm(1e4) + rho * rexp(1e4, rate = rate))

    }
    return(out)
}


##----------------------------------------------------------------------
### normexp.R ends here
