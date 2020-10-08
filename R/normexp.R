### normexp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  6 2020 (14:06) 
## Version: 
## Last-Updated: okt  8 2020 (17:22) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * exponential distribution
## cumulative distribution fuction for Z = X + \rho Y
## where X follows a standard normal distribution
##   and Y an exponential distribution with rate parameter \lambda
## denoting \Phi the cumulative distribution function of the standard normal distribution:
## F_Z(z) = \Prob[X + \rho Y < z]
##        = \int f(x,y) \Ind[x + \rho y < z] dx dy
##        = \int_{x \in [-\inf;z]}  f(x) \int_{y \in [0;(z-x)/\rho] f(y) dx dy
##        = \int_{x \in [-\inf;z]}  f(x) (1-\exp(-(z-x)*\lambda/\rho)) dx
##        = \Phi(z) - \exp(-z\lambda/rho)/\sqrt{2\pi} \int_{x \in [-\inf;z]} \exp(-x^2/2)\exp(x*\lambda/\rho) dx
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

## * Weibull distribution
## cumulative distribution fuction for Z = X + \rho Y
## where X follows a standard normal distribution
##   and Y a weibull distribution with scale parameter \lambda and shape parameter k
## denoting \Phi the cumulative distribution function of the standard normal distribution:
## F_Z(z) = \Prob[X + \rho Y < z]
##        = \int f(x,y) \Ind[x + \rho y < z] dx dy
##        = \int_{x \in [-\inf;z]}  f(x) \int_{y \in [0;(z-x)/\rho] f(y) dx dy
##        = \int_{x \in [-\inf;z]}  f(x) (1-\exp(-(z-x)^k/(\lambda\rho)^k)) dx
##        = \Phi(z) - \int_{x \in [-\inf;z]} \exp(-x^2/2-(z-x)^k/(\lambda\rho)^k)\sqrt{2\pi} dx

## examples
## c(qnormweibull(0.5, scale = 1/2, shape = 1, rho = 1.5), quantile(rnorm(1e5) + 1.5 * rweibull(1e5, scale = 1/2, shape = 1), 0.5))
## c(qnormweibull(0.95, scale = 10, shape = 1, rho = 1.5), quantile(rnorm(1e5) + 1.5 * rweibull(1e5, scale = 10, shape = 1), 0.95))
## 
## c(qnormweibull(0.5, scale = 1/2, shape = 2, rho = 1.5), quantile(rnorm(1e5) + 1.5 * rweibull(1e5, scale = 1/2, shape = 2), 0.5))
## c(qnormweibull(0.95, scale = 10, shape = 2, rho = 1.5), quantile(rnorm(1e5) + 1.5 * rweibull(1e5, scale = 10, shape = 2), 0.95))
pnormweibull <- function(q, scale, shape, rho){
    if(abs(rho)<1e-12){
        out <- stats::pnorm(q)
    }else{
        if(shape==1){
            out <- stats::pnorm(q) - exp(-(1/(scale*rho))*q+(1/(scale*rho))^2/2)*stats::pnorm(q, mean = 1/(scale*rho))
        }else{
            I <- stats::integrate(f = function(x){exp(-x^2/2)/sqrt(2*pi)*exp(-((q-x)/(rho*scale))^shape)}, lower = min(-4,q - 7^(1/shape)*rho*scale), upper = q)
            out <- stats::pnorm(q) - I$value
        }
       
    }
    return(out)
}
qnormweibull <- function(p, scale, shape, rho){
    if(abs(rho)<1e-12){
        out <- stats::qnorm(p)
    }else{
        out <- sapply(p, function(iP){
            stats::uniroot(function(x){pnormweibull(x, scale = scale, shape = shape, rho = rho) - iP},
                           lower = stats::qnorm(iP),
                           upper = (stats::qnorm(iP)+3) + (stats::qweibull(iP, scale = scale, shape = shape) + 5*scale))$root
        })
    }
    return(out)
}

##----------------------------------------------------------------------
### normexp.R ends here
