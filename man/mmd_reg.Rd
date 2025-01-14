% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reg_function.R
\name{mmd_reg}
\alias{mmd_reg}
\title{MMD regression}
\usage{
mmd_reg(y, X, model, intercept, par1, par2, kernel.y, kernel.x, bdwth.y, bdwth.x, 
        control= list())
}
\arguments{
\item{y}{Response variable. Must be a vector of length  \eqn{n\geq 3}.}

\item{X}{Design matrix. \code{X} must be either a matrix of dimension \eqn{n\times p} or a vector of size  \eqn{n},  where \eqn{n} is the size of \code{y}.}

\item{model}{Regression model to be fitted to the data. By default, the linear regression model with \eqn{\mathcal{N}_1(0,\phi^2)} error terms is used.  See details for the list of available models.}

\item{intercept}{If \code{intercept=TRUE} an intercept is added to the model, while no intercept is added if \code{intercept=FALSE}. By default, \code{intercept=TRUE}.}

\item{par1}{Values of the regression coefficients of the model used as starting values to numerically optimize the objective function. \code{par1} must be either a vector of size \eqn{p}, with \eqn{p} the number of columns of \code{X} (or with \eqn{p=1} if \code{X} is a vector), or equal to \code{"auto"}, in which case a non-robust estimate of the regression coefficients of the model is used as starting values. By default, \code{par1="auto"}.}

\item{par2}{A value for the auxilliary parameter \eqn{\phi} of the model. \code{par2} needs to be specified only if relevant (see details for the list of models having an  auxilliary parameter \eqn{\phi}).  If the model assumes that \eqn{\phi}  is known (see details) then \code{par2} must be a strictly positive real number and the model is estimated with \eqn{\phi=}\code{par2}. If the model assumes that \eqn{\phi} is unknown (see details)  then the value specified by   \code{par2}  is used as starting value to numerically optimize the objective function. For such models \code{par2} must be either a strictly positive real number or equal to \code{"auto"}, in which case a non-robust estimate of \eqn{\phi} is used as starting value. By default, \code{par2="auto"}.}

\item{kernel.y}{Kernel applied on the response variable. Available options for  \code{kernel.y} are \code{"Gaussian"} (Gaussian kernel), \code{"Laplace"} (Laplace, or exponential, kernel) and \code{"Cauchy"} (Cauchy kernel). By default, \code{kernel.y="Gaussian"} for the linear regression model and \code{kernel.y="Laplace"} for the other models.}

\item{kernel.x}{Kernel applied on the explanatory variables. Available options for  \code{kernel.x} are \code{"Gaussian"} (Gaussian kernel), \code{"Laplace"} (Laplace, or exponential, kernel) and \code{"Cauchy"} (Cauchy kernel). By default, \code{kernel.x="Laplace"}.}

\item{bdwth.y}{Bandwidth parameter for the kernel \code{kernel.y}. \code{bdwth.y} must be eiter a strictly positive real number or equal to \code{"auto"}, in which case the median heuristic  is used to select the bandwidth parameter of \code{kernel.y} (see details). By default, \code{bdwth.y="auto"}.}

\item{bdwth.x}{Bandwidth parameter for the kernel \code{kernel.x}. \code{bdwth.x} must be either a non-negative real number or equal to \code{"auto"}, in which case a rescaled version of the median heuristic is used to specify the bandwidth parameter  of \code{kernel.x} (see details). By default,   \code{bdwth.x}=0. Remark: for computational reasons, for large dataset (i.e.~when the sample size is bigger than a few thousands) it is recommended to choose \code{bdwth.x}=0 (see details).}

\item{control}{A \code{list} of control parameters for the numerical optimization of the objective function. See details.}
}
\value{
\code{MMD_reg} returns an object of \code{class} \code{"regMMD"}.

The function \code{summary}  can be used to print the results.

An object of class \code{regMMD} is a list containing the following components:
\item{coefficients}{Estimated regression coefficients.}
\item{intercept}{If \code{intercept}=TRUE an intercept has been added to model, if \code{intercept}=FALSE no intercept has been added.}
\item{phi}{If relevant (see details), either the estimated value of the \eqn{\phi} parameter of model, or the value of \eqn{\phi} used to fit the model if \eqn{\phi} is assumed to be known.}
\item{kernel.y}{Kernel applied on the response variable used to fit the model.}
\item{bdwth.y}{Value of the bandwidth for the kernel applied on the response variable used to fit the model.}
\item{kernel.x}{Kernel applied on the explanatory variables used to fit the model.}
\item{bdwth.x}{Value of the bandwidth for the kernel applied on the explanatory variables used to fit the model.}
\item{par1}{Value of the parameter \code{par1} used to fit the model.}
\item{par2}{Value of  parameter \code{par2} used to fit the model.}
\item{trace}{If the control parameter \code{trace=TRUE}, \code{trace} is a matrix containing the parameter values obtained at the end of each iteration of the optimization algorithm.}
}
\description{
Fits a regression model to the data, using the robust procedure based on maximum mean discrepancy (MMD) minimization introduced and studied in \insertCite{MMD;textual}{regMMD}.
}
\details{
Available options for \code{model} are:
\describe{
\item{\code{"linearGaussian"}}{Linear regression model with \eqn{\mathcal{N}_1(0,\phi^2)} error terms, with \eqn{\phi} unknown.}
\item{\code{"linearGaussian.loc"}}{Linear regression model with \eqn{\mathcal{N}_1(0,\phi^2)} error terms, with \eqn{\phi} known.}
\item{\code{"gamma"}}{Gamma regression model with unknown shape parameter \eqn{\phi}. The inverse function is used as link function.}
\item{\code{"gamma.loc"}}{Gamma regression model with known shape parameter \eqn{\phi}. The inverse function is used as link function.}
\item{\code{"beta"}}{Beta regression model with unknown precision parameter \eqn{\phi}. The logistic function is used as link function.}
\item{\code{"beta.loc"}}{Beta regression model with  known precision parameter \eqn{\phi}. The logistic function is used as link function.}
\item{\code{"logistic"}}{Logistic regression model.}
\item{\code{"exponential"}}{Exponential regression.}
\item{\code{"poisson"}}{Poisson regression model.}
}

When \code{bdwth.x}>0 the function \code{reg_mmd} computes the estimator \eqn{\hat{\theta}_n} introduced in \insertCite{MMD;textual}{regMMD}. When \code{bdwth.x}=0  the function \code{reg_mmd} computes the estimator \eqn{\tilde{\theta}_n} introduced in \insertCite{MMD;textual}{regMMD}. The former estimator has stronger theoretical properties but is more expensive to compute (see below).

When \code{bdwth.x}=0 and \code{model} is  \code{"linearGaussian"}, \code{"linearGaussian.loc"} or \code{"logistic"}, the objective function and its gradient can be computed on \eqn{\mathcal{O}(n)} operations, where \eqn{n} is the sample size (i.e. the dimension of \code{y}). In this case, gradient descent with backtraking line search is used to perform the minimizatiom. The algorithm stops  when the maximum number of iterations \code{maxit} is reached, or as soon as the change in the objective function is less than \code{eps_gd} times the current function value. In the former case, a warning message is generated. By defaut, \code{maxit}=\eqn{5\times 10^4} and  \code{eps_gd=sqrt(.Machine$double.eps)}, and the value of these two parameters can be changed using the \code{control} argument of \code{reg_mmd}.

When \code{bdwth.x}>0 and \code{model} is  \code{"linearGaussian"}, \code{"linearGaussian.loc"} or \code{"logistic"}, the objective function and its gradient can be computed on \eqn{\mathcal{O}(n^2)} operations. To reduce the computational cost the objective function is minimized using norm adagrad \insertCite{duchi2011adaptive}{regMMD}, an adaptive step size stochastic gradient algorithm. Each iteration of the algorithm requires \eqn{\mathcal{O}(n)} operations. However,  the algorithm has an intialization step that requires \eqn{\mathcal{O}(n^2)} operations and has a memory requirement of size \eqn{\mathcal{O}(n^2)}.

When   \code{model} is not in \code{c("linearGaussian", "linearGaussian.loc", "logistic")},  the objective function and its gradient cannot be computed explicitly and the minimization is performed using norm adagrad. The cost per iteration of the algorithm is \eqn{\mathcal{O}(n)} but, for  \code{bdwth.x}>0, the memory requirement and the initialization cost are both of size \eqn{\mathcal{O}(n^2)}.

When adagrad is used, \code{burnin} iterations  are performed as a warm-up step.  The algorithm then stops when \code{burnin}+\code{maxit} iterations are performed, or as soon as the norm of the average value of the gradient evaluations computed in all the previous iterations is less than \code{eps_sg}. A warning message is generated if the maximum number of iterations is reached. By default, \code{burnin}=\eqn{10^3}, \code{nsteps}=\eqn{5\times 10^4} and \code{eps_sg}=\eqn{10^{-5}} and  the value of these three parameters can be changed using the \code{control} argument of  \code{reg_mmd}.

If \code{bdwth.y="auto"} then the value of the bandwidth parameter of \code{kernel.y} is equal to \eqn{H_n/\sqrt{2}}  with \eqn{H_n}   the median value of the set \eqn{ \{|y_i-y_j|\}_{i,j=1}^n}, where \eqn{y_i} denote the ith component of \code{y}. This definition of \code{bdwth.y} is motivated by the results in \insertCite{garreau2017large;textual}{regMMD}. If \eqn{H_n=0} the bandwidth parameter of \code{kernel.y} is set to 1.

If \code{bdwth.x="auto"} then the value of the bandwidth parameter of \code{kernel.x} is equal to \eqn{0.01H_n/\sqrt{2}}  with \eqn{H_n} is the median value of the set \eqn{ \{\|x_i-x_j\|\}_{i,j=1}^n}, where \eqn{x_i} denote the ith row of the design matrix \code{X}. If \eqn{H_n=0} the bandwidth parameter of \code{kernel.x} is set to 1.

The \code{control} argument is a list  that can supply any of the following components:
\describe{
\item{rescale:}{If \code{rescale=TRUE} the (non-constant) columns of \code{X} are rescalled before perfoming the optimization, while if \code{rescale=FLASE} no rescaling is applied. By default \code{rescale=TRUE}.}
\item{burnin}{A non-negative integer.}
\item{eps_gd}{A non-negative real number.}
\item{eps_sg}{A non-negative real number.}
\item{maxit}{A integer strictly larger than 2.}
\item{stepsize}{Scaling constant for the step-sizes used by adagrad. \code{stepsize} must be a stictly positive number and by default \code{stepsize}=1.}
\item{trace:}{If \code{trace=TRUE} then the parameter value obtained at the end  of  each iteration  (after the burn-in perdiod for adagrad)  is returned. By default, \code{trace=TRUE} and   \code{trace} is automatically set to \code{TRUE} if the maximum number of iterations is reached.}
\item{epsilon}{Parameter used in adagrad to avoid numerical errors in the computation of the step-sizes. \code{epsilon} must be a strictly positive real number and by default \code{epsilon}=\eqn{10^{-4}}.}
\item{alpha}{Parameter for the backtraking line search. \code{alpha} must be a real number in \eqn{(0,1)} and by default \code{alpha}=0.8.}
\item{c_det}{Parameter used to control the computational cost of the algorithm when \code{gamma.x}\eqn{>0}, see the Suplementary material in \insertCite{MMD;textual}{regMMD} for mode details. \code{c_det} must be a real number in \eqn{(0,1)} and by default \code{c_det}=0.2.}
\item{c_rand}{Parameter used to control the computational cost of the algorithm when \code{bdwth.x}\eqn{>0}, see the Suplementary material in \insertCite{MMD;textual}{regMMD} for mode details. \code{c_rand} must be a real number in \eqn{(0,1)} and by default \code{c_rand}=0.1.}
}
}
\examples{
#Simulate data
n<-1000
p<-4
beta<-1:p
phi<-1
X<-matrix(data=rnorm(n*p,0,1),nrow=n,ncol=p)
data<-1+X\%*\%beta+rnorm(n,0,phi)

##Example 1: Linear Gaussian model 
y<-data
Est<-mmd_reg(y, X)
summary(Est)

##Example 2: Logisitic regression model
y<-data
y[data>5]<-1
y[data<=5]<-0
Est<-mmd_reg(y, X, model="logistic")
summary(Est)
Est<-mmd_reg(y, X, model="logistic", bdwth.x="auto")
summary(Est)
}
\references{
\insertAllCited{}
}
