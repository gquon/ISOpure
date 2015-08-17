function [loglikelihood, deriv_loglikelihood] = vv_obj(ww, sum_log_theta, D)
%function [loglikelihood, deriv_loglikelihood] = vv_obj(ww, sum_log_theta, D)
% computes the part of the likelihood function relevant to optimizing
% vv (loglikelihood), as well as the derivative of the likelihood with respect to
% vv (deriv_loglikelihood)

ww = ww';
vv = exp(ww) + 1;

loglikelihood = D * (gammaln(sum(vv)) - sum(gammaln(vv))) + ((vv - 1)*sum_log_theta');

%optimizing with respect to log(vv-1), so use chain rule to compute
%derivative
%dL/d(ww) = dL/dvv * dvv/d(ww)
%                = dL/dvv * exp(ww)
deriv_loglikelihood = (D * (psi(0,sum(vv)) - psi(0,vv)) +  sum_log_theta) .* exp(ww);

deriv_loglikelihood = deriv_loglikelihood';
%rasmussen's conjugate gradient method minimizes, so we negative the derivative and likelihood
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
