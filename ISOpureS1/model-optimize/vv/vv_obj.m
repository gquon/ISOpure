function [loglikelihood, deriv_loglikelihood] = vv_obj(ww, sum_log_theta, DD)
%function [loglikelihood, deriv_loglikelihood] = vv_obj(ww, sum_log_theta, DD)
%
% computes the part of the likelihood function relevant to optimizing vv
% (loglikelihood), as well as the derivative of the likelihood with respect to
% vv (deriv_loglikelihood)

ww = ww';
a = exp(ww) + 1;

loglikelihood = DD * (gammaln(sum(a)) - sum(gammaln(a))) + ((a - 1)*sum_log_theta');

deriv_loglikelihood = (DD * (psi(0,sum(a)) - psi(0,a)) +  sum_log_theta) .* exp(ww);

%we negative likelihood and derivatives because we are using a minimizer
deriv_loglikelihood = deriv_loglikelihood';
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
