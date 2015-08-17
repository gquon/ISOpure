function loglikelihood = vv_compute_loglikelihood(vv, sum_log_theta, D)
%function loglikelihood = vv_compute_loglikelihood(vv, sum_log_theta, D)
% computes the part of the loglikelihood function relevant to optimizing
% vv

loglikelihood = D * (gammaln(sum(vv)) - sum(gammaln(vv))) + ((vv - 1)*sum_log_theta');
