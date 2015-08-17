function loglikelihood = vv_compute_loglikelihood(vv, sum_log_theta, DD)
% function loglikelihood = vv_compute_loglikelihood(vv, sum_log_theta, DD)
% computes the part of the loglikelihood function relevant to optimizing vv

loglikelihood = DD * (gammaln(sum(vv)) - sum(gammaln(vv))) + ((vv - 1)*sum_log_theta');
