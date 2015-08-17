function loglikelihood = theta_compute_loglikelihood(theta, tumordata, dd, model)
% computes portion of log likelihood function that is relevant to
% optimization of theta
loglikelihood = (model.vv-1) * log(theta');
W = size(model.log_all_rates,2);
K = size(model.log_all_rates,1);

log_P_t_given_theta = logsum(   repmat( log(theta), W, 1)' + model.log_all_rates    , 1);
loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));
