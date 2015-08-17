function loglikelihood = theta_compute_loglikelihood(theta, tumordata, dd, model)
%function loglikelihood = theta_compute_loglikelihood(theta, tumordata, dd, model)
% computes the part of the loglikelihood function relevant to optimizing
% theta

loglikelihood = (model.vv-1) * log(theta');
log_all_rates = [model.log_BBtranspose; model.log_cc(dd,:)];
W = size(log_all_rates,2);
K = size(log_all_rates,1);

%add ln P(t_d|\theta,\beta)
log_P_t_given_theta = logsum(   repmat( log(theta), W, 1)' + log_all_rates    , 1);
loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));
