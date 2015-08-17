function [loglikelihood, deriv_loglikelihood] = theta_obj(ww, tumordata, dd, model)

expww = exp(ww');
theta = expww ./ sum(expww);

loglikelihood = (model.vv-1) * log(theta');
W = size(model.log_all_rates,2);
K = size(model.log_all_rates,1);

log_P_t_given_theta = logsum(   repmat( log(theta), W, 1)' + model.log_all_rates    , 1);
loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));


dLdtheta = ((model.vv-1)./theta) +   sum(  exp(model.log_all_rates - repmat(log_P_t_given_theta, K,1))' .* repmat(tumordata(:,dd),1,K)         );


bb = expww ./sum(expww);
deriv_loglikelihood = real(exp(logsum(log(dLdtheta) + log(bb),2) + log(bb)));
deriv_loglikelihood = -deriv_loglikelihood + (dLdtheta.*bb);

deriv_loglikelihood = deriv_loglikelihood';

deriv_loglikelihood(1) = 0;

loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
