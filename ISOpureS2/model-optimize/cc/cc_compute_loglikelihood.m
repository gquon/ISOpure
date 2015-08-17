function [loglikelihood] = cc_compute_loglikelihood(log_cancer_rates, tumordata, dd, model)
%function [loglikelihood] = cc_compute_loglikelihood(log_cancer_rates, tumordata, dd, model)
% computes the part of the loglikelihood function relevant to optimizing
% cc

kappaomegaPP = model.omega(dd,:) * model.PPtranspose .* model.kappa(dd);


W = size(model.log_BBtranspose,2);
K = size(model.theta,2);
D = size(tumordata,2);
log_all_rates = [model.log_BBtranspose; log_cancer_rates];


loglikelihood = (kappaomegaPP-1) * log_cancer_rates';


log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + log_all_rates    , 1);
loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));