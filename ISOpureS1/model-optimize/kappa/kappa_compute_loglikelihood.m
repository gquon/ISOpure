function loglikelihood = kappa_compute_loglikelihood(kappa, tumordata, model)
% function loglikelihood = kappa_compute_loglikelihood(kappa, tumordata, model)
% computes the part of the loglikelihood function relevant to optimizing kappa
% (loglikelihood)

omega = model.omega;
kappaomegaPP = kappa * omega' * model.PPtranspose;

loglikelihood = gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + (kappaomegaPP-1)*(model.log_all_rates(end,:)');
