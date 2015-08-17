function loglikelihood = omega_compute_loglikelihood(omega, tumordata, model)
% function loglikelihood = omega_compute_loglikelihood(omega, tumordata, model)
% computes the part of the loglikelihood function relevant to optimizing
% omega 

kappaomegaPP = model.kappa * omega' * model.PPtranspose;

loglikelihood = gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + (kappaomegaPP-1)*(model.log_all_rates(end,:)');
