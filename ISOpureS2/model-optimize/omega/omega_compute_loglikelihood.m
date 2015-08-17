function loglikelihood = omega_compute_loglikelihood(omega, tumordata, dd, model)
%function loglikelihood = omega_compute_loglikelihood(omega, tumordata, dd, model)
% computes the part of the loglikelihood function relevant to optimizing
% omega
kappa = model.kappa(dd);
kappaomegaPP = kappa * omega * model.PPtranspose;

loglikelihood = gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + (kappaomegaPP-1)*(model.log_cc(dd,:)');
