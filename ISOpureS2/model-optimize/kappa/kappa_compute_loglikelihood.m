function loglikelihood = kappa_compute_loglikelihood(kappa, model)
%function loglikelihood = kappa_compute_loglikelihood(kappa, model)
% computes the part of the loglikelihood function relevant to optimizing
% kappa

kappaomegaPP = kappa .* (model.omega*model.PPtranspose);


D = size(model.log_cc,1);
loglikelihood = 0;

for dd=1:D
    loglikelihood = loglikelihood + gammaln(sum(kappaomegaPP(dd,:))) - sum(gammaln(kappaomegaPP(dd,:)))   +   (kappaomegaPP(dd,:)-1)*model.log_cc(dd,:)';
end
