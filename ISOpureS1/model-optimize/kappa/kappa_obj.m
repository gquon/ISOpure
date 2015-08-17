function [loglikelihood, deriv_loglikelihood] = kappa_obj(log_kappa, tumordata, model)
%function [loglikelihood, deriv_loglikelihood] = kappa_obj(log_kappa, tumordata, model)
%
% computes the part of the likelihood function relevant to optimizing kappa
% (loglikelihood), as well as the derivative of the likelihood with respect to
% kappa (deriv_loglikelihood)


kappa = exp(log_kappa) + model.MIN_KAPPA;
omega = model.omega;
kappaomegaPP = kappa * omega' * model.PPtranspose;
omegaPP = omega' * model.PPtranspose;

loglikelihood = gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + (kappaomegaPP-1)*(model.log_all_rates(end,:)');

deriv_loglikelihood = sum(omegaPP)*psi(0,sum(kappaomegaPP)) - (omegaPP*psi(0,kappaomegaPP')) + model.log_all_rates(end,:)*omegaPP';
deriv_loglikelihood = deriv_loglikelihood * exp(log_kappa);

%we negative likelihood and derivatives because we are using a minimizer
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
