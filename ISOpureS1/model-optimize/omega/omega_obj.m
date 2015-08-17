function [loglikelihood, deriv_loglikelihood] = omega_obj(ww, tumordata, model)
% function [loglikelihood, deriv_loglikelihood] = omega_obj(ww, tumordata, model)
%
% computes the part of the likelihood function relevant to optimizing
% omega (loglikelihood), as well as the derivative of the likelihood with respect to
% omega (deriv_loglikelihood)

expww = exp(ww);
omega = expww'./ sum(expww);

kappa = model.kappa;
kappaomegaPP = kappa * omega * model.PPtranspose;
omegaPP = omega * model.PPtranspose;
kappaPP = kappa * model.PPtranspose;

loglikelihood = gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + (kappaomegaPP-1)*(model.log_all_rates(end,:)');


dLdomega = sum(kappaPP,2)*psi(0,sum(kappaomegaPP)) - (kappaPP*psi(0,kappaomegaPP')) + kappaPP * model.log_all_rates(end,:)';

deriv_loglikelihood = real(exp(logsum(log(dLdomega') + log(omega),2) + log(omega)));
deriv_loglikelihood = -deriv_loglikelihood + (dLdomega'.*omega);

deriv_loglikelihood = deriv_loglikelihood';

%we fix the first element of the derivative to zero, to fix the scale of
%the unconstrained variables
deriv_loglikelihood(1) = 0;

%we negative likelihood and derivatives because we are using a minimizer
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
