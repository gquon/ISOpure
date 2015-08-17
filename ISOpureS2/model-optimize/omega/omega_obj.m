function [loglikelihood, deriv_loglikelihood] = omega_obj(ww, tumordata, dd, model)
%function [loglikelihood, deriv_loglikelihood] = omega_obj(ww, tumordata, dd, model)
% computes the part of the likelihood function relevant to optimizing
% omega (loglikelihood), as well as the derivative of the likelihood with respect to
% omega (deriv_loglikelihood)

expww = exp(ww');
omega = expww./ sum(expww);
kappa = model.kappa(dd);
kappaomegaPP = kappa .* (omega * model.PPtranspose);
kappaPP = kappa .* model.PPtranspose;


loglikelihood = gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + (kappaomegaPP-1)*(model.log_cc(dd,:)');

dLdomega = sum(kappaPP,2)*psi(0,sum(kappaomegaPP)) - (kappaPP*psi(0,kappaomegaPP')) + kappaPP * model.log_cc(dd,:)';

deriv_loglikelihood = real(exp(logsum(log(dLdomega') + log(omega),2) + log(omega)));

%change of variables
deriv_loglikelihood = -deriv_loglikelihood + (dLdomega'.*omega);

deriv_loglikelihood = deriv_loglikelihood';

%set the first derivative to be zero, to set the scale of the w's
deriv_loglikelihood(1) = 0;

%rasmussen's conjugate gradient method minimizes, so we negative the
%derivative and log likelihood
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
