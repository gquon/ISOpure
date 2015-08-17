function [loglikelihood, deriv_loglikelihood] = kappa_obj(log_kappa, model)
%function [loglikelihood, deriv_loglikelihood] = kappa_obj(log_kappa, model)
% computes the part of the likelihood function relevant to optimizing
% kappa (loglikelihood), as well as the derivative of the likelihood with respect to
% kappa (deriv_loglikelihood)


kappa = exp(log_kappa') + model.MIN_KAPPA;
expww = exp(log_kappa);

omegaPP = model.omega*model.PPtranspose;
kappaomegaPP = omegaPP .* repmat(kappa',1, size(model.PPtranspose,2));

D = size(model.log_cc,1);
W = size(model.log_cc,2);

loglikelihood= 0;

for dd=1:D
    loglikelihood = loglikelihood + gammaln(sum(kappaomegaPP(dd,:))) - sum(gammaln(kappaomegaPP(dd,:)))   +   (kappaomegaPP(dd,:)-1)*model.log_cc(dd,:)';
end

deriv_loglikelihood = zeros(D,1);

for dd=1:D
    deriv_loglikelihood(dd) = (   psi(0,kappa(dd)) - (omegaPP(dd,:)*psi(0,kappaomegaPP(dd,:)')) )       +  ( omegaPP(dd,:)*model.log_cc(dd,:)');
end

%optimizing with respect to log(kappa-(MIN_KAPPA)), so use chain rule to compute
%derivative
%dL/d(ww) = dL/dkappa * dkappa/d(ww)
%                = dL/dkappa * exp(ww)
deriv_loglikelihood = deriv_loglikelihood .* expww;

%rasmussen's conjugate gradient method minimizes, so we negative the derivative and log likelihood
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
