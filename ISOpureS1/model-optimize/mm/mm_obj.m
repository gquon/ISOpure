function [loglikelihood, deriv_loglikelihood] = mm_obj(ww, tumordata, model)
%function [loglikelihood, deriv_loglikelihood] = mm_obj(ww, tumordata, model)
%
% computes the part of the likelihood function relevant to optimizing
% mm (loglikelihood), as well as the derivative of the likelihood with respect to
% mm (deriv_loglikelihood)

log_cancer_rates = ww' - logsum(ww,1);
expww = exp(ww');

kappaomegaPP = model.kappa * model.omega' * model.PPtranspose;
W = size(model.log_BBtranspose,2);
K = size(model.theta,2);
D = size(tumordata,2);
log_all_rates = [model.log_BBtranspose; log_cancer_rates];


loglikelihood = (kappaomegaPP-1) * log_cancer_rates';

for dd=1:D
        log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + log_all_rates    , 1);
        loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));
end


%derivative is computed not with respect to mm directly, but with respect
%to unconstrained variables via change of variables
dLdb = zeros(size(log_cancer_rates)) + real(exp(log(kappaomegaPP-1)-log_cancer_rates));


for dd=1:D
    log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + log_all_rates    , 1);
	dLdb = dLdb +  exp(  log(tumordata(:,dd)') + log(model.theta(dd,end)) - log_P_t_given_theta     );
end

bb = expww ./sum(expww);
deriv_loglikelihood = real(exp(logsum(log(dLdb) + log(bb),2) + log(bb))); 
deriv_loglikelihood = -deriv_loglikelihood + (dLdb.*bb);

%we fix the first element of the derivative to zero, to fix the scale of
%the unconstrained variables
deriv_loglikelihood(1) = 0;
deriv_loglikelihood = deriv_loglikelihood';

if (nnz(imag(deriv_loglikelihood)) > 0 | nnz(imag(loglikelihood)) > 0)
        error('imaginary number returned from rminimize in opt_mm');
end

%we negative likelihood and derivatives because we are using a minimizer
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
