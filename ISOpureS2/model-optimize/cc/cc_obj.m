function [loglikelihood, deriv_loglikelihood] = cc_obj(ww, tumordata, dd, model)
%function [loglikelihood, deriv_loglikelihood] = cc_obj(ww, tumordata, dd, model)
% computes the part of the likelihood function relevant to optimizing
% cc (loglikelihood), as well as the derivative of the likelihood with respect to
% cc (deriv_loglikelihood)

log_cancer_rates = ww' - logsum(ww,1);
expww = exp(ww');
kappaomegaPP = model.kappa(dd) .* model.omega(dd,:) * model.PPtranspose;

W = size(model.log_BBtranspose,2);
K = size(model.theta,2);
D = size(tumordata,2);
log_all_rates = [model.log_BBtranspose; log_cancer_rates];


loglikelihood = (kappaomegaPP-1) * log_cancer_rates';

log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + log_all_rates    , 1);
loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));

%change of variables
%need to wrap exp() with the real() because log(kappaomegaPP-1) could be imaginary.
dLdb = zeros(size(log_cancer_rates)) + real(exp(log(kappaomegaPP-1)-log_cancer_rates));


log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + log_all_rates    , 1);
dLdb = dLdb +  exp(  log(tumordata(:,dd)') + log(model.theta(dd,end)) - log_P_t_given_theta     );

bb = expww ./sum(expww);
deriv_loglikelihood = real(exp(logsum(log(dLdb) + log(bb),2) + log(bb)));
deriv_loglikelihood = -deriv_loglikelihood + (dLdb.*bb);

%set the first derivative to zero to set the scale
deriv_loglikelihood(1) = 0;
deriv_loglikelihood = deriv_loglikelihood';

if (nnz(imag(deriv_loglikelihood)) > 0 | nnz(imag(loglikelihood)) > 0)
    error('imaginary number returned from rminimize in opt_beta_cancer');
end


%rasmussen's conjugate gradient method minimizes, so we negative the
%derivative and loglikelihood
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
