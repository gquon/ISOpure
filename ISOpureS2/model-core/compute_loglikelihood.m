function loglikelihood = compute_loglikelihood(tumordata, model)
%function loglikelihood = compute_loglikelihood(tumordata, model)
% computes the full loglikelihood function that is being optimized in ISOpure
% step two

loglikelihood = 0;

%D is # tumor samples, W is # genes/transcripts
D = size(tumordata,2);
W = size(model.log_BBtranspose,2);

%add ln P(\theta_d|\vv)
loglikelihood = loglikelihood + D*(gammaln(sum(model.vv)) - sum(gammaln(model.vv))   ) + sum((model.vv-1) * log(model.theta)');

kappaomegaPP = model.omega * model.PPtranspose .* repmat(model.kappa',1, size(model.PPtranspose,2));

for dd=1:D
    log_all_rates = [model.log_BBtranspose; model.log_cc(dd,:)];
    %add ln P(t_d|\theta,BB,cc_i)
    log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + log_all_rates    , 1);
    loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));
    
    %add ln P(cc_i | \kappa, \omega, PP)
    loglikelihood = loglikelihood + gammaln(sum(kappaomegaPP(dd,:))) - sum(gammaln(kappaomegaPP(dd,:))) + ((kappaomegaPP(dd,:)-1) * model.log_cc(dd,:)');
end
