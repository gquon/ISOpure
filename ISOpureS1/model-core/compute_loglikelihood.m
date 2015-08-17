function loglikelihood = compute_loglikelihood(tumordata, model)
%function loglikelihood = compute_loglikelihood(tumordata, model)
%
% computes the full loglikelihood function that is being optimized in ISOpure
% step one

loglikelihood = 0;

%D is # tumor samples
D = size(tumordata,2);
%W is # transcripts/genes
W = size(model.log_all_rates,2);
%kappaomegaPP is the set of parameters of the Dirichlet prior over the
%reference cancer profile
kappaomegaPP = model.omega' * model.PPtranspose .* model.kappa;

%loglikelihood of reference cancer profile
loglikelihood = loglikelihood + gammaln(sum(kappaomegaPP)) - sum(gammaln(kappaomegaPP)) + ((kappaomegaPP-1) * model.log_all_rates(end,:)');

%loglikelihood of thetas
loglikelihood = loglikelihood + D*(gammaln(sum(model.vv)) - sum(gammaln(model.vv))   ) + sum((model.vv-1) * log(model.theta)');

for dd=1:D
	%loglikelihood of observed tumor profiles t_i
    log_ptgt = logsum(   repmat( log(model.theta(dd,:)), W, 1)' + model.log_all_rates    , 1);
	loglikelihood = loglikelihood + (log_ptgt * tumordata(:,dd));
end
