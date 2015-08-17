function [model loglikelihood] = opt_model(tumordata, model)
%function [model loglikelihood] = opt_model(tumordata, model)
% optimizes the ISOpure parameters cyclically until convergence

loglikelihood_old = -Inf;
change_ll_frac=Inf;

NUM_ITERATIONS_RMINIMIZE=20 ;
NUM_GRID_SEARCH_ITERATIONS=0; %used to do grid searching on parameters
NUM_ITERATIONS=35;

iter=1;

while (change_ll_frac > 0.0000001 && iter < 100)
    %optimize c_i variables (tumor-specific cancer profiles)
    disp('--- optimizing cc...');
	model = opt_cc(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

    %optimize mixing proportions theta
    disp('--- optimizing theta...');
	model = opt_theta(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

    %optimize priors over mixing proportions vv
    disp('--- optimizing vv...');
	model = opt_vv(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

    %optimize tumor-specific strength parameters kappa
    disp('--- optimizing kappa...');
	model = opt_kappa(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

    %the code to optimize omega is here (for future extension), but not
    %needed because there's only one component to the prior over the
    %tumor-specific cancer profiles (the reference cancer profile)
    %disp('--- optimizing omega...');
	%model = opt_omega(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

	loglikelihood = compute_loglikelihood(tumordata, model);
	change_ll = loglikelihood-loglikelihood_old;
    change_ll_frac = abs(change_ll/loglikelihood);

	disp(['iter: ' num2str(iter) '/' num2str(NUM_ITERATIONS) ', log likelihood: ' num2str(loglikelihood) ', frac change: ' num2str(change_ll_frac)]);


	loglikelihood_old = loglikelihood;
    iter=iter+1;
end
