function model = opt_kappa(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_kappa(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
%
% optimize kappa, the strength parameter in the prior over the reference
% cancer profile

kappa = model.kappa;

%note that we don't directly optimize kappa because it has constraints
%(must be greater than the minimum determined in directoptlearn).  
init_xx = log(kappa - model.MIN_KAPPA)';

[xx tmp num_iter] = rminimize(init_xx, 'kappa_obj', NUM_ITERATIONS_RMINIMIZE,tumordata, model);

%convert optimized value back into kappa
model.kappa = exp(xx') + model.MIN_KAPPA;


%do some random restarts for kappa optimization -- this isn't an expensive
%operation.

if (iter <= NUM_GRID_SEARCH_ITERATIONS)
	loglikelihood = kappa_compute_loglikelihood(model.kappa, tumordata, model);
	
	MIN_POW_KAPPA = ceil(log10(model.MIN_KAPPA));

	for pow=MIN_POW_KAPPA:15
		scales = 10^pow .* ones(size(model.kappa'));
		init_xx = log(scales - model.MIN_KAPPA);

		[xx tmp num_iter] = rminimize(init_xx, 'kappa_obj', NUM_ITERATIONS_RMINIMIZE, tumordata, model);
		newkappa = exp(xx) + model.MIN_KAPPA;

		newll = kappa_compute_loglikelihood(newkappa, tumordata, model);
		
		if (newll > loglikelihood)
			model.kappa = newkappa';
			loglikelihood = newll;
		end	
	end
end
