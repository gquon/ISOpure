function model = opt_kappa(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_kappa(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)

%initialize vv to what we learned it on previous step
kappa = model.kappa;

%we're actually optimizing log(kappa), to keep it >=model.MIN_KAPPA
init_xx = log(kappa - model.MIN_KAPPA)';

% optimize
[xx tmp num_iter] = rminimize(init_xx, 'kappa_obj', NUM_ITERATIONS_RMINIMIZE,model);

%re-exponentiate
model.kappa = exp(xx') + model.MIN_KAPPA;