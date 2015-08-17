function model = opt_vv(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_vv(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)

%make sure vv is of the right dimensions
vv = reshape(model.vv, length(model.vv), 1);

%optimizing log(vv) to keep the components >=1
init_ww = [log(vv-1)];

sum_log_theta = sum(log(model.theta),1);

%optimize
ww = rminimize(init_ww, 'vv_obj', NUM_ITERATIONS_RMINIMIZE, sum_log_theta, size(tumordata,2));

%re-exponentiate
vv = exp(ww)+1;

%make sure vv is a column vector
model.vv = reshape(vv, 1, length(vv));