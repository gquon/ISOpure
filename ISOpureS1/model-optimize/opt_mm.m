function model = opt_mm(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_mm(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
% optimize the reference cancer profile mm
%

%because mm is constrained (must be parameters of multinomial/discrete
%distribution), we don't directly optimize the likelihood function w.r.t.
%mm, but we perform change of variables to do unconstrained
%optimization.  We therefore store these unconstrained variables in the
%field "mm_weights", and update these variables
if ~isfield(model, 'mm_weights')
        model.mm_weights = model.log_all_rates(end,:);
end

init_xx = model.mm_weights';


[xx tmp num_iter] = rminimize(init_xx, 'mm_obj', NUM_ITERATIONS_RMINIMIZE,tumordata, model);

%convert optimized values back to mm
model.mm_weights = xx';
model.log_all_rates(end,:) = xx'-logsum(xx,1);
