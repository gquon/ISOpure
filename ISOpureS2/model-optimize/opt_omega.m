function model = opt_omega(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_omega(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
% not really used right now, but could in the future if we implement a
% mixture model for multiple reference cancer profiles

P = size(model.omega,2);
D = size(model.omega,1);

%because omega is constrained (parameters of multinomial/discrete
%distribution), we don't directly optimize the likelihood function w.r.t.
%omega, but we perform change of variables to do unconstrained
%optimization.  We therefore store these unconstrained variables in the
%field "omega_weights", and update these variables
if ~isfield(model, 'omega_weights')
    model.omega_weights = log(model.omega);
end

for dd = 1:D
    init_xx = model.omega_weights(dd,:)';
    
    %initial starting point
    [xx tmp num_iter] = rminimize(init_xx, 'omega_obj', NUM_ITERATIONS_RMINIMIZE, tumordata, dd, model);
    
    
    %re-exponentiate
    model.omega(dd,:) = exp(xx)'/sum(exp(xx));
    model.omega_weights(dd,:) = xx';
end

end
