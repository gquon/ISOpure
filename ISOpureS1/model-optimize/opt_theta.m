function model = opt_theta(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_theta(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
% optimize thetas


K = size(model.theta,2);
D = size(model.theta,1);

%because thetas are constrained (must be parameters of multinomial/discrete
%distribution), we don't directly optimize the likelihood function w.r.t.
%theta, but we perform change of variables to do unconstrained
%optimization.  We therefore store these unconstrained variables in the
%field "theta_weights", and update these variables
if ~isfield(model, 'theta_weights')
        model.theta_weights = log(model.theta);
end


%update each theta_d separately
for dd=1:D
	init_xx = model.theta_weights(dd,:)';
	[xx tmp num_iter] = rminimize(init_xx, 'theta_obj', NUM_ITERATIONS_RMINIMIZE,tumordata, dd, model);

    %convert from unconstrained variables to theta
	model.theta(dd,:) = exp(xx)'/sum(exp(xx));
	model.theta_weights(dd,:) = xx';
end

end
