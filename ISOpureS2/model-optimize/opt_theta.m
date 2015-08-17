function model = opt_theta(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_theta(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)


%rminimize seems to only work on vectors - they are independent anyways

K = size(model.theta,2);
D = size(model.theta,1);

%because theta are constrained (each theta_d are parameters of multinomial/discrete
%distribution), we don't directly optimize the likelihood function w.r.t.
%theta, but we perform change of variables to do unconstrained
%optimization.  We therefore store these unconstrained variables in the
%field "theta_weights", and update these variables.
%Furthermore ,note that we are fixing the tumor purities (last column of
%theta), so we are only storing/updating the remaining columns of theta,
%and optimizing them to sum to 1-alpha_i for tumor i
if ~isfield(model, 'theta_weights')
    model.theta_weights = log(model.theta(:,1:(end-1)));
end


for dd=1:D
    ww = model.theta_weights(dd,:)';
    %remaining = 1-model.theta(dd,end); %numerically unstable when % cancer close to 1
    remaining = sum(model.theta(dd,1:(end-1)));
    %initial starting point
    init_xx = [ww];
    [xx tmp num_iter] = rminimize(init_xx, 'theta_obj', NUM_ITERATIONS_RMINIMIZE,tumordata, dd, model);
    
    %re-exponentiate
    model.theta(dd,1:(end-1)) = remaining .* exp(xx)'./sum(exp(xx));
    model.theta_weights(dd,:) = xx';    
end
