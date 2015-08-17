function model = opt_cc(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_cc(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
% optimize the tumor-specific cancer profiles

D = size(model.log_cc,1);
W = size(model.log_BBtranspose,2);
K = size(model.theta,2);

%because cc is constrained (each cc_i are parameters of multinomial/discrete
%distribution), we don't directly optimize the likelihood function w.r.t.
%cc, but we perform change of variables to do unconstrained
%optimization.  We therefore store these unconstrained variables in the
%field "cc_weights", and update these variables
if ~isfield(model, 'cc_weights')
    model.cc_weights = model.log_cc;
end

for dd=1:D
    init_xx = model.cc_weights(dd,:)';
    
    %optimize
    [xx tmp num_iter] = rminimize(init_xx, 'cc_obj', NUM_ITERATIONS_RMINIMIZE,tumordata, dd, model);
    
    %re-exponentiate to get cc multinomial parameters
    model.cc_weights(dd,:) = xx';
    model.log_cc(dd,:) = xx'-logsum(xx,1);
end
