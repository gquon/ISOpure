function model = opt_omega(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_omega(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
% optimize omega, the convex mixing weights that govern prior over the
% reference cancer profile.

%P is the number of normals in this "Site of origin panel"
P = size(model.omega,1);

%because omega variables are constrained (must be parameters of multinomial/discrete
%distribution), we don't directly optimize the likelihood function w.r.t.
%omega, but we perform change of variables to do unconstrained
%optimization.  We therefore store these unconstrained variables in the
%field "omega_weights", and update these variables
if ~isfield(model, 'omega_weights')
        model.omega_weights = log(model.omega);
end


init_xx = model.omega_weights;

[xx tmp num_iter] = rminimize(init_xx, 'omega_obj', NUM_ITERATIONS_RMINIMIZE,tumordata, model);

%convert from unconstrained variables to omega
model.omega = exp(xx)/sum(exp(xx));
model.omega_weights = xx;

loglikelihood = omega_compute_loglikelihood(model.omega, tumordata, model);


%try a few other random initializations of omega to see if we can get a
%better loglikelihood value -- not that expensive an operation.
if (iter <= NUM_GRID_SEARCH_ITERATIONS)

        for i=0:1:10
        init_xx = rand(P,1); 

        [xx tmp num_iter] = rminimize(log(init_xx), 'omega_obj', NUM_ITERATIONS_RMINIMIZE, tumordata, model);
        newomega = exp(xx) ./ sum(exp(xx));
        newll = omega_compute_loglikelihood(newomega, tumordata, model);
        
        if (newll > loglikelihood)
            model.omega_weights = xx;
            model.omega = newomega;
            loglikelihood = newll;
        end
        
        end
end

end
