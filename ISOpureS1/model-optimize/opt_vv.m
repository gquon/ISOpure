function model = opt_vv(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%function model = opt_vv(tumordata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)
%
% optimize vv, the prior over the mixing proportions theta
%

%make sure these are of specific dimensions
vv = reshape(model.vv, length(model.vv), 1);



%note that we don't directly optimize vvs because it has constraints
%(must be >=1 to guarantee real-valued likelihoods).  
init_ww = [log(vv-1)];

sum_log_theta = sum(log(model.theta),1);
 
ww = rminimize(init_ww, 'vv_obj', NUM_ITERATIONS_RMINIMIZE, sum_log_theta, size(tumordata,2));

%convert back into vv 
vv = exp(ww)+1;
model.vv = reshape(vv, 1, length(vv));

%do some re-starts because they're not computationally expensive
if (iter <= NUM_GRID_SEARCH_ITERATIONS)
        loglikelihood = vv_compute_loglikelihood(model.vv, sum_log_theta, size(tumordata,2));

        for cancer_vv=10:10:100
                newvv = ones(size(model.vv'));  newvv(end) = cancer_vv;

                [xx tmp num_iter] = rminimize(log(newvv-1), 'vv_obj', NUM_ITERATIONS_RMINIMIZE, sum_log_theta, size(tumordata,2));
                newvv = exp(xx')+1;

                newll = vv_compute_loglikelihood(newvv, sum_log_theta, size(tumordata,2));
                
		if (newll > loglikelihood)
                        model.vv = newvv;
                        loglikelihood = newll;
                end
        end
end

