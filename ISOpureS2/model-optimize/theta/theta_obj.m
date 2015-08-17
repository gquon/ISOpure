function [loglikelihood, deriv_loglikelihood] = theta_obj(ww, tumordata, dd, model)
%function [loglikelihood, deriv_loglikelihood] = theta_obj(ww, tumordata, dd, model)
% computes the part of the likelihood function relevant to optimizing
% theta (loglikelihood), as well as the derivative of the likelihood with respect to
% theta (deriv_loglikelihood)

%remaining = 1-model.theta(dd,end);  %when %cancer almost is 1, MATLAB will round remaining down to exactly 0.
remaining = sum(model.theta(dd,1:(end-1)));

expww = exp(ww');
theta = remaining .* expww ./ sum(expww);
alltheta = [theta model.theta(dd,end)];
loglikelihood = (model.vv-1) * log(alltheta');
W = size(model.log_BBtranspose,2);
K = size(model.theta,2);

log_all_rates = [model.log_BBtranspose; model.log_cc(dd,:)];

log_P_t_given_theta = logsum(   repmat( log(alltheta), W, 1)' + log_all_rates    , 1);
loglikelihood = loglikelihood + (log_P_t_given_theta * tumordata(:,dd));


dLdtheta = ((model.vv(1:(end-1))-1)./theta) +   sum(  exp(log_all_rates(1:(end-1),:) - repmat(log_P_t_given_theta, K-1,1))' .* repmat(tumordata(:,dd),1,K-1)         );

%change of variables
dthetadw = -(expww'*expww)./(sum(expww)^2);
dthetadw = dthetadw + diag(expww/sum(expww));
dthetadw = dthetadw.*remaining;

deriv_loglikelihood = (dLdtheta*dthetadw)';

%set the first derivative to be zero, to set the scale of the w's
deriv_loglikelihood(1) = 0;

%rasmussen's conjugate gradient method minimizes, so we negative the
%derivative and log likelihood
loglikelihood = -loglikelihood;
deriv_loglikelihood = -deriv_loglikelihood;
