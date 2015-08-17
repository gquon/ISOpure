function model = new_model(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose)
%function model = new_model(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose)
% Generates a new model structure to hold parameters

model.vv = INITIAL_VV;
model.log_BBtranspose = log(BBtranspose);
model.PPtranspose = PPtranspose;

D = size(tumordata,2);
K = length(model.vv);

model.kappa = ones(1,D) .* kappa;

P = size(PPtranspose, 1);
model.omega = zeros(D, P);
model.omega(:,:) = 1/P;

%initialize tumor-specific cancer profiles to the reference cancer profile
model.log_cc = model.omega*PPtranspose;
