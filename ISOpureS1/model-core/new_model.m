function model = new_model(tumordata, kappa, INITIAL_VV, PPtranspose, BBtranspose)
%function model = new_model(tumordata, kappa, INITIAL_VV,PPtranspose, BBtranspose)
% Generates a new model structure to hold parameters

model.vv = INITIAL_VV;

%log_BBtranspose = log(B^T) = log( [b_1 ... b_K]^T )
model.log_BBtranspose = log(BBtranspose); %log BBtranspose
model.PPtranspose = PPtranspose;

model.kappa = kappa;

%D = # tumor samples
%K = # components of each tumor sample (# normal profiles + 1 cancer
%profile)
D = size(tumordata,2);
K = length(model.vv);

%randomly initialize theta but give higher weight to cancer component
%initially
model.theta = zeros(D,K);
model.theta(:,end) = 0.5;
model.theta(:,1:(end-1)) = 0.5/(K-1);

model.theta = model.theta ./ repmat(sum(model.theta,2), 1, K);

P = size(PPtranspose, 1);
model.omega = zeros(P, 1);
model.omega(:,:) = 1/P;

%model.log_all_rates = [B mm]^T , where B is the matrix of normal profiles
%[b_1 ... b_K] , and mm is the reference cancer profile
model.log_all_rates = [model.log_BBtranspose; log(model.omega'* model.PPtranspose)];
