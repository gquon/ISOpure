% This is a sample script to run ISOLATE.

%load lung adenocarcinoma expression data from Beer et al., 2002.
load beer02natmed.mat;

%add path to ISOpureS1 code so we can run all subroutines
addpath(genpath('ISOpureS1'));



disp('---------------------');
disp('running ISOpure step 1...');
%run ISOpure S1
[ISOpureS1model loglikelihood] = learnmodel(tumordata, normaldata);

%now, remove path to ISOpureS1 code.  We need to do this because ISOpureS2 has some functions with the same name as those in ISOpureS1 (but different code).
rmpath(genpath('ISOpureS1'));



%add path to ISOpureS2 code so we can run all subroutines
addpath(genpath('ISOpureS2'));

disp('---------------------');
disp('running ISOpure step 2...');
%run ISOpure S2
[ISOpureS2model loglikelihood] = learnmodel(tumordata, normaldata, ISOpureS1model);

%save output
save ISOpure-beer.mat ISOpureS1model ISOpureS2model;

%pull out purified cancer profiles
cancerdata = ISOpureS2model.cc_cancerprofiles;

%now, use cancerdata in the same pipelines as we would have used tumordata.
