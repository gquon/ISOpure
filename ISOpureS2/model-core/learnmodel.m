function [ISOpureS2model loglikelihood] = learnmodel(tumordata, BB, ISOpureS1model, MIN_KAPPA)
%function [ISOpureS2model loglikelihood] = learnmodel(tumordata, BB, ISOpureS1model, MIN_KAPPA)
% INPUT:
%  tumordata: (same as for ISOpureS1) a GxD matrix representing gene
%  expression profiles of heterogeneous (mixed) tumor samples, where G is
%  the number of genes, D is the number of tumor samples.
%
%  BB: (same as for ISOpureS1) represents B = [b_1 ... b_K] matrix (from
%  Genome Medicine paper). a GxK matrix, where K is the number of normal profiles
%      (\beta_1,...,\beta_K in the paper), G is the number of genes.  These
%      are the normal profiles representing normal cells that contaminate
%      the tumor samples  (i.e. normal samples from the same tissue
%      location as the tumor).
%
%  ISOpureS1model: output model structure from ISOpureS1 code
%
%  MIN_KAPPA: (optional) The minimum value allowed for the strength parameters kappa_d placed
%  over the individual cancer profiles c_n (see Quon et al, 2013).  By default, this is set
%  to 1/min(m) (where m is the reference cancer profile) such that the log likelihood of the model is always finite.  However,
%  when the min(m) is very small, this forces MIN_KAPPA to be very large, and can sometimes
%  cause the cancer profiles to look too similar to the reference profile m
%  If this is the case, you can try setting MIN_KAPPA=1, or some other small value.  For reference, for the data
%  presented in Quon et al., 2013, MIN_KAPPA is on the order of 10^5 - 10^6.
%
%
% OUTPUT:
%
% loglikelihood: log likelihood of the final model
%
% ISOpureS2model: a structure with the following important fields:
%
%  theta: a DxK matrix, giving the fractional composition of each tumor
%  sample.  Each row represents a tumor sample that was part of the input,
%  and the first K-1 columns correspond to the fractional composition with
%  respect to the Source Panel contaminants.  The last column represents
%  the fractional composition of the pure cancer cells.  In other words,
%  each row sums to 1, and element (i,j) of the matrix denotes the fraction 
%  of tumor i attributable to component j (where the last column refers to 
%  cancer cells, and the first K-1 columns refer to different 'normal cell'
%  components).  The "% cancer", or tumor purity, estimate of each tumor is
%  simply the last column of theta. ffds
%
%  alphapurities:  (same as ISOpureS1) tumor purities (alpha_i in paper),
%  same as the last column of the theta variable, pulled out for user 
%  convenience. 
%
%  cc_cancerprofiles: purified cancer profiles.  This matrix is of the same
%  dimensionality as tumordata, and is also on the same scale (i.e.
%  although ISOpureS2 treats purified cancer profiles as parameters of a
%  multinomial distribution, we re-scale them to be on the same scale as
%  the input tumor profiles -- see Genome Medicine paper).  column ii of
%  cc_cancerprofiles corresponds to column ii of tumordata.



%make sure BB is a proper intensity/read count matrix (no negative elements)
if (min(min(BB))<0)
    error('negative elements found in input matrix BB.')
end

%make sure minimum value is not 0 (i.e. all genes need to have some probability of being observed in a sample)
if (min(min(BB))==0)
    nzix=find(BB>0);
    mymin=min(BB(nzix));
    BB(find(BB==0))=mymin;
    warning(['minimum element in input matrix BB is 0 -- setting all zeros to smallest non-zero element ' num2str(mymin)]);
end

%make sure data is not log transformed
if (max(max(BB)) < 30)
    warning('maximum element in matrix BB is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
end

%make sure tumor data is not log transformed
if (max(max(tumordata)) < 30)
warning('maximum element in matrix tumordata is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
end



%initial value of kappa of 10^4 for all tumors seems to work well.  This is optimized
%later.
kappa = 10^4;

%we work with transpose of BB
BBtranspose=BB';
%initial starting point of theta will be where we left off at in ISOpureS1
INITIAL_THETA=ISOpureS1model.theta;
%will start VV variable where we left off in ISOpureS1, except add some
%small random amount to each component, because many components may be
%exactly 1 (lower bound set on VV)
INITIAL_VV = ISOpureS1model.vv+rand(size(ISOpureS1model.vv));;
%prior on the tumor-specific cancer profiles is just the reference cancer
%profile learned in ISOpureS1
PPtranspose = exp(ISOpureS1model.log_all_rates(end,:));

NTOPICS=size(BBtranspose,1)+1;

%make sure PP and BB are scaled to sum to 1
PPtranspose = PPtranspose ./ repmat(sum(PPtranspose,2),1, size(PPtranspose,2));
BBtranspose = BBtranspose ./ repmat(sum(BBtranspose,2),1, size(BBtranspose,2));

disp('initializing..');

%identify the minimum value of kappa such that the Dirichlet prior over
%cancer profiles will give real-valued likelihoods

if (nargin < 4)
	MIN_KAPPA = 1/min(min(PPtranspose));
end
disp(['MIN_KAPPA set to ' num2str(MIN_KAPPA)]);

INIT_MODEL = new_model(tumordata, max(kappa, 10*MIN_KAPPA), INITIAL_VV, PPtranspose, BBtranspose);
INIT_MODEL.MIN_KAPPA = MIN_KAPPA;

%set initial theta to values learned in ISOpureS1
INIT_MODEL.theta = INITIAL_THETA;

%estimate parameters/hidden variables
[ISOpureS2model loglikelihood] = optmodel(tumordata, INIT_MODEL);

%copy over some important variables
ISOpureS2model.alphapurities=ISOpureS1model.alphapurities;

%transpose log_cc to get same dimensions as tumordata, then scale
%multinomial parameters to be on the same scale as tumordata
ISOpureS2model.cc_cancerprofiles=exp(ISOpureS2model.log_cc');
ISOpureS2model.cc_cancerprofiles=ISOpureS2model.cc_cancerprofiles .* repmat(sum(tumordata,1),size(ISOpureS2model.cc_cancerprofiles,1),1);
