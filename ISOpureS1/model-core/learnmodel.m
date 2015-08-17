function [ISOpureS1model loglikelihood] = learnmodel(tumordata, BB, PP, MIN_KAPPA)
%function [ISOpureS1model loglikelihood] = learnmodel(tumordata, BB, PP, MIN_KAPPA)
%
% INPUT:
%  tumordata: a GxD matrix representing gene expression profiles of
%  heterogeneous (mixed) tumor samples, where G is the number of genes, D
%  is the number of tumor samples.
%
%  BB: represents B = [b_1 ... b_K] matrix (from Genome Medicine paper)
%      a GxK matrix, where K is the number of normal profiles
%      (\beta_1,...,\beta_K in the paper), G is the number of genes.  These
%      are the normal profiles representing normal cells that contaminate
%      the tumor samples  (i.e. normal samples from the same tissue
%      location as the tumor).  The minimum element of BB must be greater than 0 --
%      i.e. every gene/transcript must be observed on some level in each normal sample.
%
%  PP: (optional) a GxM matrix, representing the expression profiles whose
%  convex combination form the prior over the purified cancer profile
%  learned.   If only primary tumors from the same site of origin are
%  represented in tumordata, then this is the same as BB (default
%  behavior).   This parameter is for backwards compatibility and replacing
%  the original ISOLATE code from the 2009 Bioinformatics paper, and can
%  represent potential sites of origins of the metastatic tumor (in which
%  case tumordata represents one or more expression profiles of the
%  secondary tumor).  Set PP=BB for default behavior, and if you need to specify
%  MIN_KAPPA.
%
%  MIN_KAPPA: (optional) The minimum value allowed for the strength parameter kappa' placed
%  over the reference cancer profile m (see Quon et al, 2013).  By default, this is set
%  to 1/min(BB), such that the log likelihood of the model is always finite.  However,
%  when the min(BB) is very small, this forces MIN_KAPPA to be very large, and can sometimes
%  cause the reference profile m to look too much like a 'normal profile' (and therefore
%  you may observe the tumor samples having low % cancer content estimates).  If this is the case,
%  you can try setting MIN_KAPPA=1, or some other small value.  For reference, for the data
%  presented in Quon et al., 2013, MIN_KAPPA is on the order of 10^5. 
%
% OUTPUT:
%
% loglikelihood: log likelihood of the final model
%
% ISOpureS1model: a structure with the following important fields:
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
%  simply the last column of theta.
%
%  alphapurities:  tumor purities (alpha_i in paper), same as the last
%  column of the theta variable, pulled out for user convenience.
%
%  mm: reference cancer profile, in the form of parameters of a multinomial
%  or discrete distribution (sum of elements is 1).  This is the same as
%  the purified cancer profile that ISOLATE was designed to learn.
%
%  omega: a Mx1 vector describing the convex combination weights learned by
%  ISOpure step 1 over the PPtranspose matrix, that when applied to the
%  Site of Origin Panel, forms the prior over the reference cancer profile.
%  When ISOpure step 1 is used in a similar fashion to the ISOLATE
%  algorithm, entry ii indicates the "probability" that the normal profile
%  in the (ii)th column of PP is the site of origin of the secondary tumors
%  stored in tumordata.  


%by default, we are looking at primary tumors from the same site, so the
%"Source Panel" (profiles that form the components of the prior over the reference cancer profile) is the same as
%BB
if nargin < 3
    PP=BB;
end

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

%do same checked for the PP matrix
if (min(min(PP))<0)
error('negative elements found in input matrix PP.')
end

if (min(min(PP))==0)
nzix=find(PP>0);
mymin=min(PP(nzix));
PP(find(PP==0))=mymin;
warning(['minimum element in input matrix PP is 0 -- setting all zeros to smallest non-zero element ' num2str(mymin)]);
end

if (max(max(PP)) < 30)
warning('maximum element in matrix PP is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
end


%make sure tumor data is not log transformed
if (max(max(tumordata)) < 30)
warning('maximum element in matrix tumordata is less than 30 -- make sure data is in normal (not log) space, otherwise output is wrong.');
end


%we work with these BB and PP in their transpose
BBtranspose=BB';
PPtranspose=PP';

%NTOPICS is the total number of component profiles (# normals + 1 for the
%reference cancer profile)
NTOPICS=size(BBtranspose,1)+1;

%initial value of kappa of 10^4 seems to work well.  This is optimized
%later.
kappa = 10^4;

%also randomly initialize prior over mixing proportions, theta
INITIAL_VV = rand(1,NTOPICS) + 1;  INITIAL_VV(end) = INITIAL_VV(end)+5;

%we standardize the normal profiles to sum to 1, so that we can interpret
%them as parameters of a discrete or multinomial distribution
PPtranspose = PPtranspose ./ repmat(sum(PPtranspose,2),1, size(PPtranspose,2));
BBtranspose = BBtranspose ./ repmat(sum(BBtranspose,2),1, size(BBtranspose,2));

loglikelihood = -Inf;

disp('-----------------');
disp(['Initializing...']);

%MIN_KAPPA represents the minimum value of kappa that we enforce during
%optimization, so that the Dirichlet distribution gives real valued
%loglikelihood values

if (nargin < 4)
	MIN_KAPPA = 1/min(min(PPtranspose));
end
disp(['MIN_KAPPA set to ' num2str(MIN_KAPPA)]);

%initialize the model structure that holds the parameters
INIT_MODEL = new_model(tumordata, max(kappa, 10*MIN_KAPPA), INITIAL_VV, PPtranspose, BBtranspose);
INIT_MODEL.MIN_KAPPA = MIN_KAPPA;

%optimize.
[ISOpureS1model newloglikelihood] = optmodel(tumordata, INIT_MODEL);

%explicitly save reference cancer profile mm, cancer purities alpha for
%ease of extraction by user
ISOpureS1model.mm = exp(ISOpureS1model.log_all_rates(end,:)');
ISOpureS1model.alphapurities=ISOpureS1model.theta(:,end);

loglikelihood = newloglikelihood;
