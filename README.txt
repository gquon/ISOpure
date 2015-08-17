ISOpure
===========

ISOpure is a two-step algorithm, as outlined in the Genome Medicine 2013
manuscript.   example.m gives an example how to run both steps.  Below is
a summary of the inputs and outputs of both steps -- the main function
being called is learnmodel.m, and is stored in each step's respective
directories  in model-core/learnmodel.m.

Note in general, that the variable names are matched to the Genome
Medicine paper notation as closely as possible.  Some variables, like c_i
and v_i, are actually named cc and vv in the code, so that they can be
searched for easily. 

Note on relationship to ISOLATE
===========
ISOpure step 1 (ISOpureS1) is functionally equivalent to the ISOLATE model
we published in 2009 in Bioinformatics, except it essentially performs
inference much more quickly, and uses a Dirichlet prior over the reference
(average purified) cancer profile mm, instead of a set of gamma priors
governing the 'perturbations' that cancer induces on normal cells.  We
recommend that ISOpureS1 be used instead of ISOLATE, as the input
requirements are the same but inference will be much quicker.

 
==================================================
1) ISOpure step 1 (ISOpureS1) 
==================================================
There's a lot of subdirectories (made for code organization), so when
you run ISOpureS1, you have to run a command like
"addpath(genpath('ISOpureS1'))" to make the entire directory
tree accessible for the code to run.

The main function to run is in learnmodel.m, and takes by default
2 matrices for input (3rd is optional):

 INPUT:
  tumordata: a GxD matrix representing gene expression profiles of
  heterogeneous (mixed) tumor samples, where G is the number of genes, D
  is the number of tumor samples.

  BB: represents B = [b_1 ... b_K] matrix (from Genome Medicine paper)
      a GxK matrix, where K is the number of normal profiles
      (\beta_1,...,\beta_K in the paper), G is the number of genes.  These
      are the normal profiles representing normal cells that contaminate
      the tumor samples  (i.e. normal samples from the same tissue
      location as the tumor).

  PP: (optional) a GxM matrix, representing the expression profiles whose
  convex combination form the prior over the purified cancer profile
  learned.   If only primary tumors from the same site of origin are
  represented in tumordata, then this is the same as BB (default
  behavior).   This parameter is for backwards compatibility and replacing
  the original ISOLATE code from the 2009 Bioinformatics paper, and can
  represent potential sites of origins of the metastatic tumor (in which
  case tumordata represents one or more expression profiles of the
  secondary tumor).

 OUTPUT:

 loglikelihood: log likelihood of the final model

 ISOpureS1model: a structure with the following important fields:

  theta: a DxK matrix, giving the fractional composition of each tumor
  sample.  Each row represents a tumor sample that was part of the input,
  and the first K-1 columns correspond to the fractional composition with
  respect to the Source Panel contaminants.  The last column represents
  the fractional composition of the pure cancer cells.  In other words,
  each row sums to 1, and element (i,j) of the matrix denotes the fraction
  of tumor i attributable to component j (where the last column refers to
  cancer cells, and the first K-1 columns refer to different 'normal cell'
  components).  The " cancer", or tumor purity, estimate of each tumor is
  simply the last column of theta.

  alphapurities:  tumor purities (alpha_i in paper), same as the last
  column of the theta variable, pulled out for user convenience.

  mm: reference cancer profile, in the form of parameters of a multinomial
  or discrete distribution (sum of elements is 1).  This is the same as
  the purified cancer profile that ISOLATE was designed to learn.

  omega: a Mx1 vector describing the convex combination weights learned by
  ISOpure step 1 over the PPtranspose matrix, that when applied to the
  Site of Origin Panel, forms the prior over the reference cancer profile.
  When ISOpure step 1 is used in a similar fashion to the ISOLATE
  algorithm, entry ii indicates the "probability" that the normal profile
  in the (ii)th column of PP is the site of origin of the secondary tumors
  stored in tumordata.  


==================================================
1) ISOpure step 2 (ISOpureS2) 
==================================================
Once again, there's a lot of subdirectories, so when
you run ISOpure step 2, you have to run a command like
"addpath(genpath('ISOpureS2'))" to make the entire directory
tree accessible.

The main function to run is in learnmodel.m, and should be
adequately documented on syntax.  ISOpureS2 needs the ISOpureS1model
structure from Step 1 to run properly.  Its input requirements are:

 INPUT:
  tumordata: (same as for ISOpureS1) a GxD matrix representing gene
  expression profiles of heterogeneous (mixed) tumor samples, where G is
  the number of genes, D is the number of tumor samples.

  BB: (same as for ISOpureS1) represents B = [b_1 ... b_K] matrix (from
  Genome Medicine paper). a GxK matrix, where K is the number of normal profiles
      (\beta_1,...,\beta_K in the paper), G is the number of genes.  These
      are the normal profiles representing normal cells that contaminate
      the tumor samples  (i.e. normal samples from the same tissue
      location as the tumor).

  ISOpureS1model: output model structure from ISOpureS1 code

 OUTPUT:

 loglikelihood: log likelihood of the final model

 ISOpureS2model: a structure with the following important fields:

  theta: a DxK matrix, giving the fractional composition of each tumor
  sample.  Each row represents a tumor sample that was part of the input,
  and the first K-1 columns correspond to the fractional composition with
  respect to the Source Panel contaminants.  The last column represents
  the fractional composition of the pure cancer cells.  In other words,
  each row sums to 1, and element (i,j) of the matrix denotes the fraction 
  of tumor i attributable to component j (where the last column refers to 
  cancer cells, and the first K-1 columns refer to different 'normal cell'
  components).  The " cancer", or tumor purity, estimate of each tumor is
  simply the last column of theta. ffds

  alphapurities:  (same as ISOpureS1) tumor purities (alpha_i in paper),
  same as the last column of the theta variable, pulled out for user 
  convenience. 

  cc_cancerprofiles: purified cancer profiles.  This matrix is of the same
  dimensionality as tumordata, and is also on the same scale (i.e.
  although ISOpureS2 treats purified cancer profiles as parameters of a
  multinomial distribution, we re-scale them to be on the same scale as
  the input tumor profiles -- see Genome Medicine paper).  column ii of
  cc_cancerprofiles corresponds to column ii of tumordata.

Sample Scripts
=================
For convenience we have included  a sample script running ISOpure on the
Beer et al., 2002, lung adenocarcinoma data.  We pre-processed the data
with the RMA algorithm (see Genome Medicine paper) and converted back
into microarray intensity space afterwards. 
