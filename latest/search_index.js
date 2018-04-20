var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: JWAS)JWAS is a well-documented software platform based on Julia and an interactive Jupyter notebook for analyses of general univariate and multivariate Bayesian mixed effects models.  These models are especially useful for, but not limited to, routine single-trait and multi-trait genomic prediction and genome-wide association studies using either complete or incomplete genomic data (\"single-step\" methods). Currently, JWAS provides broad scope of analyses, e.g., a wide collection of Bayesian methods for whole-genome analyses, including shrinkage estimation and variable selection methods. The features of JWAS include:No limitations on fixed effects (e.g. herd-year, age, sex)                                                                    \nRandom effects other than markers (e.g. litter, pen)                                  \nRandom effects using pedigree information                                                                                \nRandom permanent environmental effects  \nSingle-trait analyses                                            \nMulti-trait analyses                                                                  \nUse of genomic information                                                                                \nComplete genomic data                                      		\nIncomplete genomic data\nCorrelated residuals		"
},

{
    "location": "theory/theory.html#",
    "page": "Some Theory",
    "title": "Some Theory",
    "category": "page",
    "text": ""
},

{
    "location": "theory/theory.html#Some-Theory-in-JWAS-1",
    "page": "Some Theory",
    "title": "Some Theory in JWAS",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory.html#A-Table-for-Bayesian-Linear-Mixed-Models-(BLMM)-1",
    "page": "Some Theory",
    "title": "A Table for Bayesian Linear Mixed Models (BLMM)",
    "category": "section",
    "text": "(Image: BLMM)"
},

{
    "location": "theory/theory.html#Models-1",
    "page": "Some Theory",
    "title": "Models",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory.html#Complete-Genomic-Data-1",
    "page": "Some Theory",
    "title": "Complete Genomic Data",
    "category": "section",
    "text": "The general form of the multivariate (univariate) mixed effects model for individual i from n individuals with complete genomic data in JWAS is\nmathbfy_i\n =sum_j=1^p_betaX_ijboldsymbolbeta_j+sum_k=1^p_uZ_ikmathbfu_k\n +sum_l=1^pM_ilboldsymbolalpha_l+mathbfe_i(1)where mathbfy_i is a vector of phenotypes of t traits for individual i; X_ij is the incidence matrix covariate corresponding to the jth fixed effect for individual i; boldsymbolbeta_j is a vector of jth fixed effects for the t traits; Z_ik is the incidence matrix covariate corresponding to the kth random effect for individual i; boldsymbolu_k is a vector of the kth random effects of t traits; M_il is the genotype covariate at locus l for individual i, p is the number of genotyped loci (each coded as 0,1,2), boldsymbolalpha_l is a vector of allele substitution effects or marker effects of t traits for locus j, and mathbfe_i is the vector of random residual effects of t traits for individual i. The JWAS implementation of this model involves missing phenotypes being imputed at each iteration of MCMC \\cite{sorensenGianolaBook} so that all individuals have observations for all traits. Note that when the number of traits t=1, the general form above simplifies to the single-trait  mixed effects model, and all vectors of effects in equation (1) become scalars."
},

{
    "location": "theory/theory.html#Incomplete-Genomic-Data-1",
    "page": "Some Theory",
    "title": "Incomplete Genomic Data",
    "category": "section",
    "text": "The general form of the multivariate (univariate) mixed effects model with incomplete genomic data (\"single-step\" methods) for non-genotyped individuals ismathbfy_i\n=sum_j=1^p_betaX_ijboldsymbolbeta_j+sum_k=1^p_uZ_ikmathbfu_k+\nsum_l=1^phatM_ilboldsymbolalpha_l+sum_m=1^p_epsilonZ_nimboldsymbolepsilon_m+boldsymbole_i (2)where mathbfy_i is a vector of phenotypes of t traits for non-genotyped individual i;  hatM_il is the imputed genotype covariate at locus l for non-genotyped individual i, Z_nim is the incidence matrix covariate corresponding to the mth imputation residual for individual i and boldsymbolepsilon_i is a vector of imputation residuals. W_im is the incidence matrix covariate corresponding to the mth random effect for individual i. That vector of imputation residuals, boldsymbolepsilon=beginbmatrixboldsymbolepsilon_1^T  boldsymbolepsilon_2^T  ldots  endbmatrix^T, are a priori assumed to be Nleft(0(mathbfA_nn-mathbfA_ngmathbfA_gg^-1mathbfA_gn)otimesmathbfG_gright), where mathbfA_nn is the partition of the numerator relationship matrix  mathbfA that corresponds to non-genotyped individuals, mathbfA_ng or its transpose mathbfA_gn are partitions of mathbfA corresponding to relationships between non-genotyped and genotyped individuals or vice versa,  mathbfA_gg is the  partition of mathbfA that corresponds to genotyped animals, and mathbfG_g is the additive genetic covariance matrix. All the other variables are the same as in equation (1)."
},

{
    "location": "theory/theory.html#Priors-1",
    "page": "Some Theory",
    "title": "Priors",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory.html#Priors-for-effects-other-than-markers-1",
    "page": "Some Theory",
    "title": "Priors for effects other than markers",
    "category": "section",
    "text": "The fixed effects are assigned flat priors. The vector of random effects, mathbfu=beginbmatrixmathbfu_1^T  mathbfu_2^T  ldots  mathbfu_p_2^Tendbmatrix^T, are a priori assumed to be Nleft(0mathbfAotimesmathbfGright) with various options for mathbfA. For example, mathbfA could be an identity matrix if boldsymbolu_k is assumed to be independently and identically distributed. mathbfA can be  the numerator relationship matrix, when boldsymbolu is a vector of polygenic effects and mathbfG represents the additive-genetic variance not explained by molecular markers. Note that boldsymbolu can also be a concatenation of vectors of different types of random effects, such as litter, pen, polygenic and maternal effects. The vector boldsymbole_i of residuals are a priori assumed to be independently and identically following multivariate normal distributions with null mean and covariance matrix mathbfR, which in turn is a priori assumed to have an inverse Wishart prior distribution, W_t^-1left(mathbfS_enu_eright). Note that when number of traits t=1, the priors for mathbfG and mathbfR in single-trait analyses follow scaled inverted chi-square distributions."
},

{
    "location": "theory/theory.html#Priors-for-marker-effects-1",
    "page": "Some Theory",
    "title": "Priors for marker effects",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory.html#single-trait-BayesA-1",
    "page": "Some Theory",
    "title": "single-trait BayesA",
    "category": "section",
    "text": "The prior assumption is that marker effects have identical and independent univariate-t distributions each with a null mean, scale parameter S^2_alpha and nu degrees of freedom. This is equivalent to assuming that the marker effect at locus i has a univariate normal with null mean and unknown, locus-specific variance sigma^2_i, which in turn is assigned a scaled inverse chi-square prior with scale parameter S^2_alpha and nu_alpha degrees of freedom."
},

{
    "location": "theory/theory.html#single-trait-BayesB-1",
    "page": "Some Theory",
    "title": "single-trait BayesB",
    "category": "section",
    "text": "In BayesB, the prior assumption is that marker effects have identical and independent mixture distributions, where each has a point mass at zero with probability pi and a univariate-t distribution with probability 1-pi having a null mean, scale parameter S^2_alpha and nu degrees of freedom. Thus, BayesA is a special case of BayesB with pi=0. Further, as in BayesA, the t-distribution in BayesB is equivalent to a univariate normal with null mean and unknown, locus-specific variance, which in turn is assigned a scaled inverse chi-square prior with scale parameter S^2_alpha and nu_alpha degrees of freedom. (A fast and efficient Gibbs sampler was implemented for BayesB in JWAS.)"
},

{
    "location": "theory/theory.html#single-trait-BayesC-and-BayesC\\pi-1",
    "page": "Some Theory",
    "title": "single-trait BayesC and BayesCpi",
    "category": "section",
    "text": "In BayesC, the prior assumption is that marker effects have identical and independent mixture distributions, where each has a point mass at zero with probability pi and a univariate-normal distribution with probability 1-pi having a null mean and variance sigma^2_alpha, which in turn has a scaled inverse chi-square prior with scale parameter S^2_alpha and nu_alpha degrees of freedom. In addition to the above assumptions, in BayesC pi, pi is treated as unknown with a uniform prior."
},

{
    "location": "theory/theory.html#multiple-trait-BayesABC-1",
    "page": "Some Theory",
    "title": "multiple-trait BayesABC",
    "category": "section",
    "text": "In multi-trait BayesCPi, the prior for alpha_lk, the marker effect of trait k for locus l, is a mixture with a point mass at zero and a univariate normal distribution conditional on sigma_k^2:beginalign*\nalpha_lkmidpi_ksigma_k^2  begincases\nsim Nleft(0sigma_k^2right)  probability(1-pi_k)\n0  probabilitypi_k\nendcases\nendalign*and the covariance between effects for traits k and k at the same locus, i.e., alpha_lk and alpha_lk^ isbeginalign*\ncovleft(alpha_lkalpha_lk^midsigma_kk^right)=begincases\nsigma_kk^  ifbothalpha_lkneq0andalpha_lk^neq0\n0  otherwise\nendcases\nendalign*The vector of marker effects at a particular locus boldsymbolalpha_l is written as boldsymbolalpha_l=boldsymbolD_lboldsymbolbeta_l, where boldsymbolD_l is a diagonal matrix with elements diagleft(boldsymbolD_lright)=boldsymboldelta_l=left(delta_l1delta_l2delta_l3ldotsdelta_ltright), where delta_lk is an indicator variable indicating whether the marker effect of locus l for trait k is zero or non-zero, and the vector boldsymbolbeta_l follows a multivariate normal distribution with null mean and covariance matrix boldsymbolG. The covariance matrix boldsymbolG is a priori assumed to follow an inverse Wishart distribution, W_t^-1left(mathbfS_betanu_betaright).In the most general case, any marker effect might be zero for any possible combination of t traits resulting in 2^t possible combinations of boldsymboldelta_l. For example, in a t=2 trait model, there are 2^t=4 combinations for  boldsymboldelta_l: (00), (01), (10), (11). Suppose in general we use numerical labels \"1\", \"2\",ldots, \"l\" for the 2^t possible outcomes for  boldsymboldelta_l, then the prior for  boldsymboldelta_l is a categorical distributionbeginalign*\n  pleft(boldsymboldelta_l=iright)\n=  Pi_1Ileft(boldsymboldelta_l=1right)+Pi_2Ileft(boldsymboldelta_l=2right)++Pi_lIleft(boldsymboldelta_l=lright)\nendalign*where sum_i=1^lPi_i=1 with Pi_i being the prior probability that the vector boldsymboldelta_l corresponds to the vector labelled i. A Dirichlet distribution with all parameters equal to one, i.e., a uniform distribution, can be used for the prior for boldsymbolPi=left(Pi_1Pi_2Pi_lright).   The differences in multi-trait BayesB method is that the prior for boldsymbolbeta_l is a multivariate t distribution, rather than a multivariate normal distribution. This is equivalent to assuming boldsymbolbeta_l has a multivariate normal distribution with null mean and locus-specific covariance matrix boldsymbolG_l, which is assigned an inverse Wishart prior, W_t^-1left(mathbfS_betanu_betaright). Multi-trait BayesA method is a special case of multi-trait BayesB method where boldsymboldelta_l is always a vector of ones.referencesMeuwissen T, Hayes B, Goddard M. Prediction of total genetic value using genome-wide dense marker maps. Genetics, 157:1819–1829, 2001.\nFernando R, Garrick D. Bayesian methods applied to GWAS. Methods Mol Biol. 2013;1019:237–274.\nCheng H, Garrick D, Fernando R. A fast and efficient Gibbs sampler for BayesB in whole- genome analyses. Genet Sel Evol, 2015, 47:80.\nFernando R, Dekkers J,Garrick D. A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses. Genetics Selection Evolution, 2015 46(1), 50.\nCheng H, Kizilkaya K, Zeng J, Garrick D, Fernando R. Genomic Prediction from Multiple-trait Bayesian Regression Methods using Mixture Priors. Genetics. 2018"
},

{
    "location": "man/guide.html#",
    "page": "Guide",
    "title": "Guide",
    "category": "page",
    "text": ""
},

{
    "location": "man/guide.html#guide-1",
    "page": "Guide",
    "title": "guide",
    "category": "section",
    "text": "abc"
},

{
    "location": "man/examples.html#",
    "page": "example",
    "title": "example",
    "category": "page",
    "text": ""
},

{
    "location": "man/examples.html#example-1",
    "page": "example",
    "title": "example",
    "category": "section",
    "text": ""
},

{
    "location": "man/examples.html#Get-Started-1",
    "page": "example",
    "title": "Get Started",
    "category": "section",
    "text": ""
},

{
    "location": "man/examples.html#JWAS.runMCMC",
    "page": "example",
    "title": "JWAS.runMCMC",
    "category": "function",
    "text": "runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,starting_value=false,printout_frequency=100,missing_phenotypes=false,constraint=false,methods=\"conventional (no markers)\",output_samples_frequency::Int64 = 0)\n\nRun MCMC (marker information included or not) with sampling of variance components.\n\navailable methods include \"conventional (no markers)\", \"BayesC0\", \"BayesC\", \"BayesCC\",\"BayesB\".\nmissing_phenotypes\nPi for single-trait analyses is a number; Pi for multi-trait analyses is a dictionary such as Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1),\nif Pi (Π) is not provided in multi-trait analysis, it will be generated assuming all markers have effects on all traits.\nsave MCMC samples every output_samples_frequency iterations\nstarting_value can be provided as a vector for all location parameteres except marker effects.\nprint out the monte carlo mean in REPL with printout_frequency\nconstraint=true if constrain residual covariances between traits to be zero.\n\n\n\n"
},

{
    "location": "man/examples.html#example-1-1",
    "page": "example",
    "title": "example 1",
    "category": "section",
    "text": "runMCMClink to JWAS.jl Documentation\nlink to add_genotypes"
},

{
    "location": "man/examples.html#JWAS.add_genotypes",
    "page": "example",
    "title": "JWAS.add_genotypes",
    "category": "function",
    "text": "add_genotypes(mme::MME,file,G;separator=\' \',header=true,center=true,G_is_marker_variance=false,df=4.0)\n\nGet marker informtion from a genotype file (same order as the phenotype file).\nG defaults to the genetic variance with degree of freedom df=4.0.\nFile format:\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,2,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n"
},

{
    "location": "man/examples.html#JWAS.add_markers",
    "page": "example",
    "title": "JWAS.add_markers",
    "category": "function",
    "text": "same to add_genotypes\n\n\n\n"
},

{
    "location": "man/examples.html#example-2-1",
    "page": "example",
    "title": "example 2",
    "category": "section",
    "text": "add_genotypes\nadd_markers"
},

{
    "location": "man/examples.html#Test-these-examples-1",
    "page": "example",
    "title": "Test these examples",
    "category": "section",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#Bayesian-Linear-Mixed-Models-(Genomic-Data)-1",
    "page": "Contributing",
    "title": "Bayesian Linear Mixed Models (Genomic Data)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#Univariate-Linear-Mixed-Model-(Genomic-data)-1",
    "page": "Contributing",
    "title": "Univariate Linear Mixed Model (Genomic data)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#Multivariate-Linear-Mixed-Model-(Genomic-data)-1",
    "page": "Contributing",
    "title": "Multivariate Linear Mixed Model (Genomic data)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/conventionalBLMM.html#",
    "page": "Linear Mixed Model (conventional)",
    "title": "Linear Mixed Model (conventional)",
    "category": "page",
    "text": ""
},

{
    "location": "examples/conventionalBLMM.html#Bayesian-Linear-Mixed-Models-1",
    "page": "Linear Mixed Model (conventional)",
    "title": "Bayesian Linear Mixed Models",
    "category": "section",
    "text": ""
},

{
    "location": "examples/conventionalBLMM.html#Univariate-Linear-Mixed-Model-(conventional)-1",
    "page": "Linear Mixed Model (conventional)",
    "title": "Univariate Linear Mixed Model (conventional)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/conventionalBLMM.html#Multivariate-Linear-Mixed-Model-(conventional)-1",
    "page": "Linear Mixed Model (conventional)",
    "title": "Multivariate Linear Mixed Model (conventional)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/LinearAdditiveGeneticModel.html#",
    "page": "Linear Additive Genetic Model",
    "title": "Linear Additive Genetic Model",
    "category": "page",
    "text": ""
},

{
    "location": "examples/LinearAdditiveGeneticModel.html#Bayesian-Linear-Additive-Genetic-Model-1",
    "page": "Linear Additive Genetic Model",
    "title": "Bayesian Linear Additive Genetic Model",
    "category": "section",
    "text": ""
},

{
    "location": "examples/LinearAdditiveGeneticModel.html#Univariate-Linear-Additive-Genetic-Model-1",
    "page": "Linear Additive Genetic Model",
    "title": "Univariate Linear Additive Genetic Model",
    "category": "section",
    "text": ""
},

{
    "location": "examples/LinearAdditiveGeneticModel.html#Multivariate-Linear-Additive-Genetic-Model-1",
    "page": "Linear Additive Genetic Model",
    "title": "Multivariate Linear Additive Genetic Model",
    "category": "section",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#",
    "page": "Linear Mixed Model (Genomic data)",
    "title": "Linear Mixed Model (Genomic data)",
    "category": "page",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#Bayesian-Linear-Mixed-Models-(Genomic-Data)-1",
    "page": "Linear Mixed Model (Genomic data)",
    "title": "Bayesian Linear Mixed Models (Genomic Data)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#Univariate-Linear-Mixed-Model-(Genomic-data)-1",
    "page": "Linear Mixed Model (Genomic data)",
    "title": "Univariate Linear Mixed Model (Genomic data)",
    "category": "section",
    "text": ""
},

{
    "location": "examples/genomicBLMM.html#Multivariate-Linear-Mixed-Model-(Genomic-data)-1",
    "page": "Linear Mixed Model (Genomic data)",
    "title": "Multivariate Linear Mixed Model (Genomic data)",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public.html#Some-Theory-in-JWAS-1",
    "page": "Public",
    "title": "Some Theory in JWAS",
    "category": "section",
    "text": ""
},

{
    "location": "lib/internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals.html#internals-1",
    "page": "Internals",
    "title": "internals",
    "category": "section",
    "text": ""
},

]}
