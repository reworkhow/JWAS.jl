var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: JWAS)JWAS is a well-documented software platform based on Julia and an interactive Jupyter notebook for analyses of general univariate and multivariate Bayesian mixed effects models.  These models are especially useful for, but not limited to, routine single-trait and multi-trait genomic prediction and genome-wide association studies using either complete or incomplete genomic data (\"single-step\" methods). Currently, JWAS provides broad scope of analyses, e.g., a wide collection of Bayesian methods for whole-genome analyses, including shrinkage estimation and variable selection methods. The features of JWAS include:Univariate (single-trait) analysis\nMultivariate (multi-trait) analysis  \nNo limitations on fixed effects (e.g. herd-year, age, sex)                                                                    \nRandom effects other than markers (e.g. litter, pen)                                  \nRandom effects using pedigree information\nAdditive genetic effects\nMaternal effects\nRandom permanent environmental effects  \nCorrelated residuals		\nCorrelated random effects\nUnknown (or known) variance components\nCorrelated marker effects                                                                \nUse of genomic information                                                                                \nComplete genomic data                                      		\nIncomplete genomic data (singe-step)"
},

{
    "location": "index.html#Supporting-and-Citing-1",
    "page": "Home",
    "title": "Supporting and Citing",
    "category": "section",
    "text": "We hope the friendly user interface and fast computing speed of JWAS will provide power and convenience for users in both industry and academia to analyze large datasets. Further, as a well-documented open-source software tool, we hope JWAS will also be used by a group of active community members, who will contribute to the source code and help maintain the project. Junior scientists can understand and learn the methodologies for whole-genome analyses by using JWAS and reading the tutorials and source code.If you would like to help support JWAS, please star the repository on the upper right corner here as such statistic will help to demonstrate the active involvement of the community. If you use JWAS for your research, teaching, or other activities, we would be grateful if you could cite our work following this citation guideline."
},

{
    "location": "index.html#The-trouble,-the-error-and-the-new-feature-1",
    "page": "Home",
    "title": "The trouble, the error and the new feature",
    "category": "section",
    "text": "If you have trouble using JWAS, want new features or find errors in JWAS, please open an issue or contact <qtlcheng@ucdavis.edu>."
},

{
    "location": "index.html#Tutorials-1",
    "page": "Home",
    "title": "Tutorials",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Theory-1",
    "page": "Home",
    "title": "Theory",
    "category": "section",
    "text": "Pages = [\n  \"theory/theory.md\"\n]\nDepth = 4"
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\n  \"manual/getstarted.md\",\n  \"manual/workflow.md\",\n  \"manual/public.md\",\n  \"manual/internals.md\",\n]\nDepth = 3"
},

{
    "location": "index.html#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "Pages = [\n  \"examples/conventionalBLMM.md\"\n  \"examples/genomicBLMM.md\"\n  \"examples/LinearAdditiveGeneticModel.md\"\n]\nDepth = 4"
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
    "location": "theory/theory.html#multiple-trait-Bayesian-Alphabet-1",
    "page": "Some Theory",
    "title": "multiple-trait Bayesian Alphabet",
    "category": "section",
    "text": "In multi-trait BayesCPi, the prior for alpha_lk, the marker effect of trait k for locus l, is a mixture with a point mass at zero and a univariate normal distribution conditional on sigma_k^2:beginalign*\nalpha_lkmidpi_ksigma_k^2  begincases\nsim Nleft(0sigma_k^2right)  probability(1-pi_k)\n0  probabilitypi_k\nendcases\nendalign*and the covariance between effects for traits k and k at the same locus, i.e., alpha_lk and alpha_lk^ isbeginalign*\ncovleft(alpha_lkalpha_lk^midsigma_kk^right)=begincases\nsigma_kk^  ifbothalpha_lkneq0andalpha_lk^neq0\n0  otherwise\nendcases\nendalign*The vector of marker effects at a particular locus boldsymbolalpha_l is written as boldsymbolalpha_l=boldsymbolD_lboldsymbolbeta_l, where boldsymbolD_l is a diagonal matrix with elements diagleft(boldsymbolD_lright)=boldsymboldelta_l=left(delta_l1delta_l2delta_l3ldotsdelta_ltright), where delta_lk is an indicator variable indicating whether the marker effect of locus l for trait k is zero or non-zero, and the vector boldsymbolbeta_l follows a multivariate normal distribution with null mean and covariance matrix boldsymbolG. The covariance matrix boldsymbolG is a priori assumed to follow an inverse Wishart distribution, W_t^-1left(mathbfS_betanu_betaright).In the most general case, any marker effect might be zero for any possible combination of t traits resulting in 2^t possible combinations of boldsymboldelta_l. For example, in a t=2 trait model, there are 2^t=4 combinations for  boldsymboldelta_l: (00), (01), (10), (11). Suppose in general we use numerical labels \"1\", \"2\",ldots, \"l\" for the 2^t possible outcomes for  boldsymboldelta_l, then the prior for  boldsymboldelta_l is a categorical distributionbeginalign*\n  pleft(boldsymboldelta_l=iright)\n=  Pi_1Ileft(boldsymboldelta_l=1right)+Pi_2Ileft(boldsymboldelta_l=2right)++Pi_lIleft(boldsymboldelta_l=lright)\nendalign*where sum_i=1^lPi_i=1 with Pi_i being the prior probability that the vector boldsymboldelta_l corresponds to the vector labelled i. A Dirichlet distribution with all parameters equal to one, i.e., a uniform distribution, can be used for the prior for boldsymbolPi=left(Pi_1Pi_2Pi_lright).   The differences in multi-trait BayesB method is that the prior for boldsymbolbeta_l is a multivariate t distribution, rather than a multivariate normal distribution. This is equivalent to assuming boldsymbolbeta_l has a multivariate normal distribution with null mean and locus-specific covariance matrix boldsymbolG_l, which is assigned an inverse Wishart prior, W_t^-1left(mathbfS_betanu_betaright). Multi-trait BayesA method is a special case of multi-trait BayesB method where boldsymboldelta_l is always a vector of ones.referencesMeuwissen T, Hayes B, Goddard M. Prediction of total genetic value using genome-wide dense marker maps. Genetics, 157:1819–1829, 2001.\nFernando R, Garrick D. Bayesian methods applied to GWAS. Methods Mol Biol. 2013;1019:237–274.\nCheng H, Garrick D, Fernando R. A fast and efficient Gibbs sampler for BayesB in whole- genome analyses. Genet Sel Evol, 2015, 47:80.\nFernando R, Dekkers J,Garrick D. A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses. Genetics Selection Evolution, 2015 46(1), 50.\nCheng H, Kizilkaya K, Zeng J, Garrick D, Fernando R. Genomic Prediction from Multiple-trait Bayesian Regression Methods using Mixture Priors. Genetics. 2018"
},

{
    "location": "manual/getstarted.html#",
    "page": "Get Started",
    "title": "Get Started",
    "category": "page",
    "text": ""
},

{
    "location": "manual/getstarted.html#Get-Started-1",
    "page": "Get Started",
    "title": "Get Started",
    "category": "section",
    "text": ""
},

{
    "location": "manual/getstarted.html#Installation-1",
    "page": "Get Started",
    "title": "Installation",
    "category": "section",
    "text": "To install julia, please go to the offical Julia website. Please see platform specific instructions if you have trouble installing Julia.To install the package, use the following command inside the Julia REPL (or IJulia Notebook):Pkg.add(\"JWAS\")To load the JWAS package, use the following command inside the Julia REPL (or IJulia Notebook):using JWASThe command Pkg.add(\"JWAS\") will add the registered official JWAS.jl and dependencies.To use the latest/beta features under development, run Pkg.checkout(\"JWAS\") to get the newest unofficial JWAS. Run Pkg.free(\"JWAS\") to go back to the offical one."
},

{
    "location": "manual/getstarted.html#IJulia-notebook-1",
    "page": "Get Started",
    "title": "IJulia notebook",
    "category": "section",
    "text": "If you prefer “reproducible research”, an interactive Jupyter notebook interface is available for Julia (and therefore JWAS). The Jupyter notebook is an open-source web application for creating and sharing documents that contain live code, equations, visualizations and explanatory text. To install IJulia, please go to IJulia."
},

{
    "location": "manual/getstarted.html#Standalone-application-1",
    "page": "Get Started",
    "title": "Standalone application",
    "category": "section",
    "text": "A fully self-contained application for JWAS (no installation required) will come out this year."
},

{
    "location": "manual/getstarted.html#Access-documentation-1",
    "page": "Get Started",
    "title": "Access documentation",
    "category": "section",
    "text": "To show the basic information (README file) of JWAS in REPL or IJulia notebook using ?JWAS and press enter (Please load the JWAS package at first).For help on a specific function, type ? followed by its name, e.g. ?runMCMC and press enter in REPL or IJulia notebook (Please load the JWAS package at first).The full documentation is available here."
},

{
    "location": "manual/getstarted.html#run-your-analysis-1",
    "page": "Get Started",
    "title": "run your analysis",
    "category": "section",
    "text": "There are several ways to run you analysis."
},

{
    "location": "manual/workflow.html#",
    "page": "Workflow",
    "title": "Workflow",
    "category": "page",
    "text": ""
},

{
    "location": "manual/workflow.html#Workflow-1",
    "page": "Workflow",
    "title": "Workflow",
    "category": "section",
    "text": ""
},

{
    "location": "manual/workflow.html#Data-format-1",
    "page": "Workflow",
    "title": "Data format",
    "category": "section",
    "text": "#data.txt\nid,y1,y2,y3,x1,x2,x3,x4,dam\n\n#pedigree\nid,sire,dam\n\n#genotype\nid,m1,m2,m3,m4,m5\n"
},

{
    "location": "manual/workflow.html#Read-data-1",
    "page": "Workflow",
    "title": "Read data",
    "category": "section",
    "text": "data = CSV.read(\"data.txt\")"
},

{
    "location": "manual/workflow.html#Build-Model-Equations-1",
    "page": "Workflow",
    "title": "Build Model Equations",
    "category": "section",
    "text": "model_equation = \"y1 = x1 + x2 + x3 + x4;\n                  y2 = x1 + x2 + x3*x4\"\nmodel=build_model(model_equation)link to build_model"
},

{
    "location": "manual/workflow.html#Set-Factors-or-Covariate-1",
    "page": "Workflow",
    "title": "Set Factors or Covariate",
    "category": "section",
    "text": "set_covariate(\"x1\")link to set_covariate"
},

{
    "location": "manual/workflow.html#Set-Random-or-Fixed-Effects-1",
    "page": "Workflow",
    "title": "Set Random or Fixed Effects",
    "category": "section",
    "text": "set_random(\"x2\",)link to set_random"
},

{
    "location": "manual/workflow.html#Use-Pedigree-Information-1",
    "page": "Workflow",
    "title": "Use Pedigree Information",
    "category": "section",
    "text": "ped=get_pedigree(\"pedigree.txt\")link to get_pedigree"
},

{
    "location": "manual/workflow.html#Use-Genomic-Information-1",
    "page": "Workflow",
    "title": "Use Genomic Information",
    "category": "section",
    "text": "add_genotypes(model,\"genotypes.txt\")link to add_genotypeslink to Workflow"
},

{
    "location": "manual/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "manual/public.html#Public-functions-1",
    "page": "Public",
    "title": "Public functions",
    "category": "section",
    "text": "Documentation for JWAS.jl\'s public interface, which are available to general users."
},

{
    "location": "manual/public.html#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]\nModules = [JWAS]Pages = [\"public.md\"]\nModules = [JWAS.misc]"
},

{
    "location": "manual/public.html#JWAS.build_model",
    "page": "Public",
    "title": "JWAS.build_model",
    "category": "function",
    "text": "build_model(model_equations::AbstractString,R;df::Float64=4.0)\n\nBuild models from model equations with residual varainces R and degree of freedom for residual variance df defaulting to 4.0.\nBy default, all variabels in model_equations are fixed and factors. Set variables to be covariates or random using functions set_covariate() or set_random().\n\n#single-trait\nmodel_equations = \"BW = intercept + age + sex\"\nR               = 6.72\nmodels          = build_model(model_equations,R);\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex;\n                   CW = intercept + litter\";\nR               = [6.72   24.84\n                   24.84  708.41]\nmodels          = build_model(model_equations,R);\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.set_covariate",
    "page": "Public",
    "title": "JWAS.set_covariate",
    "category": "function",
    "text": "set_covariate(mme::MME,variables::AbstractString...)\n\nset variables as covariates; mme is the output of function build_model().\n\n#After running build_model, variabels age and year can be set to be covariates as\nset_covariate(models,\"age\",\"year\")\n#or\nset_covariate(models,\"age year\")\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.set_random",
    "page": "Public",
    "title": "JWAS.set_random",
    "category": "function",
    "text": "set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)\n\nset variables as random polygenic effects with pedigree information ped, variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + Age + Animal\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = 1.6\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#single-trait (example 2)\nmodel_equation  = \"y = intercept + Age + Animal + Animal*Age\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = [1.6   0.2\n                   0.2  1.0]\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex + Animal\n                   CW = intercept + age + sex + Animal\"\nmodel           = build_model(model_equations,R);\nped             = get_pedigree(pedfile);\nG               = [6.72   2.84\n                   2.84  8.41]\nset_random(model,\"Animal\", ped,G)\n\n\n\nset_random(mme::MME,randomStr::AbstractString,G;df=4)\n\nset variables as i.i.d random effects with variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + litter + sex\"\nmodel           = build_model(model_equation,R)\nG               = 0.6\nset_random(model,\"litter\",G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + litter + sex\n                   CW = intercept + litter + sex\"\nmodel           = build_model(model_equations,R);\nG               = [3.72  1.84\n                   1.84  3.41]\nset_random(model,\"litter\",G)\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.get_pedigree",
    "page": "Public",
    "title": "JWAS.get_pedigree",
    "category": "function",
    "text": "get_pedigree(pedfile::AbstractString)\n\nGet pedigree informtion from a pedigree file.\nFile format:\n\na 0 0\nb 0 0\nc a b\nd a c\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.add_genotypes",
    "page": "Public",
    "title": "JWAS.add_genotypes",
    "category": "function",
    "text": "add_genotypes(mme::MME,file,G;separator=\' \',header=true,center=true,G_is_marker_variance=false,df=4.0)\n\nGet marker informtion from a genotype file (same order as the phenotype file).\nG defaults to the genetic variance with degree of freedom df=4.0.\nFile format:\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,2,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.add_markers",
    "page": "Public",
    "title": "JWAS.add_markers",
    "category": "function",
    "text": "same to add_genotypes\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.outputMCMCsamples",
    "page": "Public",
    "title": "JWAS.outputMCMCsamples",
    "category": "function",
    "text": "outputMCMCsamples(mme::MME,trmStr::AbstractString...)\n\nGet samples for specific variables.\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.showMME",
    "page": "Public",
    "title": "JWAS.showMME",
    "category": "function",
    "text": "showMME(mme::MME,df::DataFrame)\n\nShow left-hand side and right-hand side of mixed model equations (no markers).\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.solve",
    "page": "Public",
    "title": "JWAS.solve",
    "category": "function",
    "text": "solve(mme::MME,df::DataFrame;solver=\"default\",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)\n\nSolve the mixed model equations (no marker information) without estimating variance components.\n\nAvailable solvers includes default,Jacobi,GaussSeidel,Gibbs sampler.\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.runMCMC",
    "page": "Public",
    "title": "JWAS.runMCMC",
    "category": "function",
    "text": "runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,starting_value=false,printout_frequency=100,missing_phenotypes=false,constraint=false,methods=\"conventional (no markers)\",output_samples_frequency::Int64 = 0)\n\nRun MCMC (marker information included or not) with sampling of variance components.\n\navailable methods include \"conventional (no markers)\", \"BayesC0\", \"BayesC\", \"BayesCC\",\"BayesB\".\nmissing_phenotypes\nPi for single-trait analyses is a number; Pi for multi-trait analyses is a dictionary such as Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1),\nif Pi (Π) is not provided in multi-trait analysis, it will be generated assuming all markers have effects on all traits.\nsave MCMC samples every output_samples_frequency iterations\nstarting_value can be provided as a vector for all location parameteres except marker effects.\nprint out the monte carlo mean in REPL with printout_frequency\nconstraint=true if constrain residual covariances between traits to be zero.\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.misc.get_additive_genetic_variances",
    "page": "Public",
    "title": "JWAS.misc.get_additive_genetic_variances",
    "category": "function",
    "text": "get_additive_genetic_variances(model::MME,files...;header=true)\n\nGet MCMC samples for additive genetic variances using samples of marker effects stored in files.\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.misc.get_heritability",
    "page": "Public",
    "title": "JWAS.misc.get_heritability",
    "category": "function",
    "text": "get_heritability(samples_for_genetic_variances::Array{Array{Float64,2},1},samples_for_residual_vairances::Array{Array{Float64,2},1}))\n\nGet MCMC samples for heritabilities.\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.misc.get_correlations",
    "page": "Public",
    "title": "JWAS.misc.get_correlations",
    "category": "function",
    "text": "get_correlations(samples_for_genetic_variances::Array{Array{Float64,2},1})\n\nGet MCMC samples for correlations\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.misc.report",
    "page": "Public",
    "title": "JWAS.misc.report",
    "category": "function",
    "text": "report(X::Array{Array{Float64,2},1};index=false)\n\nshow summary statistics for MCMC samples (matrices or index [i,j] of matrices)\n\n\n\nreport(X::Array{Array{Float64,1},1};index=false)\n\nshow summary statistics for MCMC samples (vectors or index i of vectors)\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.misc.QC",
    "page": "Public",
    "title": "JWAS.misc.QC",
    "category": "function",
    "text": "QC(infile,outfile;separator=\' \',header=true,missing=false,MAF=0.1)\n\nQuality control for input file infile, then write out to output file: outfile.\nDelete loci with minor allele frequency < MAF.\nmissing genotypes are replaced by column means.\nFile format (header=true,separator=\',\',missing=9):\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,9,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n"
},

{
    "location": "manual/public.html#JWAS.misc.get_breeding_values",
    "page": "Public",
    "title": "JWAS.misc.get_breeding_values",
    "category": "function",
    "text": "get_breeding_values(model::MME,files...;header=true)\n\nGet esitimated breeding values and prediction error variances using samples of marker effects stored in files.\n\n\n\n"
},

{
    "location": "manual/public.html#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": "build_model\nset_covariate\nset_random\nget_pedigree\nadd_genotypes\nadd_markers\noutputMCMCsamples\nshowMME\nsolve\nrunMCMCJWAS.misc.get_additive_genetic_variances\nJWAS.misc.get_heritability\nJWAS.misc.get_correlations\nJWAS.misc.report\nJWAS.misc.QC\nJWAS.misc.get_breeding_values"
},

{
    "location": "manual/internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "manual/internals.html#Internal-functions-1",
    "page": "Internals",
    "title": "Internal functions",
    "category": "section",
    "text": "Documentation for JWAS.jl\'s internal (private) interface, which are not available to general users. These internal functions are small blocks that public function build on."
},

{
    "location": "manual/internals.html#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internals.md\"]"
},

{
    "location": "manual/internals.html#JWAS.add_genotypes-Tuple{JWAS.MME,Any,Any}",
    "page": "Internals",
    "title": "JWAS.add_genotypes",
    "category": "method",
    "text": "add_genotypes(mme::MME,file,G;separator=\' \',header=true,center=true,G_is_marker_variance=false,df=4.0)\n\nGet marker informtion from a genotype file (same order as the phenotype file).\nG defaults to the genetic variance with degree of freedom df=4.0.\nFile format:\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,2,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.add_markers-Tuple{JWAS.MME,Any,Any}",
    "page": "Internals",
    "title": "JWAS.add_markers",
    "category": "method",
    "text": "same to add_genotypes\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.build_model-Tuple{AbstractString,Any}",
    "page": "Internals",
    "title": "JWAS.build_model",
    "category": "method",
    "text": "build_model(model_equations::AbstractString,R;df::Float64=4.0)\n\nBuild models from model equations with residual varainces R and degree of freedom for residual variance df defaulting to 4.0.\nBy default, all variabels in model_equations are fixed and factors. Set variables to be covariates or random using functions set_covariate() or set_random().\n\n#single-trait\nmodel_equations = \"BW = intercept + age + sex\"\nR               = 6.72\nmodels          = build_model(model_equations,R);\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex;\n                   CW = intercept + litter\";\nR               = [6.72   24.84\n                   24.84  708.41]\nmodels          = build_model(model_equations,R);\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.get_pedigree-Tuple{AbstractString}",
    "page": "Internals",
    "title": "JWAS.get_pedigree",
    "category": "method",
    "text": "get_pedigree(pedfile::AbstractString)\n\nGet pedigree informtion from a pedigree file.\nFile format:\n\na 0 0\nb 0 0\nc a b\nd a c\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.outputMCMCsamples-Tuple{JWAS.MME,Vararg{AbstractString,N} where N}",
    "page": "Internals",
    "title": "JWAS.outputMCMCsamples",
    "category": "method",
    "text": "outputMCMCsamples(mme::MME,trmStr::AbstractString...)\n\nGet samples for specific variables.\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.runMCMC-Tuple{Any,Any}",
    "page": "Internals",
    "title": "JWAS.runMCMC",
    "category": "method",
    "text": "runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,starting_value=false,printout_frequency=100,missing_phenotypes=false,constraint=false,methods=\"conventional (no markers)\",output_samples_frequency::Int64 = 0)\n\nRun MCMC (marker information included or not) with sampling of variance components.\n\navailable methods include \"conventional (no markers)\", \"BayesC0\", \"BayesC\", \"BayesCC\",\"BayesB\".\nmissing_phenotypes\nPi for single-trait analyses is a number; Pi for multi-trait analyses is a dictionary such as Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1),\nif Pi (Π) is not provided in multi-trait analysis, it will be generated assuming all markers have effects on all traits.\nsave MCMC samples every output_samples_frequency iterations\nstarting_value can be provided as a vector for all location parameteres except marker effects.\nprint out the monte carlo mean in REPL with printout_frequency\nconstraint=true if constrain residual covariances between traits to be zero.\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.set_covariate-Tuple{JWAS.MME,Vararg{AbstractString,N} where N}",
    "page": "Internals",
    "title": "JWAS.set_covariate",
    "category": "method",
    "text": "set_covariate(mme::MME,variables::AbstractString...)\n\nset variables as covariates; mme is the output of function build_model().\n\n#After running build_model, variabels age and year can be set to be covariates as\nset_covariate(models,\"age\",\"year\")\n#or\nset_covariate(models,\"age year\")\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.set_random-Tuple{JWAS.MME,AbstractString,Any}",
    "page": "Internals",
    "title": "JWAS.set_random",
    "category": "method",
    "text": "set_random(mme::MME,randomStr::AbstractString,G;df=4)\n\nset variables as i.i.d random effects with variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + litter + sex\"\nmodel           = build_model(model_equation,R)\nG               = 0.6\nset_random(model,\"litter\",G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + litter + sex\n                   CW = intercept + litter + sex\"\nmodel           = build_model(model_equations,R);\nG               = [3.72  1.84\n                   1.84  3.41]\nset_random(model,\"litter\",G)\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.set_random-Tuple{JWAS.MME,AbstractString,JWAS.PedModule.Pedigree,Any}",
    "page": "Internals",
    "title": "JWAS.set_random",
    "category": "method",
    "text": "set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)\n\nset variables as random polygenic effects with pedigree information ped, variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + Age + Animal\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = 1.6\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#single-trait (example 2)\nmodel_equation  = \"y = intercept + Age + Animal + Animal*Age\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = [1.6   0.2\n                   0.2  1.0]\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex + Animal\n                   CW = intercept + age + sex + Animal\"\nmodel           = build_model(model_equations,R);\nped             = get_pedigree(pedfile);\nG               = [6.72   2.84\n                   2.84  8.41]\nset_random(model,\"Animal\", ped,G)\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.showMME-Tuple{JWAS.MME,DataFrames.DataFrame}",
    "page": "Internals",
    "title": "JWAS.showMME",
    "category": "method",
    "text": "showMME(mme::MME,df::DataFrame)\n\nShow left-hand side and right-hand side of mixed model equations (no markers).\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.solve-Tuple{JWAS.MME,DataFrames.DataFrame}",
    "page": "Internals",
    "title": "JWAS.solve",
    "category": "method",
    "text": "solve(mme::MME,df::DataFrame;solver=\"default\",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)\n\nSolve the mixed model equations (no marker information) without estimating variance components.\n\nAvailable solvers includes default,Jacobi,GaussSeidel,Gibbs sampler.\n\n\n\n"
},

{
    "location": "manual/internals.html#JWAS.getMME-Tuple{JWAS.MME,DataFrames.DataFrame}",
    "page": "Internals",
    "title": "JWAS.getMME",
    "category": "method",
    "text": "Construct mixed model equations with\n\nincidence matrix: X      ; response        : ySparse; left-hand side  : mmeLhs ; right-hand side : mmeLhs ;\n\n\n\n"
},

{
    "location": "manual/internals.html#Internal-interface-1",
    "page": "Internals",
    "title": "Internal interface",
    "category": "section",
    "text": "Modules = [JWAS,JWAS.PedModule]\nOrder   = [:function, :type]"
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
    "text": "using DataFrames,CSV,JWAS,JWAS.Datasets;phenofile = Datasets.dataset(\"testMME\",\"data.txt\");\ndata      = CSV.read(phenofile,delim = \',\',header=true);\nhead(data)\n\n# 10×9 DataFrames.DataFrame\n# │ Row │ sow       │ site │ yr   │ age │ geneticCode │ parity │ nwn │ SYS           │ bw   │\n# ├─────┼───────────┼──────┼──────┼─────┼─────────────┼────────┼─────┼───────────────┼──────┤\n# │ 1   │ 100-113   │ 113  │ 2005 │ 18  │ P1          │ 1      │ 8   │ 113_2005_WNTR │ 9.0  │\n# │ 2   │ 100-113   │ 113  │ 2006 │ 18  │ P1          │ 2      │ 12  │ 113_2006_SPNG │ 8.0  │\n# │ 3   │ 100-5     │ 5    │ 2008 │ 15  │ P2          │ 1      │ 10  │ 5_2008_ATMN   │ 7.5  │\n# │ 4   │ 1000-5    │ 5    │ 2009 │ 17  │ P2          │ 1      │ 10  │ 5_2009_SPNG   │ 8.3  │\n# │ 5   │ 10000-131 │ 13   │ 2004 │ 16  │ Commercial  │ 1      │ 9   │ 13_2004_WNTR  │ 4.3  │\n# │ 6   │ 10000-131 │ 13   │ 2004 │ 18  │ Commercial  │ 2      │ 10  │ 13_2004_SMMR  │ 2.8  │model_equation    = \"nwn = intercept +parity + parity*site + yr + geneticCode + age\"\n\nresidual_variance = 2.97;\nmodel             = build_model(model_equation,residual_variance)\n\nset_covariate(model,\"age\");\n\ngeneticCode_variance = 0.26;\nset_random(model,\"geneticCode\",geneticCode_variance);outputMCMCsamples(model,\"parity\",\"age\");out=runMCMC(data,model,chain_length=50000,output_samples_frequency=100);a = 1\nb = 2\na + b"
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
    "text": "using DataFrames,CSV,JWAS,JWAS.Datasets;phenofile = Datasets.dataset(\"testMME\",\"data.txt\");\ndata      = CSV.read(phenofile,delim = \',\',header=true);\nhead(data);model_equation    = \"nwn = intercept +parity + parity*site + yr + geneticCode + age\"\n\nresidual_variance = 2.97;\nmodel             = build_model(model_equation,residual_variance)\n\nset_covariate(model,\"age\");\n\ngeneticCode_variance = 0.26;\nset_random(model,\"geneticCode\",geneticCode_variance);outputMCMCsamples(model,\"parity\",\"age\");out=runMCMC(data,model,chain_length=50000,output_samples_frequency=100)"
},

]}
