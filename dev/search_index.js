var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: JWAS)JWAS is a well-documented software platform based on Julia and an interactive Jupyter notebook for analyses of general univariate and multivariate Bayesian mixed effects models.  These models are especially useful for, but not limited to, routine single-trait and multi-trait genomic prediction and genome-wide association studies using either complete or incomplete genomic data (\"single-step\" methods). Currently, JWAS provides broad scope of analyses, e.g., a wide collection of Bayesian methods for whole-genome analyses, including shrinkage estimation and variable selection methods. The features of JWAS include:Univariate (single-trait) analysis\nMultivariate (multi-trait) analysis  \nNo limitations on fixed effects (e.g., herd, year, age, sex)\nRandom effects other than markers (e.g., litter, pen)                                  \nRandom effects using pedigree information\nAdditive genetic effects\nMaternal effects\nRandom permanent environmental effects  \nCorrelated residuals		\nCorrelated random effects\nUnknown (or known) variance components\nUse of genomic information\nComplete genomic data                                      		\nIncomplete genomic data (singe-step)"
},

{
    "location": "#Supporting-and-Citing-1",
    "page": "Home",
    "title": "Supporting and Citing",
    "category": "section",
    "text": "We hope the friendly user interface and fast computing speed of JWAS will provide power and convenience for users in both industry and academia to analyze large datasets. Further, as a well-documented open-source software tool, we hope JWAS will also be used by a group of active community members, who will contribute to the source code and help maintain the project. Junior scientists can understand and learn the methodologies for whole-genome analyses by using JWAS and reading the tutorials and source code.If you would like to help support JWAS, please star the repository on the upper right corner here as such statistic will help to demonstrate the active involvement of the community. If you use JWAS for your research, teaching, or other activities, we would be grateful if you could cite our work following Citing."
},

{
    "location": "#The-trouble,-the-error-and-the-new-feature-1",
    "page": "Home",
    "title": "The trouble, the error and the new feature",
    "category": "section",
    "text": "If you have trouble using JWAS, want new features or find errors in JWAS, please open an issue or contact <qtlcheng@ucdavis.edu>."
},

{
    "location": "#Tutorials-1",
    "page": "Home",
    "title": "Tutorials",
    "category": "section",
    "text": ""
},

{
    "location": "#Theory-1",
    "page": "Home",
    "title": "Theory",
    "category": "section",
    "text": "Pages = [\n  \"theory/theory.md\"\n]\nDepth = 3"
},

{
    "location": "#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\n  \"manual/getstarted.md\",\n  \"manual/workflow.md\",\n  \"manual/public.md\",\n  \"manual/internals.md\",\n]\nDepth = 3"
},

{
    "location": "#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "Pages = [\n  \"examples/examples.md\"\n]\nDepth = 2"
},

{
    "location": "theory/theory/#",
    "page": "Some Theory",
    "title": "Some Theory",
    "category": "page",
    "text": ""
},

{
    "location": "theory/theory/#Some-Theory-in-JWAS-1",
    "page": "Some Theory",
    "title": "Some Theory in JWAS",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory/#A-Table-for-Bayesian-Linear-Mixed-Models-(BLMM)-1",
    "page": "Some Theory",
    "title": "A Table for Bayesian Linear Mixed Models (BLMM)",
    "category": "section",
    "text": "(Image: BLMM)"
},

{
    "location": "theory/theory/#Models-1",
    "page": "Some Theory",
    "title": "Models",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory/#Complete-Genomic-Data-1",
    "page": "Some Theory",
    "title": "Complete Genomic Data",
    "category": "section",
    "text": "The general form of the multivariate (univariate) mixed effects model for individual i from n individuals with complete genomic data in JWAS is\nmathbfy_i\n =sum_j=1^p_betaX_ijboldsymbolbeta_j+sum_k=1^p_uZ_ikmathbfu_k\n +sum_l=1^pM_ilboldsymbolalpha_l+mathbfe_i(1)where mathbfy_i is a vector of phenotypes of t traits for individual i; X_ij is the incidence matrix covariate corresponding to the jth fixed effect for individual i; boldsymbolbeta_j is a vector of jth fixed effects for the t traits; Z_ik is the incidence matrix covariate corresponding to the kth random effect for individual i; boldsymbolu_k is a vector of the kth random effects of t traits; M_il is the genotype covariate at locus l for individual i, p is the number of genotyped loci (each coded as 0,1,2), boldsymbolalpha_l is a vector of allele substitution effects or marker effects of t traits for locus j, and mathbfe_i is the vector of random residual effects of t traits for individual i. The JWAS implementation of this model involves missing phenotypes being imputed at each iteration of MCMC \\cite{sorensenGianolaBook} so that all individuals have observations for all traits. Note that when the number of traits t=1, the general form above simplifies to the single-trait  mixed effects model, and all vectors of effects in equation (1) become scalars."
},

{
    "location": "theory/theory/#Incomplete-Genomic-Data-1",
    "page": "Some Theory",
    "title": "Incomplete Genomic Data",
    "category": "section",
    "text": "The general form of the multivariate (univariate) mixed effects model with incomplete genomic data (\"single-step\" methods) for non-genotyped individuals ismathbfy_i\n=sum_j=1^p_betaX_ijboldsymbolbeta_j+sum_k=1^p_uZ_ikmathbfu_k+\nsum_l=1^phatM_ilboldsymbolalpha_l+sum_m=1^p_epsilonZ_nimboldsymbolepsilon_m+boldsymbole_i (2)where mathbfy_i is a vector of phenotypes of t traits for non-genotyped individual i;  hatM_il is the imputed genotype covariate at locus l for non-genotyped individual i, Z_nim is the incidence matrix covariate corresponding to the mth imputation residual for individual i and boldsymbolepsilon_i is a vector of imputation residuals. W_im is the incidence matrix covariate corresponding to the mth random effect for individual i. That vector of imputation residuals, boldsymbolepsilon=beginbmatrixboldsymbolepsilon_1^T  boldsymbolepsilon_2^T  ldots  endbmatrix^T, are a priori assumed to be Nleft(0(mathbfA_nn-mathbfA_ngmathbfA_gg^-1mathbfA_gn)otimesmathbfG_gright), where mathbfA_nn is the partition of the numerator relationship matrix  mathbfA that corresponds to non-genotyped individuals, mathbfA_ng or its transpose mathbfA_gn are partitions of mathbfA corresponding to relationships between non-genotyped and genotyped individuals or vice versa,  mathbfA_gg is the  partition of mathbfA that corresponds to genotyped animals, and mathbfG_g is the additive genetic covariance matrix. All the other variables are the same as in equation (1)."
},

{
    "location": "theory/theory/#Priors-1",
    "page": "Some Theory",
    "title": "Priors",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory/#Priors-for-effects-other-than-markers-1",
    "page": "Some Theory",
    "title": "Priors for effects other than markers",
    "category": "section",
    "text": "The fixed effects are assigned flat priors. The vector of random effects, mathbfu=beginbmatrixmathbfu_1^T  mathbfu_2^T  ldots  mathbfu_p_2^Tendbmatrix^T, are a priori assumed to be Nleft(0mathbfAotimesmathbfGright) with various options for mathbfA. For example, mathbfA could be an identity matrix if boldsymbolu_k is assumed to be independently and identically distributed. mathbfA can be  the numerator relationship matrix, when boldsymbolu is a vector of polygenic effects and mathbfG represents the additive-genetic variance not explained by molecular markers. Note that boldsymbolu can also be a concatenation of vectors of different types of random effects, such as litter, pen, polygenic and maternal effects. The vector boldsymbole_i of residuals are a priori assumed to be independently and identically following multivariate normal distributions with null mean and covariance matrix mathbfR, which in turn is a priori assumed to have an inverse Wishart prior distribution, W_t^-1left(mathbfS_enu_eright). Note that when number of traits t=1, the priors for mathbfG and mathbfR in single-trait analyses follow scaled inverted chi-square distributions."
},

{
    "location": "theory/theory/#Priors-for-marker-effects-1",
    "page": "Some Theory",
    "title": "Priors for marker effects",
    "category": "section",
    "text": ""
},

{
    "location": "theory/theory/#single-trait-BayesA-1",
    "page": "Some Theory",
    "title": "single-trait BayesA",
    "category": "section",
    "text": "The prior assumption is that marker effects have identical and independent univariate-t distributions each with a null mean, scale parameter S^2_alpha and nu degrees of freedom. This is equivalent to assuming that the marker effect at locus i has a univariate normal with null mean and unknown, locus-specific variance sigma^2_i, which in turn is assigned a scaled inverse chi-square prior with scale parameter S^2_alpha and nu_alpha degrees of freedom."
},

{
    "location": "theory/theory/#single-trait-BayesB-1",
    "page": "Some Theory",
    "title": "single-trait BayesB",
    "category": "section",
    "text": "In BayesB, the prior assumption is that marker effects have identical and independent mixture distributions, where each has a point mass at zero with probability pi and a univariate-t distribution with probability 1-pi having a null mean, scale parameter S^2_alpha and nu degrees of freedom. Thus, BayesA is a special case of BayesB with pi=0. Further, as in BayesA, the t-distribution in BayesB is equivalent to a univariate normal with null mean and unknown, locus-specific variance, which in turn is assigned a scaled inverse chi-square prior with scale parameter S^2_alpha and nu_alpha degrees of freedom. (A fast and efficient Gibbs sampler was implemented for BayesB in JWAS.)"
},

{
    "location": "theory/theory/#single-trait-BayesC-and-BayesC\\pi-1",
    "page": "Some Theory",
    "title": "single-trait BayesC and BayesCpi",
    "category": "section",
    "text": "In BayesC, the prior assumption is that marker effects have identical and independent mixture distributions, where each has a point mass at zero with probability pi and a univariate-normal distribution with probability 1-pi having a null mean and variance sigma^2_alpha, which in turn has a scaled inverse chi-square prior with scale parameter S^2_alpha and nu_alpha degrees of freedom. In addition to the above assumptions, in BayesC pi, pi is treated as unknown with a uniform prior."
},

{
    "location": "theory/theory/#multiple-trait-Bayesian-Alphabet-1",
    "page": "Some Theory",
    "title": "multiple-trait Bayesian Alphabet",
    "category": "section",
    "text": "In multi-trait BayesCPi, the prior for alpha_lk, the marker effect of trait k for locus l, is a mixture with a point mass at zero and a univariate normal distribution conditional on sigma_k^2:beginalign*\nalpha_lkmidpi_ksigma_k^2  begincases\nsim Nleft(0sigma_k^2right)  probability(1-pi_k)\n0  probabilitypi_k\nendcases\nendalign*and the covariance between effects for traits k and k at the same locus, i.e., alpha_lk and alpha_lk^ isbeginalign*\ncovleft(alpha_lkalpha_lk^midsigma_kk^right)=begincases\nsigma_kk^  ifbothalpha_lkneq0andalpha_lk^neq0\n0  otherwise\nendcases\nendalign*The vector of marker effects at a particular locus boldsymbolalpha_l is written as boldsymbolalpha_l=boldsymbolD_lboldsymbolbeta_l, where boldsymbolD_l is a diagonal matrix with elements diagleft(boldsymbolD_lright)=boldsymboldelta_l=left(delta_l1delta_l2delta_l3ldotsdelta_ltright), where delta_lk is an indicator variable indicating whether the marker effect of locus l for trait k is zero or non-zero, and the vector boldsymbolbeta_l follows a multivariate normal distribution with null mean and covariance matrix boldsymbolG. The covariance matrix boldsymbolG is a priori assumed to follow an inverse Wishart distribution, W_t^-1left(mathbfS_betanu_betaright).In the most general case, any marker effect might be zero for any possible combination of t traits resulting in 2^t possible combinations of boldsymboldelta_l. For example, in a t=2 trait model, there are 2^t=4 combinations for  boldsymboldelta_l: (00), (01), (10), (11). Suppose in general we use numerical labels \"1\", \"2\",ldots, \"l\" for the 2^t possible outcomes for  boldsymboldelta_l, then the prior for  boldsymboldelta_l is a categorical distributionbeginalign*\n  pleft(boldsymboldelta_l=iright)\n=  Pi_1Ileft(boldsymboldelta_l=1right)+Pi_2Ileft(boldsymboldelta_l=2right)++Pi_lIleft(boldsymboldelta_l=lright)\nendalign*where sum_i=1^lPi_i=1 with Pi_i being the prior probability that the vector boldsymboldelta_l corresponds to the vector labelled i. A Dirichlet distribution with all parameters equal to one, i.e., a uniform distribution, can be used for the prior for boldsymbolPi=left(Pi_1Pi_2Pi_lright).   The differences in multi-trait BayesB method is that the prior for boldsymbolbeta_l is a multivariate t distribution, rather than a multivariate normal distribution. This is equivalent to assuming boldsymbolbeta_l has a multivariate normal distribution with null mean and locus-specific covariance matrix boldsymbolG_l, which is assigned an inverse Wishart prior, W_t^-1left(mathbfS_betanu_betaright). Multi-trait BayesA method is a special case of multi-trait BayesB method where boldsymboldelta_l is always a vector of ones.referencesMeuwissen T, Hayes B, Goddard M. Prediction of total genetic value using genome-wide dense marker maps. Genetics, 157:1819–1829, 2001.\nFernando R, Garrick D. Bayesian methods applied to GWAS. Methods Mol Biol. 2013;1019:237–274.\nCheng H, Garrick D, Fernando R. A fast and efficient Gibbs sampler for BayesB in whole- genome analyses. Genet Sel Evol, 2015, 47:80.\nFernando R, Dekkers J,Garrick D. A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses. Genetics Selection Evolution, 2015 46(1), 50.\nCheng H, Kizilkaya K, Zeng J, Garrick D, Fernando R. Genomic Prediction from Multiple-trait Bayesian Regression Methods using Mixture Priors. Genetics. 2018"
},

{
    "location": "citing/citing/#",
    "page": "Citing",
    "title": "Citing",
    "category": "page",
    "text": ""
},

{
    "location": "citing/citing/#Citing-1",
    "page": "Citing",
    "title": "Citing",
    "category": "section",
    "text": ""
},

{
    "location": "citing/citing/#Software-1",
    "page": "Citing",
    "title": "Software",
    "category": "section",
    "text": "Cheng, H., Fernando, R. L., and Garrick, D. J. 2018 JWAS: Julia implementation of whole-genome analysis software. Proceedings of the World Congress on Genetics Applied to Livestock Production,11.859. Auckland, New Zealand."
},

{
    "location": "citing/citing/#Methods-1",
    "page": "Citing",
    "title": "Methods",
    "category": "section",
    "text": "Meuwissen T, Hayes B, Goddard M. 2001 Prediction of total genetic value using genome-wide dense marker maps. Genetics, 157:1819–1829.\nFernando R, Garrick D. 2013 Bayesian methods applied to GWAS. Methods Mol Biol., 1019:237–274.\nHabier D, Fernando R, Kizilkaya K, Garrick D. 2011 Extension of the bayesian alphabet for genomic selection. BMC Bioinformatics, 12(1), 186.\nCheng H, Garrick D, Fernando R. 2015 A fast and efficient Gibbs sampler for BayesB in whole-genome analyses. Genet Sel Evol, 47:80.\nCheng H, Kizilkaya K, Zeng J, Garrick D, Fernando R 2018 Genomic Prediction from Multiple-Trait Bayesian Regression Methods Using Mixture Priors. Genetics, 209(1): 89-103.\nFernando R, Dekkers J,Garrick D. 2015 A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses. Genetics Selection Evolution, 46(1), 50."
},

{
    "location": "manual/getstarted/#",
    "page": "Get Started",
    "title": "Get Started",
    "category": "page",
    "text": ""
},

{
    "location": "manual/getstarted/#Get-Started-1",
    "page": "Get Started",
    "title": "Get Started",
    "category": "section",
    "text": ""
},

{
    "location": "manual/getstarted/#Installation-1",
    "page": "Get Started",
    "title": "Installation",
    "category": "section",
    "text": "To install julia, please go to the offical Julia website. Please see platform specific instructions if you have trouble installing Julia.To install the package, use the following command inside the Julia REPL (or IJulia Notebook):Pkg.add(\"JWAS\")To load the JWAS package, use the following command inside the Julia REPL (or IJulia Notebook):using JWASThe command Pkg.add(\"JWAS\") will add the registered official JWAS package and dependencies.To use the latest/beta features under development, run Pkg.add(PackageSpec(name=\"JWAS\", rev=\"master\")) to get the newest unofficial JWAS. Run Pkg.free(\"JWAS\") to go back to the official one."
},

{
    "location": "manual/getstarted/#Jupyter-Notebook-1",
    "page": "Get Started",
    "title": "Jupyter Notebook",
    "category": "section",
    "text": "If you prefer “reproducible research”, an interactive Jupyter Notebook interface is available for Julia (and therefore JWAS). The Jupyter Notebook is an open-source web application for creating and sharing documents that contain live code, equations, visualizations and explanatory text. To install IJulia for Jupyter Notebook, please go to IJulia."
},

{
    "location": "manual/getstarted/#Docker-1",
    "page": "Get Started",
    "title": "Docker",
    "category": "section",
    "text": "note: Jupyter Notebooks with JWAS via Docker\nDocker provides a straightforward way to install Jupyter Notebooks with JWAS.Install Docker from here for your platform.\nFrom a terminal (on Mac or Linux), run the command:docker run -it --rm -p 8888:8888 qtlrocks/jwas-dockerThis will start a Jupyter Notebook server listening for HTTP connections on port 8888 with a randomly generated authentication token. Examples for JWAS can be accessed from the notebook: notebooks/0_index.ipynb.The directories and files created within the Docker container will be lost when the container is stopped. To save your work on the host machine, a directory on the host machine can be mounted as a folder in the container with the -v option. After cd into your working directory on your local machine or a server, run the commanddocker run -it --rm -p 8888:8888 -v `pwd`:/home/jovyan/work qtlrocks/jwas-dockerThis command creates a Docker container with the folder /home/jovyan/work with the contents of pwd of the host machine. Files and directories that are in the folder pwd will not be lost when the container is stopped.  After running this command, it is expected to prompt something like[I 10:41:54.774 NotebookApp] Writing notebook server cookie secret to /home/ubuntu/.local/share/jupyter/runtime/notebook_cookie_secret\n[I 10:41:54.920 NotebookApp] Serving notebooks from local directory: /home/ubuntu\n[I 10:41:54.920 NotebookApp] 0 active kernels\n[I 10:41:54.920 NotebookApp] The Jupyter Notebook is running at:\n[I 10:41:54.920 NotebookApp] http://0.0.0.0:8888/?token=75ad671f75b4c47be70591f46bec604997d8a9bd9dd51f0d\n[I 10:41:54.920 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).\n[C 10:41:54.921 NotebookApp]\n\n    Copy/paste this URL into your browser when you connect for the first time,\n    to login with a token:\n        http://0.0.0.0:8888/?token=75ad671f75b4c47be70591f46bec604997d8a9bd9dd51f0dThen, open the url in an internet browser (IE, Firefox, Chrome, Safari, etc) if JWAS-docker is launched on your local machine.If you prefer running scripts using linux commands in Bash instead of Jupyter Notebook, please run the commanddocker run -it --rm -v `pwd`:/home/jovyan/work qtlrocks/jwas-docker bash"
},

{
    "location": "manual/getstarted/#Standalone-application-1",
    "page": "Get Started",
    "title": "Standalone application",
    "category": "section",
    "text": "note: standalone application (no installation required)\nA fully self-contained application for JWAS (no installation required) will come out next year."
},

{
    "location": "manual/getstarted/#Access-documentation-1",
    "page": "Get Started",
    "title": "Access documentation",
    "category": "section",
    "text": "warning: Warning\nPlease load the JWAS package at first.To show the basic information (README file) of JWAS in REPL or IJulia notebook using ?JWAS and press enter.For help on a specific function, type ? followed by its name, e.g. ?runMCMC and press enter in REPL or IJulia notebook.The full documentation is available here."
},

{
    "location": "manual/getstarted/#Run-your-analysis-1",
    "page": "Get Started",
    "title": "Run your analysis",
    "category": "section",
    "text": "There are several ways to run your analysis.(1) The easiest way to run analysis in Julia is by starting an interactive session (REPL) by double-clicking the Julia executable or running julia from the command line (e.g., terminal) asjulia> 1+2\n3\n\njulia> 3*4\n12To evaluate code written in a file script.jl in REPL, write and runjulia> include(\"script.jl\").To exit the interactive session, type ^D – the control key together with the d key or type quit().(2) To run code in a file non-interactively from the command line (e.g.,termial), you can give it as the first argument to the julia command:julia script.jlIf you want to pass arguments to your script, run it asjulia script.jl arg1 arg2where arguments arg1 and arg2 are passed to your script as ARGS[1] and ARGS[2] of type String. Please see julia docs for more options.(3) To run code in Jupyter Notebook, please see IJulia.(4) To run code in Jupyter Notebook via Docker, please see Docker."
},

{
    "location": "manual/workflow/#",
    "page": "Workflow",
    "title": "Workflow",
    "category": "page",
    "text": ""
},

{
    "location": "manual/workflow/#Workflow-1",
    "page": "Workflow",
    "title": "Workflow",
    "category": "section",
    "text": "A step by step workflow for how to run JWAS is shown in this section. The workflow below is used to demonstrate a three-trait Bayesian linear mixed model fitting fixed effects (x1, x3), random effects (x2), direct genetic effects (ID), maternal genetic effects (dam) and genomic information.Pages = [\n  \"workflow.md\"\n]\nDepth = 2"
},

{
    "location": "manual/workflow/#Available-Models-1",
    "page": "Workflow",
    "title": "Available Models",
    "category": "section",
    "text": "Given the data and model equations, several different types of models are available in JWAS as shown below. In the table below, \"X\" denotes the type of available data, and \"Y<=A\" denotes that Y individuals is a subset of A individuals.  Linear Mixed Models (LMM) phenotypes (Y) pedigree (A) genotypes (G) notes\nConventional LMM X   \nPedigree-based LMM X X  Y<=A\nComplete Genomic LMM X maybe X Y<=G\nIncomplete Genomic LMM X X X Y<=A,G<=Anote: Note\nIncomplete Genomic LMM is also called \"single-step\" methods in animal breeding.\nPedigree information may be used in Complete Genomic LMM for extra polygenic effects to account for genetic variance not explained by the genomic data (e.g., SNPs).\nPedigree-based LMM (none of the individuals in the pedigree are genotyped) and Complete Genomic LMM (all individuals in the pedigree are genotyped) are special cases of Incomplete Genomic LMM (part of the individuals in the pedigree are genotyped)."
},

{
    "location": "manual/workflow/#Get-Data-Ready-1",
    "page": "Workflow",
    "title": "Get Data Ready",
    "category": "section",
    "text": "By default, input data files are comma-separated values (CSV) files, where each line of the file consists of one or more fields, separated by commas. Other field separators such as space (\' \') or tab (\'\\t\') can be used if you supply the keyword argument, e.g, CSV.read(...,delim=\'\\t\') or add_genotypes(...,separator=\'\\t\')Click on the buttons inside the tabbed menu to see the data:<head>\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n<style>\nbody {font-family: Arial;}\n\n/* Style the tab */\n.tab {\n    overflow: hidden;\n    border: 1px solid #ccc;\n    background-color: #f1f1f1;\n}\n\n/* Style the buttons inside the tab */\n.tab button {\n    background-color: inherit;\n    float: left;\n    border: none;\n    outline: none;\n    cursor: pointer;\n    padding: 14px 16px;\n    transition: 0.3s;\n    font-size: 17px;\n}\n\n/* Change background color of buttons on hover */\n.tab button:hover {\n    background-color: #ddd;\n}\n\n/* Create an active/current tablink class */\n.tab button.active {\n    background-color: #ccc;\n}\n\n/* Style the tab content */\n.tabcontent {\n    display: none;\n    padding: 6px 12px;\n    border: 1px solid #ccc;\n    border-top: none;\n}\n</style>\n</head>\n<body>\n\n<div class=\"tab\">\n  <button class=\"tablinks\" onclick=\"openCity(event, \'phenotypes\')\">phenotypes.txt</button>\n  <button class=\"tablinks\" onclick=\"openCity(event, \'pedigree\')\">pedigree.txt</button>\n  <button class=\"tablinks\" onclick=\"openCity(event, \'genotypes\')\">genotypes.txt</button>\n  <button class=\"tablinks\" onclick=\"openCity(event, \'map file\')\">map.txt</button>\n</div>\n\n<div id=\"phenotypes\" class=\"tabcontent\">\n<p>ID,y1,y2,y3,x1,x2,x3,dam</p>\n<p>a1,-0.06,3.58,-1.18,0.9,2,m,0</p>\n<p>a2,-0.6,4.9,0.88,0.3,1,f,0</p>\n<p>a3,-2.07,3.19,0.73,0.7,2,f,0</p>\n<p>a4,-2.63,6.97,-0.83,0.6,1,m,a2</p>\n<p>a5,2.31,3.5,-1.52,0.4,2,m,a2</p>\n<p>a6,0.93,4.87,-0.01,05,2,f,a3</p>\n<p>a7,-0.69,3.1,-1.47,0.5,2,f,a3</p>\n<p>a8,-4.69,7.31,-1.09,0.3,2,m,a6</p>\n<p>a9,-2.81,7.18,0.76,0.4,2,m,a6</p>\n<p>a10,1.92,1.78,-0.88,0.2,1,m,a7</p>\n</div>\n\n<div id=\"pedigree\" class=\"tabcontent\">\n<p>ID,Sire,Dam</p>\n<p>a1,0,0</p>\n<p>a2,0,0</p>\n<p>a3,0,0</p>\n<p>a4,a1,a2</p>\n<p>a5,a1,a2</p>\n<p>a6,a1,a3</p>\n<p>a7,a1,a3</p>\n<p>a8,a4,a6</p>\n<p>a9,a4,a6</p>\n<p>a10,a5,a7</p>\n</div>\n\n<div id=\"genotypes\" class=\"tabcontent\">\n<p>ID,m1,m2,m3,m4,m5</p>\n<p>a1,1,2,1,1,0</p>\n<p>a2,2,1,1,1,1</p>\n<p>a3,1,1,0,1,1</p>\n<p>a4,2,2,0,1,0</p>\n<p>a5,1,1,2,1,1</p>\n<p>a6,2,1,0,0,0</p>\n<p>a7,0,2,1,2,1</p>\n<p>a8,2,2,0,0,0</p>\n<p>a9,2,1,0,1,0</p>\n<p>a10,0,2,2,2,1</p>\n</div>\n\n<div id=\"map file\" class=\"tabcontent\">\n<p>markerID,chromosome,position</p>\n<p>m1,1,16977</p>\n<p>m2,1,434311</p>\n<p>m3,1,1025513</p>\n<p>m4,2,70350</p>\n<p>m5,2,101135</p>\n</div>\n\n\n<script>\nfunction openCity(evt, cityName) {\n    var i, tabcontent, tablinks;\n    tabcontent = document.getElementsByClassName(\"tabcontent\");\n    for (i = 0; i < tabcontent.length; i++) {\n        tabcontent[i].style.display = \"none\";\n    }\n    tablinks = document.getElementsByClassName(\"tablinks\");\n    for (i = 0; i < tablinks.length; i++) {\n        tablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n    }\n    document.getElementById(cityName).style.display = \"block\";\n    evt.currentTarget.className += \" active\";\n}\n</script>\n</body>"
},

{
    "location": "manual/workflow/#Step-1:-Load-Packages-1",
    "page": "Workflow",
    "title": "Step 1: Load Packages",
    "category": "section",
    "text": "using JWAS,CSV,DataFramesThe JWAS package is loaded, as well as the CSV and DataFrame packages for reading text files."
},

{
    "location": "manual/workflow/#Step-2:-Read-data-1",
    "page": "Workflow",
    "title": "Step 2: Read data",
    "category": "section",
    "text": "phenotypes = CSV.read(\"phenotypes.txt\",delim = \',\',header=true)\npedigree   = get_pedigree(\"pedigree.txt\",separator=\",\",header=true)\nhead(phenotypes)output:6×8 DataFrames.DataFrame\n│ Row │ ID │ y1    │ y2   │ y3    │ x1  │ x2 │ x3 │ dam │\n├─────┼────┼───────┼──────┼───────┼─────┼────┼────┼─────┤\n│ 1   │ a1 │ -0.06 │ 3.58 │ -1.18 │ 0.9 │ 2  │ m  │ 0   │\n│ 2   │ a2 │ -0.6  │ 4.9  │ 0.88  │ 0.3 │ 1  │ f  │ 0   │\n│ 3   │ a3 │ -2.07 │ 3.19 │ 0.73  │ 0.7 │ 2  │ f  │ 0   │\n│ 4   │ a4 │ -2.63 │ 6.97 │ -0.83 │ 0.6 │ 1  │ m  │ a2  │\n│ 5   │ a5 │ 2.31  │ 3.5  │ -1.52 │ 0.4 │ 2  │ m  │ a2  │\n│ 6   │ a6 │ 0.93  │ 4.87 │ -0.01 │ 5.0 │ 2  │ f  │ a3  │link to documentation for get_pedigreeThe phenotypic data is read on line 1, and the pedigree data is read on line 2. On line 3, the first several rows of data are shown."
},

{
    "location": "manual/workflow/#Step-3:-Build-Model-Equations-1",
    "page": "Workflow",
    "title": "Step 3: Build Model Equations",
    "category": "section",
    "text": "model_equation = \"y1 = intercept + x1 + x3 + ID + dam\n                  y2 = intercept + x1 + x2 + x3 + ID  \n                  y3 = intercept + x1 + x1*x3 + x2 + ID\"\nmodel=build_model(model_equation, R)link to documentation for build_modelThe non-genomic part of the model equation for a 3-trait analysis is defined on the first 3 lines.The effects fitted in the model for trait y1 are the intercept, x1, x3, direct genetic effects (ID) and maternal genetic effects (dam).\nThe effects fitted in the model for trait y2 are the intercept, x1, x2, x3 and direct genetic effects (ID).\nThe effects fitted in the model for trait y3 are the intercept, x1, the interaction between x1 and x3, x2 and direct genetic effects (ID).On the last line, the model is built given the model equation and residual variance R (a 3x3 matrix). By default, all effects are treated as fixed and classed as factors (categorical variables) rather than covariates (quantitative variables)."
},

{
    "location": "manual/workflow/#Step-4:-Set-Factors-or-Covariate-1",
    "page": "Workflow",
    "title": "Step 4: Set Factors or Covariate",
    "category": "section",
    "text": "set_covariate(model,\"x1\")link to documentation for set_covariateOn line 1, the effect x1 is defined to be a covariate rather than class effect."
},

{
    "location": "manual/workflow/#Step-5:-Set-Random-or-Fixed-Effects-1",
    "page": "Workflow",
    "title": "Step 5: Set Random or Fixed Effects",
    "category": "section",
    "text": "set_random(model,\"x2\",G1)\nset_random(model,\"ID dam\",pedigree,G2)link to documentation for set_randomOn line 1, the x2 class effect is defined as random with variance G1(a 2x2 matrix). On line 2, direct genetic effects and maternal genetic effects are fitted as ID and dam with G2 (a 4x4 matrix) and the inverse of the numerator relationship matrix defined from pedigree."
},

{
    "location": "manual/workflow/#Step-6:-Use-Genomic-Information-1",
    "page": "Workflow",
    "title": "Step 6: Use Genomic Information",
    "category": "section",
    "text": "add_genotypes(model,\"genotypes.txt\",G3,separator=\',\')link to documentation for add_genotypesOn line 1, the genomic part of the model is defined with the genotype file and variance G3 (a 3x3 matrix)."
},

{
    "location": "manual/workflow/#Step-7:-Run-Bayesian-Analysis-1",
    "page": "Workflow",
    "title": "Step 7: Run Bayesian Analysis",
    "category": "section",
    "text": "outputMCMCsamples(model,\"x2\")\nout=runMCMC(model,phenotypes,methods=\"BayesC\",output_samples_frequency=100)link to documentation for outputMCMCsamples\nlink to documentation for runMCMCOn line 1, MCMC samples from runMCMC for x2 is saved to a file, where each row represents one sample from the MCMC. On line 2, a multi-trait BayesC analysis is performed with model and phenotypes as had been defined in step 1-6. MCMC samples for marker effects, location parameters specified on line 1, and all variance components from this analysis are saved every output_samples_frequency iterations to files.Several steps above can be skipped if no related information is available, e.g., step 6 is skipped for pedigree-based LMM. Several detailed examples are available in the examples section. Here is the link to documentation for all Public functions."
},

{
    "location": "manual/workflow/#check-results-1",
    "page": "Workflow",
    "title": "check results",
    "category": "section",
    "text": "Posterior means of location parameters, most variance components, and marker effects are saved in out. They can be listed and obtained askeys(out)\n\n# output:\n#\n# Base.KeyIterator for a Dict{Any,Any} with 7 entries. Keys:\n#   \"Posterior mean of polygenic effects covariance matrix\"\n#   \"Model frequency\"\n#   \"Posterior mean of residual covariance matrix\"\n#   \"Posterior mean of marker effects\"\n#   \"Posterior mean of marker effects covariance matrix\"\n#   \"Posterior mean of location parameters\"\n#   \"Posterior mean of Pi\"\n\nout[\"Posterior mean of residual covariance matrix\"]\n\n# output:\n#\n# 3×3 Array{Float64,2}:\n#   0.674651   -0.103877   0.0834044\n#  -0.103877    0.828135  -0.121798\n#   0.0834044  -0.121798   0.720751\nMCMC samples for marker effects, location parameters specified in step 7, and all variance components are saved to text files in your working directory. They can be obtained asres=readdlm(\"MCMC_samples_marker_effects_y1.txt\",\',\',header=true)"
},

{
    "location": "manual/public/#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "manual/public/#Public-functions-1",
    "page": "Public",
    "title": "Public functions",
    "category": "section",
    "text": "Documentation for JWAS.jl\'s public interface. Below are functions available to general users."
},

{
    "location": "manual/public/#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]\nModules = [JWAS]Pages = [\"public.md\"]\nModules = [JWAS.misc]"
},

{
    "location": "manual/public/#JWAS.build_model",
    "page": "Public",
    "title": "JWAS.build_model",
    "category": "function",
    "text": "build_model(model_equations::AbstractString,R;df::Float64=4.0)\n\nBuild models from model equations with residual varainces R and degree of freedom for residual variance df defaulting to 4.0.\nBy default, all variabels in modelequations are fixed and factors. Set variables to be covariates or random using functions `setcovariate()orset_random()`.\n\n#single-trait\nmodel_equations = \"BW = intercept + age + sex\"\nR               = 6.72\nmodels          = build_model(model_equations,R);\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex\n                   CW = intercept + litter\";\nR               = [6.72   24.84\n                   24.84  708.41]\nmodels          = build_model(model_equations,R);\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.set_covariate",
    "page": "Public",
    "title": "JWAS.set_covariate",
    "category": "function",
    "text": "set_covariate(mme::MME,variables::AbstractString...)\n\nset variables as covariates; mme is the output of function build_model().\n\n#After running build_model, variabels age and year can be set to be covariates as\nset_covariate(models,\"age\",\"year\")\n#or\nset_covariate(models,\"age year\")\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.set_random",
    "page": "Public",
    "title": "JWAS.set_random",
    "category": "function",
    "text": "set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)\n\nset variables as random polygenic effects with pedigree information ped, variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + Age + Animal\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = 1.6\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#single-trait (example 2)\nmodel_equation  = \"y = intercept + Age + Animal + Animal*Age\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = [1.6   0.2\n                   0.2  1.0]\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex + Animal\n                   CW = intercept + age + sex + Animal\"\nmodel           = build_model(model_equations,R);\nped             = get_pedigree(pedfile);\nG               = [6.72   2.84\n                   2.84  8.41]\nset_random(model,\"Animal\", ped,G)\n\n\n\n\n\nset_random(mme::MME,randomStr::AbstractString,G;df=4)\n\nset variables as i.i.d random effects with variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + litter + sex\"\nmodel           = build_model(model_equation,R)\nG               = 0.6\nset_random(model,\"litter\",G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + litter + sex\n                   CW = intercept + litter + sex\"\nmodel           = build_model(model_equations,R);\nG               = [3.72  1.84\n                   1.84  3.41]\nset_random(model,\"litter\",G)\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.get_pedigree",
    "page": "Public",
    "title": "JWAS.get_pedigree",
    "category": "function",
    "text": "get_pedigree(pedfile::AbstractString;header=false,separator=\',\')\n\nGet pedigree informtion from a pedigree file with header defaulting to false and separator defaulting to ,.\nPedigree file format:\n\na,0,0\nc,a,b\nd,a,c\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.add_genotypes",
    "page": "Public",
    "title": "JWAS.add_genotypes",
    "category": "function",
    "text": "add_genotypes(mme::MME,file,G;separator=\' \',header=true,center=true,G_is_marker_variance=false,df=4.0)\n\nGet marker informtion from a genotype file (same order as the phenotype file).\nG defaults to the genetic variance with degree of freedom df=4.0.\nFile format:\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,2,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.add_markers",
    "page": "Public",
    "title": "JWAS.add_markers",
    "category": "function",
    "text": "same to add_genotypes\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.runMCMC",
    "page": "Public",
    "title": "JWAS.runMCMC",
    "category": "function",
    "text": "runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,burnin = 0,starting_value=false,printout_frequency=100,\nmissing_phenotypes=false,constraint=false,methods=\"conventional (no markers)\",output_samples_frequency::Int64 = 0,\nprintout_model_info=true,outputEBV=false)\n\nRun MCMC (marker information included or not) with sampling of variance components.\n\navailable methods include \"conventional (no markers)\", \"RR-BLUP\", \"BayesB\", \"BayesC\".\nsave MCMC samples every outputsamplesfrequency iterations to files output_file defaulting to MCMC_samples.\nthe first burnin iterations are discarded at the beginning of an MCMC run\nPi for single-trait analyses is a number; Pi for multi-trait analyses is a dictionary such as Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1),\nif Pi (Π) is not provided in multi-trait analysis, it will be generated assuming all markers have effects on all traits.\nstarting_value can be provided as a vector for all location parameteres except marker effects.\nprint out the monte carlo mean in REPL with printout_frequency\nconstraint=true if constrain residual covariances between traits to be zeros.\nIndividual EBVs are returned if outputEBV=true.\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.outputMCMCsamples",
    "page": "Public",
    "title": "JWAS.outputMCMCsamples",
    "category": "function",
    "text": "outputMCMCsamples(mme::MME,trmStr::AbstractString...)\n\nGet MCMC samples for specific location parameters.\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.showMME",
    "page": "Public",
    "title": "JWAS.showMME",
    "category": "function",
    "text": "showMME(mme::MME,df::DataFrame)\n\nShow left-hand side and right-hand side of mixed model equations (no markers).\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.solve",
    "page": "Public",
    "title": "JWAS.solve",
    "category": "function",
    "text": "solve(mme::MME,df::DataFrame;solver=\"default\",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)\n\nSolve the mixed model equations (no marker information) without estimating variance components.\n\nAvailable solvers includes default,Jacobi,GaussSeidel,Gibbs sampler.\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.GWAS",
    "page": "Public",
    "title": "JWAS.misc.GWAS",
    "category": "function",
    "text": "GWAS(marker_effects_file;header=false)\n\nCompute the model frequency for each marker (the probability the marker is included in the model) using samples of marker effects stored in markereffectsfile.\n\n\n\n\n\nGWAS(marker_effects_file,map_file,model;header=false,window_size=\"1 Mb\",threshold=0.001)\n\nrun genomic window-based GWAS\n\nMCMC samples of marker effects are stored in markereffectsfile\nmap_file has the marker position information\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.get_additive_genetic_variances",
    "page": "Public",
    "title": "JWAS.misc.get_additive_genetic_variances",
    "category": "function",
    "text": "get_additive_genetic_variances(model::MME,files...;header=true)\n\nGet MCMC samples for additive genetic variances using samples for marker effects stored in files.\nReturn a vector for single-trait analysis and an array of matrices for multi-trait analysis\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.get_heritability",
    "page": "Public",
    "title": "JWAS.misc.get_heritability",
    "category": "function",
    "text": "get_heritability(samples_for_genetic_variances::Array{Array{Float64,2},1},samples_for_residual_variances::Array{Array{Float64,2},1}))\n\nGet MCMC samples for heritabilities using MCMC samples for genetic variances and residual variances for multi-trait analysis\n\n\n\n\n\nget_heritability(samples_for_genetic_variances::Array{Float64,1},samples_for_residual_variances::Array{Float64,1}))\n\nGet MCMC samples for heritabilities using MCMC samples for genetic variances and residual variances for single-trait analysis\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.get_correlations",
    "page": "Public",
    "title": "JWAS.misc.get_correlations",
    "category": "function",
    "text": "get_correlations(samples_for_genetic_variances::Array{Array{Float64,2},1})\n\nGet MCMC samples for correlations using MCMC samples for covariance matrces\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.get_breeding_values",
    "page": "Public",
    "title": "JWAS.misc.get_breeding_values",
    "category": "function",
    "text": "get_breeding_values(model)\n\nGet esitimated breeding values and prediction error variances using samples of marker effects stored in files   for individuals defined by outputEBV(model,IDs::Array{String,1}), defaulting to all phenotyped individuals.\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.reformat",
    "page": "Public",
    "title": "JWAS.misc.reformat",
    "category": "function",
    "text": "reformat(G::Array{Float64,2},ntraits=false)\n\nconvert the format from a Matrix to a Array of Matrices, when the matrix is\n\nread from text file storing samples for covariance matrices or heritabilities\n\nthe size of the Matrix is\nntraits x nsamples for MCMC samples for heritbilities\n(ntraits^2) x nsamples for MCMC samples for covariance matrices\nthe Array of Matrices is an array (length=nsamples) of\nmatrices of size ntraits X ntraits for covariance matrices\nvectors of length ntraits for heritabilities\n\n\n\n\n\nreformat(G::Array{Array{Float64,2},1},1})\n\nconvert the format from a Array of Matrices to a Matrix, then the Matrix chance\n\nbe write to text files to save samples for covariance matrices or heritabilities\n\nthe Array of Matrices is an array (length=nsamples) of\nmatrices of size ntraits X ntraits for covariance matrices (Array{Array{Float64,2},1})\nvectors of length ntraits for heritabilities (Array{Array{Float64,1},1})\nthe size of the Matrix is\nntraits x nsamples for MCMC samples for heritbilities\n(ntraits^2) x nsamples for MCMC samples for covariance matrices\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.report",
    "page": "Public",
    "title": "JWAS.misc.report",
    "category": "function",
    "text": "report(X::Array{Array{Float64,2},1};index=false)\n\nshow summary statistics for MCMC samples (matrices or index [i,j] of matrices)\n\n\n\n\n\nreport(X::Array{Array{Float64,1},1};index=false)\n\nshow summary statistics for MCMC samples (vectors or index i of vectors)\n\n\n\n\n\n"
},

{
    "location": "manual/public/#JWAS.misc.QC",
    "page": "Public",
    "title": "JWAS.misc.QC",
    "category": "function",
    "text": "QC(infile,outfile;separator=\' \',header=true,missing=false,MAF=0.1)\n\nQuality control for input file infile, then write out to output file: outfile.\nDelete loci with minor allele frequency < MAF.\nmissing genotypes are replaced by column means.\nFile format (header=true,separator=\',\',missing=9):\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,9,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n\n\n"
},

{
    "location": "manual/public/#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": "build_model\nset_covariate\nset_random\nget_pedigree\nadd_genotypes\nadd_markers\nrunMCMC\noutputMCMCsamples\nshowMME\nsolveJWAS.misc.GWAS\nJWAS.misc.get_additive_genetic_variances\nJWAS.misc.get_heritability\nJWAS.misc.get_correlations\nJWAS.misc.get_breeding_values\nJWAS.misc.reformat\nJWAS.misc.report\nJWAS.misc.QC"
},

{
    "location": "manual/internals/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "manual/internals/#Internal-functions-1",
    "page": "Internals",
    "title": "Internal functions",
    "category": "section",
    "text": "Documentation for JWAS.jl\'s internal (private) interface, which are not available to general users. These internal functions are small blocks that public function build on.<!–-"
},

{
    "location": "manual/internals/#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internals.md\"]"
},

{
    "location": "manual/internals/#JWAS.add_genotypes-Tuple{JWAS.MME,Any,Any}",
    "page": "Internals",
    "title": "JWAS.add_genotypes",
    "category": "method",
    "text": "add_genotypes(mme::MME,file,G;separator=\' \',header=true,center=true,G_is_marker_variance=false,df=4.0)\n\nGet marker informtion from a genotype file (same order as the phenotype file).\nG defaults to the genetic variance with degree of freedom df=4.0.\nFile format:\n\nAnimal,marker1,marker2,marker3,marker4,marker5\nS1,1,0,1,1,1\nD1,2,0,2,2,1\nO1,1,2,0,1,0\nO3,0,0,2,1,1\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.add_markers-Tuple{JWAS.MME,Any,Any}",
    "page": "Internals",
    "title": "JWAS.add_markers",
    "category": "method",
    "text": "same to add_genotypes\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.build_model-Tuple{AbstractString,Any}",
    "page": "Internals",
    "title": "JWAS.build_model",
    "category": "method",
    "text": "build_model(model_equations::AbstractString,R;df::Float64=4.0)\n\nBuild models from model equations with residual varainces R and degree of freedom for residual variance df defaulting to 4.0.\nBy default, all variabels in modelequations are fixed and factors. Set variables to be covariates or random using functions `setcovariate()orset_random()`.\n\n#single-trait\nmodel_equations = \"BW = intercept + age + sex\"\nR               = 6.72\nmodels          = build_model(model_equations,R);\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex\n                   CW = intercept + litter\";\nR               = [6.72   24.84\n                   24.84  708.41]\nmodels          = build_model(model_equations,R);\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.get_pedigree-Tuple{AbstractString}",
    "page": "Internals",
    "title": "JWAS.get_pedigree",
    "category": "method",
    "text": "get_pedigree(pedfile::AbstractString;header=false,separator=\',\')\n\nGet pedigree informtion from a pedigree file with header defaulting to false and separator defaulting to ,.\nPedigree file format:\n\na,0,0\nc,a,b\nd,a,c\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.outputEBV-Tuple{Any,Any}",
    "page": "Internals",
    "title": "JWAS.outputEBV",
    "category": "method",
    "text": "outputEBV(model,IDs::Array;PEV=false)\n\nOutput estimated breeding values and prediction error variances (defaulting to false) for IDs.\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.outputMCMCsamples-Tuple{JWAS.MME,Vararg{AbstractString,N} where N}",
    "page": "Internals",
    "title": "JWAS.outputMCMCsamples",
    "category": "method",
    "text": "outputMCMCsamples(mme::MME,trmStr::AbstractString...)\n\nGet MCMC samples for specific location parameters.\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.runMCMC-Tuple{Any,Any}",
    "page": "Internals",
    "title": "JWAS.runMCMC",
    "category": "method",
    "text": "runMCMC(mme,df;Pi=0.0,estimatePi=false,chain_length=1000,burnin = 0,starting_value=false,printout_frequency=100,\nmissing_phenotypes=false,constraint=false,methods=\"conventional (no markers)\",output_samples_frequency::Int64 = 0,\nprintout_model_info=true,outputEBV=false)\n\nRun MCMC (marker information included or not) with sampling of variance components.\n\navailable methods include \"conventional (no markers)\", \"RR-BLUP\", \"BayesB\", \"BayesC\".\nsave MCMC samples every outputsamplesfrequency iterations to files output_file defaulting to MCMC_samples.\nthe first burnin iterations are discarded at the beginning of an MCMC run\nPi for single-trait analyses is a number; Pi for multi-trait analyses is a dictionary such as Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1),\nif Pi (Π) is not provided in multi-trait analysis, it will be generated assuming all markers have effects on all traits.\nstarting_value can be provided as a vector for all location parameteres except marker effects.\nprint out the monte carlo mean in REPL with printout_frequency\nconstraint=true if constrain residual covariances between traits to be zeros.\nIndividual EBVs are returned if outputEBV=true.\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.set_covariate-Tuple{JWAS.MME,Vararg{AbstractString,N} where N}",
    "page": "Internals",
    "title": "JWAS.set_covariate",
    "category": "method",
    "text": "set_covariate(mme::MME,variables::AbstractString...)\n\nset variables as covariates; mme is the output of function build_model().\n\n#After running build_model, variabels age and year can be set to be covariates as\nset_covariate(models,\"age\",\"year\")\n#or\nset_covariate(models,\"age year\")\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.set_random-Tuple{JWAS.MME,AbstractString,Any}",
    "page": "Internals",
    "title": "JWAS.set_random",
    "category": "method",
    "text": "set_random(mme::MME,randomStr::AbstractString,G;df=4)\n\nset variables as i.i.d random effects with variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + litter + sex\"\nmodel           = build_model(model_equation,R)\nG               = 0.6\nset_random(model,\"litter\",G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + litter + sex\n                   CW = intercept + litter + sex\"\nmodel           = build_model(model_equations,R);\nG               = [3.72  1.84\n                   1.84  3.41]\nset_random(model,\"litter\",G)\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.set_random-Tuple{JWAS.MME,AbstractString,JWAS.PedModule.Pedigree,Any}",
    "page": "Internals",
    "title": "JWAS.set_random",
    "category": "method",
    "text": "set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)\n\nset variables as random polygenic effects with pedigree information ped, variances G whose degree of freedom df defaults to 4.0.\n\n#single-trait (example 1)\nmodel_equation  = \"y = intercept + Age + Animal\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = 1.6\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#single-trait (example 2)\nmodel_equation  = \"y = intercept + Age + Animal + Animal*Age\"\nmodel           = build_model(model_equation,R)\nped             = get_pedigree(pedfile)\nG               = [1.6   0.2\n                   0.2  1.0]\nset_random(model,\"Animal Animal*Age\", ped,G)\n\n#multi-trait\nmodel_equations = \"BW = intercept + age + sex + Animal\n                   CW = intercept + age + sex + Animal\"\nmodel           = build_model(model_equations,R);\nped             = get_pedigree(pedfile);\nG               = [6.72   2.84\n                   2.84  8.41]\nset_random(model,\"Animal\", ped,G)\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.showMME-Tuple{JWAS.MME,DataFrames.DataFrame}",
    "page": "Internals",
    "title": "JWAS.showMME",
    "category": "method",
    "text": "showMME(mme::MME,df::DataFrame)\n\nShow left-hand side and right-hand side of mixed model equations (no markers).\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.solve-Tuple{JWAS.MME,DataFrames.DataFrame}",
    "page": "Internals",
    "title": "JWAS.solve",
    "category": "method",
    "text": "solve(mme::MME,df::DataFrame;solver=\"default\",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)\n\nSolve the mixed model equations (no marker information) without estimating variance components.\n\nAvailable solvers includes default,Jacobi,GaussSeidel,Gibbs sampler.\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.add_term-Tuple{Any,AbstractString}",
    "page": "Internals",
    "title": "JWAS.add_term",
    "category": "method",
    "text": "add to model an extra term: imputation_residual\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.getMME-Tuple{JWAS.MME,DataFrames.DataFrame}",
    "page": "Internals",
    "title": "JWAS.getMME",
    "category": "method",
    "text": "Construct mixed model equations with\n\nincidence matrix: X      ; response        : ySparse; left-hand side  : mmeLhs ; right-hand side : mmeLhs ;\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#JWAS.get_outputX_others-Tuple{Any,Any}",
    "page": "Internals",
    "title": "JWAS.get_outputX_others",
    "category": "method",
    "text": "get_outputX_others(model)\n\nMake incidence matrices for effects involve in EBV inclung J, ϵ, pedTrmVec except marker covariates\n\n\n\n\n\n"
},

{
    "location": "manual/internals/#Internal-interface-1",
    "page": "Internals",
    "title": "Internal interface",
    "category": "section",
    "text": "Modules = [JWAS,JWAS.PedModule]\nOrder   = [:function, :type]–>"
},

{
    "location": "examples/examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "link to Jupyter Notebook for this example"
},

{
    "location": "examples/examples/#Bayesian-Linear-Mixed-Models-(conventional)-1",
    "page": "Examples",
    "title": "Bayesian Linear Mixed Models (conventional)",
    "category": "section",
    "text": "Univariate Linear Mixed Model (conventional)\nMultivariate Linear Mixed Model (conventional)"
},

{
    "location": "examples/examples/#Bayesian-Linear-Additive-Genetic-Model-1",
    "page": "Examples",
    "title": "Bayesian Linear Additive Genetic Model",
    "category": "section",
    "text": "Univariate Linear Additive Genetic Model\nMultivariate Linear Additive Genetic Model"
},

{
    "location": "examples/examples/#Bayesian-Linear-Mixed-Models-(Genomic-Data)-1",
    "page": "Examples",
    "title": "Bayesian Linear Mixed Models (Genomic Data)",
    "category": "section",
    "text": "Univariate Linear Mixed Model (Genomic data)\nMultivariate Linear Mixed Model (Genomic data)"
},

{
    "location": "examples/examples/#Single-step-Bayesian-Linear-Mixed-Models-(Genomic-Data)-1",
    "page": "Examples",
    "title": "Single-step Bayesian Linear Mixed Models (Genomic Data)",
    "category": "section",
    "text": "Univariate Single-step Bayesian Linear Mixed Models (Genomic Data)\nMultivariate Single-step Bayesian Linear Mixed Models (Genomic Data)"
},

]}
