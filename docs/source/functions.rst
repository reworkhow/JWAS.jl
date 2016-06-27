Functions
===========

An example to run BayesC is shown below. Let's start by simulating a dataset (a more complicated simulation
package ``XSim`` can be found at `QTL.rocks <http://QTL.rocks>`_.) Here a naive simulation is performed using Distributions.jl.

.. code-block:: julia

	using(Distributions)
	d = Binomial(2,0.5)

	nObs     = 10
	nMarkers = 100
	X        = float(rand(d,(nObs,nMarkers)))
	α        = randn(nMarkers)
	a        = X*α
	stdGen   = std(a)
	a        = a/stdGen
	y        = a + randn(nObs)

Though ``JWAS`` can fit any fixed effects and weighted phenotypes, this example below showed how to run BayesC with only
population mean as fixed effects.

.. code-block:: julia

	using JWAS

	myOption=Dict()
	myOption["run"]           = "BayesC"
	myOption["seed"]          = 314
	myOption["chainLength"]   = 5000
	myOption["probFixed"]     = 0.95
	myOption["estimatePi"]    = "yes"
	myOption["estimateScale"] = "yes"
	myOption["varGenotypic"]  = 1
	myOption["varResidual"]   = 1

	output = runJWAS(myOption,X,y)

Posterior samples for all parameters of interest are saved in the dictionary ``output``. A plot of the result is shown below.

####MTJWAS.buildModel
.. code-block:: julia

	buildModel(model_equations,R)

Arguments

* _**model_equations::AbstractString**_: model equations
* _**R::Array{Float64,2}**_: residual covariance matrix

Result
* _**::MME**_:

Example
```julia
#string for model equations
model_equations = "BW = intercept + age + sex;
                  CW = intercept + age + sex";

#residual covariance matrix
R=[6.72   24.84
   24.84  708.41]
#model building
models = buildModel(model_equations,R);
```

***
####MTJWAS.setAsCovariate
```julia
setAsCovariate(models,covStr)
```

Arguments

* _**models::MME**_: When models:MME are made with _**buildModel()**_, all variables are default to of type factors.
* _**covStr::AbstractString**_: variable set to be covariates.

Example
```julia
#set as covariates
setAsCovariate(models,"age")
```

***
####MTJWAS.setAsRandom
```julia
setAsRandom(models,randomStr,ped,G)
```

Arguments

* _**models::MME**_: When models:MME are made with _**buildModel()**_, all variables are default to fixed effects.
* _**randomStr::AbstractString**_: variable set to be random.
* **_ped_**::Pedigree:
* _**G**_:: covariace matrix

Example
```julia
pedfile="pedigree.txt"
model_equations = "BW = intercept + age + sex + Animal;
                   CW = intercept + age + sex + Animal";
model           = buildModel(model_equations,R);
ped             = PedModule.mkPed(pedfile);

setAsRandom(model2,"Animal", ped,G)
```

***
####MTJWAS.addMarkers
```julia
addMarkers(models,file,G,[keyword options])
```

Arguments

* _**models::MME**_: When models:MME are made with _**buildModel()**_, no marker information are provided.
* _**file::AbstractString**_: genotype file
* _**G::Array{Float64,2}**_: genetic variances explained by markers

Keyword Arguments

* **_separator::Char_** : Assume that fields are split by the **_separator_** character, default to ' '.
* **_header::Bool_**: Use the information from the file's header line to determine Marker IDs. Defaults to true.


Example
```julia
addMarkers(models,genofile,G,separator=',',header=false);
```




***
####<font,color="red"> MTJWAS.runMCMC</font>
```julia
runMCMC(models,data,[keyword options])
```

Arguments

* _**data::DataFrames**_:
* _**models::MME**_:


Keyword Arguments

* **_chain_length::Int64_** : the length of MCMC chain, default to 1000,
* **_starting_value::Array{Float64,1}_**: starting values for samples for location parametres, default to zeros.
* _**printout_frequency::Int64**_: print out the monte carlo mean with printout_frequency, default to 100.
* **_thin::Int64_**: save samples of marker effects every _**thin**_ iterations to files with filename _**output_files**_, default to 100.
* _**output_files::AbstractString**_: save samples of marker effects every _**thin**_ iterations to files with filename _**output_files**_, default to "marker_effects".
