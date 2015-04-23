Get Started
===========

An example to run BayesB:

	using GenSel
	
	myOption=Dict()
	myOption["run"]          = "BayesB"
	myOption["seed"]         = 10	
	myOption["chainLength"]  = 2000
	myOption["probFixed"]    = 0.5 
	myOption["estimatePi"]   = "no"
	myOption["dfEffectVar"]  = 4
	myOption["nuRes"]        = 4
	myOption["varGenotypic"] = 1  
	myOption["varResidual"]  = 1
	output = runJWAS(parm,X,y)
	
To read in X, y,


