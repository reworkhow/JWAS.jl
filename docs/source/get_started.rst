Get Started
===========

An example to run BayesB:

.. code-block:: julia

	using GenSel
	
	myOption=Dict()
	myOption["run"]          = "BayesB"
	myOption["seed"]         = 10	
	myOption["chainLength"]  = 2000
	myOption["probFixed"]    = 0.5 
	myOption["estimatePi"]   = "no"
	output = runJWAS(parm,X,y)
	
To read in X, y,


