Models
======

The starting point for these methods is a mixed linear model of the form:

.. math::

  \boldsymbol{y=X\beta+Z\alpha+e}

where  :math:`y` is an :math:`n\times 1` vector of trait
phenotypic values, :math:`\boldsymbol{X}` is an :math:`n\times p` incidence matrix relating
the vector :math:`\boldsymbol{\beta}` of non-genetic fixed effects to :math:`\boldsymbol{y}`,  :math:`\boldsymbol{Z}`
is an :math:`n\times k` matrix of genotype covariates (coded as 0, 1 or 2)
for :math:`k` SNP markers, :math:`\boldsymbol{\alpha}` is a :math:`k\times 1` vector of random
partial regression coefficients of the :math:`k` SNPs (which are more
commonly referred to as the marker effects), and :math:`\boldsymbol{e}` is a
vector of residuals. 

To proceed with Bayesian regression, prior distributions must be
specified for :math:`\beta`, :math:`\alpha` and :math:`e`. In all the models
considered here a flat prior is used for
:math:`\beta`, and conditional on the residual variance, :math:`\sigma^2_e`, a
normal distribution with null mean and covariance matrix
:math:`\sigma^2_e` is used for the vector of residuals, where :math:`R`
is a diagonal matrix. Further, :math:`\sigma^2_e` is treated as an unknown
with a scaled inverse chi-square prior. The alternative methods differ 
only in the prior used for :math:`\alpha`.



BayesB
^^^^^^

description of model

BayesCPi
^^^^^^^^

description of model

BayesCPi Dominant
^^^^^^^^^^^^^^^^^

description of model

BayesR
^^^^^^

Coming soon!

BayesN
^^^^^^

Coming soon!
