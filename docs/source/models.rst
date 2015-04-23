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

BayesA
^^^^^^

The prior assumption is that marker effects have identical
and independent univariate-t distributions each with a null mean,
scale parameter :math:`S^2_{\alpha}` and :math:`\nu` degrees of freedom.
This is equivalent to assuming that the marker effect at locus :math:`i` has a univariate normal
with null mean and unknown, locus-specific variance :math:`\sigma^2_i`,
which in turn is assigned a scaled inverse chi-square prior with scale
parameter :math:`S^2_{\alpha}` and :math:`\nu_{\alpha}` degrees of freedom. 

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
