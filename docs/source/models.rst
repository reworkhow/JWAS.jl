Models
======

The starting point for these methods is a mixed linear model of the form:

.. math::

  \boldsymbol{y=X\beta+Z\alpha+e}

where :math:`\boldsymbol{y}` is an nx1 vector of trait phenotypic values,  :math:`\boldsymbol{X}` is an nxp incidence matrix relating 
the vector  :math:`boldsymbol{\beta}` of non-genetic fixed effects to  

where :math:`boldsymbol{y}` is an :math:`n\times 1` vector of trait
phenotypic values, :math:`\boldsymbol{X}` is an :math:`n\times p` incidence matrix relating
the vector :math:`\boldsymbol{\beta}` of non-genetic fixed effects to :math:`\boldsymbol{y}`, :math:`\boldsymbol{Z}`
is an :math:`n\times k` matrix of genotype covariates (coded as 0, 1 or 2)
for :math:`k` SNP markers, :math:`\boldsymbol{\alpha}` is a :math:`k\times 1` vector of random
partial regression coefficients of the :math:`k` SNPs (which are more
commonly referred to as the marker effects), and :math:`\boldsymbol{e}` is a
vector of residuals. 




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
