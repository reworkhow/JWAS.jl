#allow intermediate omics data in neural networks (8/31/2021)

1. add an argument "latent_traits" in build_model() to allow users to provide
   which columns in the phenotypic data will be used as the observed values for
   intermediate traits. Add latent_traits as a member in the struct mme.

2. Because multi-trait models with latent traits as responses will be used in neural
   network (between input layer and hidden layer), we reassign column names of latent
   traits to mme.lhsVec, which will be used to make the matrices in the multi-trait models.
   In this case, we still need a variable to save the (empirical) phenotypes (i.e., output),
   so mme.yobs is made to save it. Add yobs as a member in the struct mme.

3. mme.ySparse is used to same values for latent traits. In intermediate omics data,
   because some/many elements in mme.ySparse are observed, only missing values in mme.ySparse
   are sampled.

# add neural network
