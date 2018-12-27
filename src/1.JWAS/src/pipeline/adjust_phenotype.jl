#result 'y-intercept-group'
function adjust_phenotypes(model,phenotype;chain_length=50000,printout_frequency=10000)
  output=runMCMC(model,phenotype,chain_length=chain_length,printout_frequency=printout_frequency);
  animals_count= model.modelTermDict["1:Animal"].nLevels
  fixed_index  = 1:size(model.X,2)-animals_count
  fixed_X      = model.X[:,fixed_index]
  fixed_effects= map(Float64,output["Posterior mean of location parameters"][:,2][fixed_index])
  adjusted_y= vec(model.ySparse-fixed_X*fixed_effects)
  adjusted_phenotype= DataFrame(Animal=phenotype[:Animal],adjusted_phenotype=adjusted_y)
  return adjusted_phenotype
end

#get_genetic_variance
