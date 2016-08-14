#now assume genotype and phenotype have same order
#genotype: file with ids || Array{Float64,2} wihtout IDs
#phenotype: file with ids || Array{Float64,1} wihtout IDs
function runBR(input;genotype="genofile",phenotype="phenofile")
 
      srand(input.seed)

      geno  = make_genotypes(genotype,center=input.centering)
      y     = make_yVecs(phenotype) #in tool.jl
      fixed = QTL.FixedMatrix(ones(geno.nObs,1),[0]); #modify later

      if input.method=="BayesC0" && input.estimateVariance==false
        out=BayesC0_constantvariance(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesC0"
        out=BayesC0(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesB"
        out=BayesB(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesC" && input.estimateVariance==false
        out=BayesC_constantvariance(y,geno,fixed,input,outFreq=input.outFreq)
      elseif input.method=="BayesC"
        out=BayesC(y,geno,fixed,input,outFreq=input.outFreq)
      end
      return out
end
