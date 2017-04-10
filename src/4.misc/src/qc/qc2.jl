
#remove problematic loci: 1) minor allele frequence 2)fixed

maf=vec(mean(M,1)/2);
select1=0.01.<maf.<0.99;
fixed=vec(var(M,1));
select2 = fixed.!=0;
select = select1 & select2
M=M[:,select];


# show distritbuion of allele frequencies




#covert 3 parameters: 1)heritability; 2)genetic variance; 3)residual variance
function H2GR(;heritability=false,genetic_variance=false,residual_variance=false)
    if heritability == false
        heritability=genetic_variance/(genetic_variance+residual_variance)
    end
    if genetic_variance == false
        genetic_variance = residual_variance*heritability/(1-heritability)
    end
    if residual_variance == false
        residual_variance= genetic_variance/heritability- genetic_variance
    end
    println("heritability      = ",heritability)
    println("genetic variance  = ",genetic_variance)
    println("residual variance = ",residual_variance)
end


#Python magic (load script) (bash)
