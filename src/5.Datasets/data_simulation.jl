
using Random,Statistics,DelimitedFiles,XSim,DataFrames
Random.seed!(2);

#set genome information
chrLength= 1.0  #length of each chromosome 
numChr   = 1    #number of chromosomes
nLoci    = 1100 #number of loci for each chromosome
nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
build_genome(numChr,chrLength,nLoci,nQTL) #this genome information will be used for subsequent computaions

#generate founders
popSizeFounder = 500
sires = sampleFounders(popSizeFounder);
dams  = sampleFounders(popSizeFounder);

ngen,popSize = 100,500
sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);

sires2,dams2,gen2 = sampleRan(100, 1, sires1, dams1);
sires3,dams3,gen2 = sampleRan(100, 1, sires2, dams2);
sires4,dams4,gen2 = sampleRan(100, 1, sires3, dams3);

animals=concatCohorts(sires2,dams2,sires3,dams3,sires4,dams4);

M = getOurGenotypes(animals);

outputPedigree(animals, "output.txt")

writedlm("pedigree.txt",map(Int64,readdlm("output.txt.ped")),',')

[reshape(["ID";"snp".*string.(1:1000)],1,1001); map(Int64,readdlm("output.txt.gen"))[:,1:1001]]

writedlm("genotypes.txt",[reshape(["ID";"snp".*string.(1:1000)],1,1001); map(Int64,readdlm("output.txt.gen"))[:,1:1001]],',') #the last 100 SNP not included in marker panel (extra polygenic effect)

ID=map(string,map(Int64,readdlm("output.txt.gen"))[:,1]);

["ID";"snp".*string.(1:1000)]

pos = collect(100:88:1100)
BV  = [M[:,pos]*randn(length(pos)); M[:,pos]*randn(length(pos)); M[:,pos]*randn(length(pos))]
for i =1:3
    BV[(i-1)*100+1 : i*100]= BV[(i-1)*100+1 : i*100]/sqrt(var(BV[(i-1)*100+1 : i*100]))
end

c_f1 = rand(rand(100),300)     #fixed covariate 0.5
f_r1 = rand(["1","2","3"],300) #random factor
f_f1 = rand(1:2,300);

df=DataFrame(ID=ID,y1=randn(300),y2=randn(300),y3=randn(300),x1=c_f1,x2=f_f1,x3=f_r1);
first(df)

using JWAS

model_equation="y1= intercept + x1
                y2= intercept + x2 + x2*x3
                y3= intercept + x1 + x2"

model=build_model(model_equation);
set_covariate(model,"x1");
set_random(model,"x2")
set_random(model,"x2*x3")

JWAS.getMME(model, df)
res=[0.0,0.1,0.0,0.1,0.2,randn(),randn(),randn(),randn(),randn(),randn(),0.0,0.2,0.1,0.2]
noBV=model.X*res

pheno=BV+noBV+[randn(300);randn(300);randn(300)];

df[!,:y1]=pheno[1:300]
df[!,:y2]=pheno[301:600]
df[!,:y3]=pheno[601:900]

using CSV

?CSV.write

CSV.write("phenotypes.txt",df)


