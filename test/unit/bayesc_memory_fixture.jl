using Test, JWAS, DataFrames, CSV, Random, Statistics, LinearAlgebra

function bayesc_memory_case(;traits=1,sampler=:I,precision=Float32,mode="ordinary",
    seed=73,annotated=false,copying=false,upstream=false,
    n=36,t=24,p=17,outer=200,burnin=50,freq=10,folder=nothing,prediction_all=false)
    rng=MersenneTwister(20261031)
    X=precision.(rand(rng,0:2,n,p))
    ids=["animal_$(lpad(i,4,'0'))" for i in 1:n]
    markers=["rs_$(1000+7*i)" for i in 1:p]
    y1=0.3 .* Float64.(X[:,1]) .+ randn(rng,n)
    y2=0.2 .* Float64.(X[:,2]) .+ 0.4 .* y1 .+ randn(rng,n)
    annotation=annotated ? Float64.(hcat(isodd.(1:p),mod.(1:p,3).==0)) : false
    pi=traits==1 ? 0.65 : Dict([0.0,0.0]=>0.45,[1.0,0.0]=>0.20,
                                  [0.0,1.0]=>0.15,[1.0,1.0]=>0.20)
    G=traits==1 ? 1.0 : [1.0 0.2;0.2 1.0]
    R=traits==1 ? 1.0 : [1.0 0.1;0.1 1.0]
    ownership=(!upstream && copying) ? (;copy_genotypes=true) : (;)
    alignment=(!upstream && copying) ? (;shared_genotypes=false) : (;)
    global bayesc_memory_genotypes=get_genotypes(X,G;method="BayesC",Pi=pi,
        quality_control=false,double_precision=(precision===Float64),
        annotations=annotation,multi_trait_sampler=sampler,ownership...)
    g=bayesc_memory_genotypes
    g.obsID .= ids; g.markerID .= markers
    owned=g.genotypes
    ph=DataFrame(ID=ids,y1=Union{Missing,Float64}[i<=t ? y1[i] : missing for i in 1:n],
        FARM=ifelse.(iseven.(1:n),"A","B"),LITTER=["L$(mod1(i,6))" for i in 1:n])
    if traits==2
        # Row 6 has one observed trait and one missing trait, but stays in training.
        ph.y2=Union{Missing,Float64}[i<=t && i!=6 ? y2[i] : missing for i in 1:n]
    end
    equation=traits==1 ? "y1 = intercept + FARM + LITTER + bayesc_memory_genotypes" :
        "y1 = intercept + FARM + LITTER + bayesc_memory_genotypes\ny2 = intercept + FARM + LITTER + bayesc_memory_genotypes"
    model=build_model(equation,R)
    set_random(model,"LITTER",traits==1 ? 0.5 : [0.5 0.1;0.1 0.5])
    prediction_rows=prediction_all ? (1:n) : (t+1:n)
    outputEBV(model,ids[prediction_rows])
    fast=mode=="ordinary" ? false : 5
    outdir=folder===nothing ? joinpath(mktempdir(),"mcmc") : folder
    output=runMCMC(model,ph;chain_length=(fast===false ? outer : outer*fast),
        burnin,output_samples_frequency=freq,seed,fast_blocks=fast,
        independent_blocks=(mode=="independent"),double_precision=(precision===Float64),
        output_folder=outdir,output_heritability=false,printout_model_info=false,
        printout_frequency=outer+1,alignment...)
    animals=DataFrame[]
    for trait in 1:traits
        frame=CSV.read(joinpath(outdir,"EBV_y$(trait).txt"),DataFrame;types=Dict(:ID=>String))
        frame.Trait=fill("y$(trait)",nrow(frame));push!(animals,frame)
    end
    ebv=vcat(animals...)
    marker=output["marker effects bayesc_memory_genotypes"]
    @test nrow(ebv)==traits*length(prediction_rows)
    @test Set(zip(String.(ebv.Trait),String.(ebv.ID)))==
        Set(("y$(k)",id) for k in 1:traits for id in ids[prediction_rows])
    @test nrow(marker)==traits*p
    @test Set(zip(string.(marker.Trait),string.(marker.Marker_ID)))==
        Set(("y$(k)",id) for k in 1:traits for id in markers)
    @test all(isfinite,ebv.EBV) && all(isfinite,ebv.PEV)
    @test length(readlines(joinpath(outdir,"MCMC_samples_residual_variance.txt")))-1==fld(outer-burnin,freq)
    storage=upstream ? nothing : JWAS.genotype_storage_report(g)
    if !upstream && !copying
        @test owned === X
        @test storage.backing_matrix_count==1
        @test storage.backing_matrix_bytes==sizeof(X)
        @test storage.output_shares_training_parent && storage.all_blocks_share_training
        @test g.genotypes isa StridedMatrix{precision}
    end
    (;ebv,marker,storage,model,g,outdir)
end

function compare_bayesc_memory(a,b)
    @test nrow(a.ebv)==nrow(b.ebv)
    @test Set(zip(a.ebv.Trait,a.ebv.ID))==Set(zip(b.ebv.Trait,b.ebv.ID))
    animals=innerjoin(a.ebv,b.ebv;on=[:Trait,:ID],makeunique=true,validate=(true,true))
    @test nrow(animals)==nrow(a.ebv)
    @test isapprox(animals.EBV,animals.EBV_1;atol=1e-4,rtol=1e-4)
    @test isapprox(animals.PEV,animals.PEV_1;atol=1e-4,rtol=1e-4)
    x,y=copy(a.marker),copy(b.marker)
    x.Marker_ID=string.(x.Marker_ID);y.Marker_ID=string.(y.Marker_ID)
    @test Set(zip(x.Trait,x.Marker_ID))==Set(zip(y.Trait,y.Marker_ID))
    markers=innerjoin(x,y;on=[:Trait,:Marker_ID],makeunique=true,validate=(true,true))
    @test nrow(markers)==nrow(x)
    for name in (:Estimate,:SD,:Model_Frequency)
        @test isapprox(markers[!,name],markers[!,Symbol(string(name)*"_1")];atol=1e-4,rtol=1e-4)
    end
    (;animals,markers)
end
