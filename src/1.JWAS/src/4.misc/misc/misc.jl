module misc

using DataFrames,CSV,DelimitedFiles,ProgressMeter,Statistics

include("convergence_diagnosis.jl")

function lsmeans(model)
    println("The reference grid (no interatction in the model) is")
    grid=[]
    output = model.output["location parameters"]
    for term in model.modelTerms
        if term.random_type == "fixed"
            tratii = term.iTrait
            effecti = split(term.trmStr,":")[2]
            print("trait ",tratii,", ",effecti," = ")
            estimate = output[(output[!,:Trait] .== tratii) .& (output[!,:Effect] .== effecti),:Estimate]
            if term.nLevels != 1
                push!(grid,[tratii,effecti,term.names,estimate])
                print.(term.names,",")
                println()
            else
                push!(grid,[tratii,effecti,mean(term.val),mean(term.val)*estimate])
                println(mean(term.val),",")
            end
        end
    end
    println("lsmeans:")
    for vari in grid
        out = 0.0
        for i in 1:length(vari[3])
            out += vari[4][i]
            for varj in deleteat!(copy(grid),i)
                out += mean(varj[4])
            end
            println(vari[1]," ",vari[2]," ",vari[3][i]," ",out)
        end
    end
end


export traceplot, PSRF

#adjust phenotype
end
