mutable struct FixedMatrix
    C::Array{Float64,2}
    variables::Array{Any,1}
end


function make_fixed(file;ID_order=false) #ID_order deprecated, not useful when multi-obs on same ind

    myfile = open(file)
    #set number of columns and column types
    row1   = split(readline(myfile))

    ncol   = length(row1)
    etv    = Array{DataType}(ncol)
    etv[1] = String  #type for 1st column: ID

    for coli in 2:length(row1) #type for all variables(predictors)
       if row1[coli][end]=='#'
           etv[coli]=String
       else
           etv[coli]=Float64
       end
    end
    close(myfile)

    #read genotypes
    df = readtable(file, eltypes=etv, separator = ' ',header=true)

    variables=Array{Any}(length(etv)-1)
    if etv[2]==Float64
        C = convert(Array,df[:,2])
        variables[1] = "covariate"
    else
        C,variables[1]=mkmat_incidence_factor(df[:,2])
    end

    for i=3:length(etv)
      if etv[i]==Float64
          C = [C convert(Array,df[:,i])]
          variables[i-1] = "covariate"
      else
          myfactor=mkmat_incidence_factor(df[:,i])
          C,variables[i-1]=[C myfactor[1]],myfactor[2]
      end
    end

    #useful trick below anyway
    #if ID_order!= false
    #    ID_order=convert(DataFrame,reshape(ID_order,length(ID_order),1))#covert to dataframe for 2-diim
    #end
    #C=convert(DataFrame,[Array(df[:,1]) C]);
    #reorder C to make it has the row order as y: phenotypes
    #Cnew=Array{Any}(size(ID_order,1),size(C,2))
    #j=1
    #for i in ID_order[:x1]
    #    index = find(C[:x1].==i)
    #    Cnew[j,:]= convert(Array,C[index,:])
    #    j = j+1
    #end
    #
    #C=Cnew[:,2:end]

    return FixedMatrix(C,variables)
end
#export FixedMatrix
