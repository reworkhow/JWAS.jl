function checkFile(file)
  #count row numbers
  f=open(file)
  nrow = countlines(f)
  close(f)

  #check column numbers
  f=open(file)
  ncol=length(split(readline(f)))
  while !eof(f)
    if length(split(readline(f))) != ncol
      error("Number of columns are not equal!!!")
    end
  end

  #SNPs = split(readline(f))
  #for i in SNPs
  #if i in SNPs
  #println("Missing or imputed values are found")
  #count+=1
  #end

  println("Good file with ",nrow, " rows and ",ncol," columns!!!")
  close(f)
end

function QC(infile,outfile;id4col=true,missing=9.0,MAF=0.1)

    myfile = open(infile)

    #set types for each column
    ncol= length(split(readline(myfile)))
    etv = Array(DataType,ncol)
    fill!(etv,Float64)
    if id4row == true
      etv[1]=UTF8String
    end
    close(myfile)

    #read genotypes
    df = readtable(file, eltypes=etv, separator = ' ',header=id4col)
    #quality control
    missing2mean(df,id4row=id4row,missing=missing)
    deleteMAF!(df,id4row=id4row,MAF=MAF)

    #write out the new file
    writetable(outfile, df, separator=' ')
end


function missing2mean!(X::DataFrames.DataFrame;id4row=false,missing=9)
    nrow,ncol = size(X)
    start = 1
    if id4row==true
      start+=1
    end
    for i=start:ncol
        index=find(x->x==missing,X[:,i])
        cols = collect(1:nrow)
        deleteat!(cols,index)
        X[index,i]=int(mean(X[cols,i]))

        if(i%3000==0)
            println("This is line ",i)
        end
    end
end

function deleteMAF!(X;id4row=false,MAF=0.1)
    start = 1
    if id4row==true
      start+=1
    end

    sel = MAF.<[mean(df[:,i])/2 for i in start:size(df,2)].< 1-MAF

    df = df[:,sel]
    return df
end

export QC
export checkFile
