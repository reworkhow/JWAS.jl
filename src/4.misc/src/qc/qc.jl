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

"""
    QC(infile,outfile;separator=' ',header=true,missing=false,MAF=0.1)
* Quality control for input file **infile**, then write out to output file: **outfile**.
  * Delete loci with minor allele frequency < **MAF**.
  * **missing** genotypes are replaced by column means.
* File format (header=true,separator=',',missing=9):

```
Animal,marker1,marker2,marker3,marker4,marker5
S1,1,0,1,1,1
D1,2,0,9,2,1
O1,1,2,0,1,0
O3,0,0,2,1,1
```
"""
function QC(infile,outfile;separator=' ',header=true,missing=false,MAF=0.1)

    myfile = open(infile)
    id4row = true

    #set types for each column
    ncol= length(split(readline(myfile)))
    etv = Array(DataType,ncol)
    fill!(etv,Float64)
    if id4row == true
      etv[1]=String
    end
    close(myfile)

    #read genotypes
    df = readtable(infile, eltypes=etv, separator = separator,header=header)
    #quality control
    missing2mean!(df,id4row=id4row,missing=missing)
    df=deleteMAF!(df,id4row=id4row,MAF=MAF)

    #write out the new file
    writetable(outfile, df, separator=separator)
end


function missing2mean!(X::DataFrames.DataFrame;id4row=true,missing=9.0)
    if missing == false
        return
    end
    nrow,ncol = size(X)
    start = 1
    if id4row==true
      start+=1
    end
    for i=start:ncol
        index=find(x->x==missing,X[:,i])
        cols = collect(1:nrow)
        deleteat!(cols,index)
        #X[index,i]=round(Int,mean(X[cols,i]))
        X[index,i]=mean(X[cols,i])
    end
end

function deleteMAF!(df;id4row=false,MAF=0.1)
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
