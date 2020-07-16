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


export QC
export checkFile
