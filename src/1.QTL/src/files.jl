#write out to text: writedlm is Okay, (18904,50564) matrix takes 305 seconds
#write out to binary

function write2binary(x,file)
  myfile=open(file,"w")
  write(myfile,x)
  close(myfile)
end

#s = open("123.bin")
#nrows = read(s,Int64)
#ncols = read(s,Int64)
#M = Array(Int64,nrows,ncols)
#read!(s,M)
#close(s)
