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
#M = Array{Float64}(undef,nrows,ncols)
#read!(s,M)
#close(s)

#only one variable
#M = Array{Float64}(undef,nrows,ncols)
#read!("1.bin",M)

#latex table from matrix
# fileName::AbstractString = "",
# colLabels = "",
# rowLabels = ""

# ┌ Warning: Deprecated syntax `parametric method syntax latexTable{T}(A::Array{T, 2},; fileName::AbstractString = "", colLabels = "", rowLabels = "")` around /home/jovyan/qtlcheng/Github/JWAS.jl/src/4.misc/src/files.jl:27.
# │ Use `latexTable(A::Array{T, 2},; fileName::AbstractString = "", colLabels = "", rowLabels = "") where T` instead.
# └ @ nothing /home/jovyan/qtlcheng/Github/JWAS.jl/src/4.misc/src/files.jl:27
# function latexTable{T}(A::Array{T,2};
#     fileName::AbstractString = "",
#     colLabels = "",
#     rowLabels = ""
#     )
#     if fileName == ""
#         outFile = STDOUT
#     else
#         outFile = open(fileName, "w")
#     end
#     rows, cols = size(A)
#     println(outFile,"\\begin{center}")
#     print(outFile,"\\begin{tabular}{",)
#     if rowLabels!=""
#         print(outFile,"l")
#     end
#     for j=1:cols
#         print(outFile,"r")
#     end
#     println(outFile,"}\\hline")
#     if colLabels!=""
#         nCol = length(colLabels)
#         for j = 1:(nCol-1)
#             print(outFile,colLabels[j]," & ")
#         end
#         print(outFile,colLabels[nCol]," \\\\ \\hline  \n")
#     end
#     for i =1:rows
#         if rowLabels!=""
#             print(outFile,rowLabels[i]," & ")
#         end
#         for j = 1:(cols-1)
#             print(outFile,A[i,j])
#             print(outFile," & ")
#         end
#         print(outFile,A[i,cols]," \\\\ \n")
#     end
#     println(outFile,"\\hline")
#     println(outFile,"\\end{tabular}")
#     println(outFile,"\\end{center}")
#     if fileName != ""
#         close(outFile)
#     end
# end
