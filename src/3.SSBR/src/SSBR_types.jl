mutable struct Numbers
  ped::Int64      #individulals in pedigree
  pedn::Int64     #non-genotyped individuals in pedigree
  pedg::Int64     #genotyped individuals in pedigree
  y::Int64        #individuals with phenotypes
  yn::Int64       #non-genotyped individuals with phenotypes
  yg::Int64       #genotyped individuals with phenotypes
  markers::Int64  #number of markers
end

mutable struct ZMats
  full::SparseMatrixCSC{Float64,Int64}   #Z   # Z= | Zn 0 |=|Z_n Z_g|
  n::SparseMatrixCSC{Float64,Int64}      #Zn  #    | 0  Zg|
  g::SparseMatrixCSC{Float64,Int64}      #Zg  #
  _n::SparseMatrixCSC{Float64,Int64}     #Z_n # Z_n=|Zn|
  _g::SparseMatrixCSC{Float64,Int64}     #Z_g #     |0 |
end

mutable struct AiMats
  full::SparseMatrixCSC{Float64,Int64}    #Ai
  nn::SparseMatrixCSC{Float64,Int64} #Ai_nn
  ng::SparseMatrixCSC{Float64,Int64} #Ai_ng
end

mutable struct YVecs
  full #::Array{Float64,1}   #y
  n    #::Array{Float64,1}      #yn
  g    #::Array{Float64,1}      #yg  #order of ids is same to order of y
  ids  #::Array{String,1} #order of ids is nongeno then geno
end

mutable struct MMats
  full::Array{Float64,2}
  n::Array{Float64,2} #Mn
  g::Array{Float64,2} #Mg
end

mutable struct JVecs
  full::Array{Float64,2}
  n::Array{Float64,2} #Jn
  g::Array{Float64,2} #Jg
end

mutable struct XMats #sparse may be better
  full::Array{Float64,2}
  n::Array{Float64,2} #Xn
  g::Array{Float64,2} #Xg
end

mutable struct WMats
  full::Array{Float64,2}
  n::Array{Float64,2} #Wn
  g::Array{Float64,2} #Wg
end

mutable struct HybridMatrices
  Z::ZMats
  Ai::AiMats
  y::YVecs
  J::JVecs
  X::XMats
  W::WMats
  M::MMats
  num::Numbers
end

mutable struct Results
  mats::HybridMatrices
  posMeans::misc.Output
  ebv
end
