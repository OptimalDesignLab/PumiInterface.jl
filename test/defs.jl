# some temporary definitions
struct MySBP{Tsbp} <: AbstractSBP{Tsbp}
  degree::Int
  numnodes::Int
  numfacenodes::Int
end

struct MyFace{Tsbp} <: SummationByParts.AbstractFace{Tsbp}
  degree::Int
  numnodes::Int
  stencilsize::Int
end


