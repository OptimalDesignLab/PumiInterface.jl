# some temporary definitions
immutable MySBP{Tsbp} <: AbstractSBP{Tsbp}
  degree::Int
  numnodes::Int
  numfacenodes::Int
end

immutable MyFace{Tsbp} <: SummationByParts.AbstractFace{Tsbp}
  degree::Int
  numnodes::Int
  stencilsize::Int
end


