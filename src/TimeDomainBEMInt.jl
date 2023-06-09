module TimeDomainBEMInt

using CompScienceMeshes

using WiltonInts84

using LinearAlgebra

using BlockArrays

using Polynomials

include("integrals.jl")
include("inttriangleRminus4.jl")
include("contour2.jl")

export inttriangletriangle

end
