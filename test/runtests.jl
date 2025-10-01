using Test
using FrequencyMaxwell
using LinearAlgebra: norm

include("solvers/boundary_conditions.jl")
include("solvers/convergent_born.jl")
include("geometry/phantoms.jl")
include("sources/plane_wave.jl")
include("fields/electromagnetic_field.jl")
