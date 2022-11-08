using NICE
using Test

@testset "NICE.jl" begin
    # Write your tests here.
    #
    # I just checked that it runs...
    a = ReactionSystem(zeros(3,3), zeros(3))
    @show a
end
