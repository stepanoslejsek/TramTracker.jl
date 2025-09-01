using TramTracker
using Test

@testset "TramTracker.jl" begin
    # Write your tests here.
    @test tram_map = load_tram_map("./data/test_data.geojson")
end
