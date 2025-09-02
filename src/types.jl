using GeoStats

struct TramTrackPoint{T <: Meshes.Point}
    id::Int
    location::T
end

struct TramTrackerAlgo
    N::Int # Number of nearest edges to the current position of the tram
    D::Float64 # Distance threshold for discarding path
    L::Int # Number of vertices in the directed path of the graph
    A::SparseMatrixCSC
    point2int::Dict{T <: Meshes.Point, Int}
    int2point::Dict{Int, T <: Meshes.Point}
end
