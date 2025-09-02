function init_tram_tracker(A::SparseMatrixCSC; N = 4, D = 20, L = 4)
    tracker = TramTrackerAlgo(N, D, L, A)
end

function find_N_nearest_edges(
        tracker::TramTrackerAlgo, measurement::T) where {T <: Meshes.Point}
    edge_ids = findall(x -> x > 0, tracker.A)
    edge2dist = Dict(edge_ids .=> 0.0)
    for edge_id in edge_ids
        a, b = tracker.int2point[edge_id]
        edge2dist[edge_id] = dist(measurement, a, b)
    end
    idx = partialsortperm(collect(values(edge2dist)), 1:(tracker.N))
    return collect(keys(edge2dist))[idx]
end
