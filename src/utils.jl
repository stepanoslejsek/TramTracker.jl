using GeoIO
using GeoStats

function dist(p₁::T, p₂::T) where {T <: Meshes.Point}
    return haversine(
        (p₁.coords.lon.val, p₁.coords.lat.val), (p₂.coords.lon.val, p₂.coords.lat.val))
end

function dist(p::T, a::T, b::T ) where {T <: Meshes.Point}
    # IMPORTANT: Not using Haversine formula for computing the distance
    # Using straightforward formula as the error is small
    ab = b - a
    ap = p - a
    t = dot(ap, ab) / dot(ab, ab)
    if t < 0
        return norm(p - a).val
    elseif t > 1
        return norm(p - b).val
    else
        proj = a + t .* ab
        return norm(p - proj).val
    end
end

function load_tram_map(path::String)
    """
        Load a .geojson file specified by path and returns an adjacency matrix of the tram track
    """
    gt = GeoIO.load(path)
    A = create_graph(gt)
    return A
end

function create_graph(gt::GeoTable)
    # Apply pipeline and sort
    processing_pipeline = DropMissing("railway") → Filter(x -> x.railway == "tram")
    gt = sort_track_segments(gt |> processing_pipeline)

    track_nodes = unique(pointify(gt.geometry))
    N = length(track_nodes)

    # Map points to integer indices
    points_mapping = Dict(track_nodes .=> 1:N)
    point2int = Dict(track_nodes .=> 1:N)
    int2point = Dict(1:N .=> track_nodes)

    # Arrays for sparse triplets
    I = Int[]
    J = Int[]
    V = Float64[]

    # Build edges
    for rope in gt.geometry
        points = pointify(rope)
        for i in 1:(length(points) - 1)
            push!(I, points_mapping[points[i]])
            push!(J, points_mapping[points[i+1]])
            push!(V, dist(points[i], points[i+1]))
        end
    end
    # Construct sparse adjacency matrix
    return sparse(I, J, V, N, N)
end

function sort_track_segments(gt::GeoTable)
    segments = pointify.(gt.geometry)
    sorted_segments_ids = [1]
    unsorted_segments_ids = collect(2:length(segments))

    for i in 1:(length(segments) - 1)
        start = sorted_segments_ids[1]
        stop = sorted_segments_ids[end]

        best_dist = Inf
        best_idx = 0
        state = 0

        for j in 1:length(unsorted_segments_ids)
            curr = unsorted_segments_ids[j]

            # Compute all possible combination of distances between two segments
            d_preceding = dist(segments[curr][end], segments[start][1])
            d_succeeding = dist(segments[curr][1], segments[stop][end])
            d_preceding_rev = dist(segments[curr][1], segments[start][1])
            d_succeeding_rev = dist(segments[curr][end], segments[stop][end])

            distances = [d_preceding, d_succeeding, d_preceding_rev, d_succeeding_rev]

            if minimum(distances) < best_dist
                state = argmin(distances)
                best_dist = minimum(distances)
                best_idx = j
            end
        end
        if state == 1
            pushfirst!(sorted_segments_ids, unsorted_segments_ids[best_idx])
            deleteat!(unsorted_segments_ids, best_idx)
        elseif state == 2
            push!(sorted_segments_ids, unsorted_segments_ids[best_idx])
            deleteat!(unsorted_segments_ids, best_idx)
        elseif state == 3
            curr = unsorted_segments_ids[best_idx]
            segments[curr] = reverse(segments[curr])
            pushfirst!(sorted_segments_ids, unsorted_segments_ids[best_idx])
            deleteat!(unsorted_segments_ids, best_idx)
        elseif state == 4
            curr = unsorted_segments_ids[best_idx]
            segments[curr] = reverse(segments[curr])
            push!(sorted_segments_ids, unsorted_segments_ids[best_idx])
            deleteat!(unsorted_segments_ids, best_idx)
        end
    end
    return gt[sorted_segments_ids, :]
end

Base.getindex(M::SparseMatrixCSC, r::T, c::T) where {T <: Meshes.Point} = M[findall(x->x==), findall()]
Base.getindex(D::Dict{Int, T}, I::CartesianIndex) where {T <: Meshes.Point} = D[I.I[1]], D[I.I[2]]
