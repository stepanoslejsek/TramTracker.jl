using GeoIO
using GeoStats

function load_tram_map(path::String)
    """
        Load a .geojson file specified by path and returns an adjacency matrix of the tram track
    """
    gt = GeoIO.load(path)
    A = create_graph(gt)
    return A
end

function create_graph(gt::GeoTable)
    processing_pipeline = DropMissing("railway") → Filter(x -> x.railway == "tram")
end
