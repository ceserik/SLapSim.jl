using Proj

# Proj.jl implements the CoordinateTransformations.jl API.
# A Proj.Transformation needs the source and target coordinate reference systems (CRS),
# or a single pipeline.
trans = Proj.Transformation("EPSG:4326", "+proj=merc +lat_ts=49")
# The CRS can be a string or the CRS type, which also interfaces with GeoFormatTypes.jl.

# Once created, you can call this object to transform points.
# The result will be a tuple of Float64s, of length 2, 3 or 4 depending on the input length.
# The 3rd coordinate is elevation (default 0), and the 4th is time (default Inf).
# Here the (latitude, longitude) of Copenhagen is entered
trans(55, 12)

