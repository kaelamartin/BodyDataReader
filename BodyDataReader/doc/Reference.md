This document is a reference for all functions exported by BodyDataReader.

# Table of Contents

- [`boddat`](#boddat)
- [Dictionary Handling](#dictionary-handling)
- [General Data Retrieval](#general-data-retrieval)
- [Small-Body Data Retrieval](#small-body-data-retrieval)
- [Ephemeris Retrieval](#ephemeris-retrieval)
- [Orientation Retrieval](#orientation-retrieval)

# `boddat`

```julia
argout = boddat(bvar::String, bi::Vector{Integer}, dict::Dict{String,Any})

boddatephem!(X::Vector, bi::Vector{Integer}, dict::Dict{String,Any})

orientationVec = boddatorient(parameter::String, body::Integer, dict::Dict{String,Any})
```

Three different `boddat` functions are exported. The two specialized versions (ephem and orient) are generally faster than `boddat` but still not suitable for performance-centric code.

# Dictionary Handling

`gendict!(dict::Dict{String,Any})`

Populates an empty dictionary with the required structures for use with BodyDataReader. Can also be used to set the options of a populated dictionary.

`tryloaddict!(dict::Dict{String,Any})`

Populates an empty dictionary with the required structures for use with BodyDataReader thn tries to load the previous dictionary from file.

`savedict(dict::Dict{String,Any})`

Saves the dictionary to file. Saves to the project working directory at the time the dictionary was populated.

Direct use of BodyDataReader functions will not work unless the dictionary provided to them has already been formatted for use.

# General Data Retrieval

`b::Integer = getspk(body::String, dict::Dict{String,Any})`

Attempts to parse string `body` to an SPK ID.

`n::String = getname(body::Integer, dict::Dict{String,Any})`

Returns the name associated with `body`.

`gm::Float = getgm(body::Integer, dict::Dict{String,Any})`

Returns the standard gravitational parameter of `body` in km^3/s^2.

`j2::Float = getj2(body::Integer, dict::Dict{String,Any})`

Returns the second dynamic form factor (J2) of `body`. J2 is only available for select bodies.

`rx:Float = getradius(body::Integer, dict::Dict{String,Any})`

Returns the radius associated with `body` in km. Typically the mean or equitorial radius.

`tri::Array{Float} = gettriaxial(body::Integer, dict::Dict{String,Any})`

Returns the triaxial radii of `body` in km.

`r::Float = getrotperiod(body::Integer, dict::Dict{String,Any})`

Returns the rotation period of `body` in solar days.

`cb::Integer = getcb(body::Integer)`

Returns the central body of `body`.

# Small-Body Data Retrieval

All functions that work for major bodies also work for small bodies, if the requested data is available. There are also additional functions specifically for small bodies, listed here.

`h::Float = getmagnitude(body::Integer, dict::Dict{String,Any})`

Returns the absolute magnitude of `body`.

`occ::Integer = getocc(body::Integer, dict::Dict{String,Any})`

Returns the orbital condition code of `body`.

`typ::String = gettype(body::Integer, dict::Dict{String,Any})`

Returns the spectral type (Tholen/SMASII taxonomic classification) of `body`.

`approach_Data::Array = getclose(body::Integer, dict::Dict{String,Any})`

Returns an array of data on the close approaches of `body`. Each row is one data entry in the following format: \[Date of closest approach (J2000)\] \[Other body\] \[Nominal distance (AU)\]

# Ephemeris Retrieval

Ephemeris retrieval commands require several additional parameters, most importantly the `EphemType` input. `EphemType` is an enumerable with the following values:

```julia
position = 1
velocity = 2
state = 3
acceleration = 4
```

Note that data retrieved using the `acceleration` ephemeris type includes position, velocity, acceleration, and jerk.

```
ephem!(body::Integer, time::Float, type::EphemType, dict::Dict{String,Any}, X::Array{Float})
ephem!(body::Integer, timeRange::Array{Float}, type::EphemType, dict::Dict{String,Any}, X::Array{Float})
```

Basic ephemeris retrieval commands for `body` relative to the Solar System Barycenter. Output is stored in `X`, which must have size 6 for position, velocity, and state output, and size 12 for acceleration output.
If a time range is used, `X` must have size 6 or 12 by n, where n is the number of points in the time range.

```
ephem!(bodies::Array{Integer}, time::Float, type::EphemType, dict::Dict{String,Any}, X::Array{Float})
ephem!(bodies::Array{Integer}, timeRange::Array{Float}, type::EphemType, dict::Dict{String,Any}, X::Array{Float})
```

Relative ephemeris retrieval. `bodies` must be an array of size two. Ephemeris `body[1]` is retrieved relative to `body[2]`. Other arguments are identical to basic ephemeris retrieval.

```
X::Array{Float} = ephem(body::Integer, time::Float, type::EphemType, dict::Dict{String, Any})
X::Array{Float} = ephem(body::Integer, timeRange::Array{Float}, type::EphemType, dict::Dict{String, Any})
```

Basic ephemeris retrieval relative to the Solar System Barycenter. Output is returned instead of stored.

Two additional ephemeris-related utilities are provided as well.

`reference::String = getreference(body::Integer, dict::Dict{String, Any})`

Returns the reference file for the ephemeris file of `body`.

`date::Float = getdate(body::Integer, dict::Dict{String,Any})`

Returns the date of ephemeris retrieval for `body`.

# Orientation Retrieval

```
pole:Array{Float} = pole(body::Integer, time::Float, dict::Dict{String,Any})
pole:Array{Float} = pole(body::Integer, timeRange::Array{Float}, dict::Dict{String,Any})
```

Returns the pole unit vector of `body` at a specific time or over a time range.

```
pm:Array{Float} = pm(body::Integer, time::Float, dict::Dict{String,Any})
pm:Array{Float} = pm(body::Integer, timeRange::Array{Float}, dict::Dict{String,Any})
```

Returns the prime meridian unit vector of `body` at a specific time or over a time range.

```
eqx:Array{Float} = eqx(body::Integer, time::Float, dict::Dict{String,Any})
eqx:Array{Float} = eqx(body::Integer, timeRange::Array{Float}, dict::Dict{String,Any})
```

Returns the equinox unit vector of `body` at a specific time or over a time range.

Additionally, three direction cosine maxtrix retrieval commands are available, one for each orientation unit vector.

`poledcm::Array{Float64} = poledcm(body::Integer, time::Float, dict::Dict{String,Any})`

Returns the pole DCM vector of `body`, which contains the following data: DCM = \[Node;Q;Pole\] where Node = \[0;0;1\] x Pole and Q = Pole x Node.

`pmdcm::Array{Float64} = pmdcm(body::Integer, time::Float, dict::Dict{String,Any})`

Returns the prime meridian DCM vector of `body`, which contains the following data: DCM = \[Node;Q;Pole\] where Node = prime meridian and Q = Pole x Node.

`eqxdcm::Array{Float64} = eqxdcm(body::Integer, time::Float, dict::Dict{String,Any})`

Returns the equinox DCM vector of `body`, which contains the following data: DCM = \[Node;Q;Pole\] where Node = equninox and Q = Pole x Node.

Additionally, a similar set of "quick" commands are available. These commands approximate orientation parameters based off of known values at J2000 epoch.
They can be extremely inaccurate, particularly for satellite and small bodies.

`pole:Array{Float} = qpole(body::Integer, dict::Dict{String,Any})`

Quickly approximates the pole unit vector of `body`.

`pm:Array{Float} = qpm(body::Integer, dict::Dict{String,Any})`

Quickly approximates the prime meridian unit vector of `body`.

`eqx:Array{Float} = qeqx(body::Integer, dict::Dict{String,Any})`

Quickly approximates the equinox unit vector of `body`.

`poledcm::Array{Float64} = qpoledcm(body::Integer, time::Float, dict::Dict{String,Any})`

Quickly approximates the pole DCM vector of `body`.

`pmdcm::Array{Float64} = qpmdcm(body::Integer, time::Float, dict::Dict{String,Any})`

Quickly approximates the prime meridian DCM vector of `body`.

`eqxdcm::Array{Float64} = qeqxdcm(body::Integer, time::Float, dict::Dict{String,Any})`

Quickly approximates the equinox DCM vector of `body`.