# BodyDataReader.jl

BodyDataReader is an ephemeris and physical constants reader that makes it easy to retrieve data on any body in the solar system.

Data is retrieved from [JPL Horizons](https://ssd.jpl.nasa.gov/?horizons), the [JPL Small-Body Database](https://ssd.jpl.nasa.gov/sbdb.cgi), and the [NAIF SPICE system](https://naif.jpl.nasa.gov/naif/index.html).
Retrieved data is saved locally for fast access later.

## Installation

BodyDataReader can be installed through Julia's package manager:
```julia
Pkg.add("BodyDataReader")
```

## Documentation

Take a look at the [manual](doc/Manual.md) for an explanation of how to use BodyDataReader.

A listing of every available function can be found the [reference](doc/Reference.md).

## Quick Start

BodyDataReader has many functions, each specialized to retrieve a different kind of data. It also provides a single function call to access every other function, `boddat`:

```julia
julia> using BodyDataReader

julia> dict = Dict{String,Any}()

julia> boddat("R", [399], dict, [0.])
3-element Array{Float64,1}:
    -2.756663225908748e7
     1.44279062385979e8
 30226.39667461067
```

This example retrieves the position of Earth in kilometers relative to the Solar System Barycenter at J2000 epoch.

`boddat` is capable of retrieving many types of data, far too many to list in this quick start. A full listing of available data can be found using Julia's built-in help command.

Direct access to data-retrieval functions is also available and is recommended for any task requiring high performance. See the docs for more information.

## Copyright

Copyright (c) 2020 California Institute of Technology (“Caltech”) and Embry Riddle Aeronautical University. U.S. Government sponsorship acknowledged. 
All rights reserved. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 

Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
