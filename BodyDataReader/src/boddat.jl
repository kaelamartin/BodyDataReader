# Copyright (c) 2020 California Institute of Technology (“Caltech”) and Embry Riddle Aeronautical University.
# U.S. Government sponsorship acknowledged.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#  • Redistributions of source code must retain the above copyright notice, this 
#    list of conditions and the following disclaimer.
#   • Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
#  • Neither the name of Caltech nor its operating division, the Jet Propulsion 
#    Laboratory, nor the names of its contributors may be used to endorse or promote 
#    products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
# SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.






"""
    argout = boddat(bvar, bi, dict, t, opti)

Main function that can call all other functions. Takes the following inputs:
  * `bvar` is the string name of the desired data (accepted strings below).
  * `bi` is a float vector of either one or two SPK numbers. Two values are only valid for ephemeris commands, and will result in the first body relative to the second.
  * `dict` is dictionary in which data is stored between calls. It must have type Dict{String,Any} and can be empty.
  * `t` is an optional float vector which is the times after J2000 in julian days that the requested data will be retrieved at. If not provided, the current time is used.
  * `opti` is an optional boolean vector of that specifies the following three options:
    * `save` - whether data should be saved. (default true)
    * `ssd` - whether data should always be retrieved from the JPL SSD. This option may break some data requests. (default false)
    * `keep_spline` - whether spline data should be kept for later calls. (default true)

Body data is retrieved from JPL Horizons and the JPL Solar System Dynamics
website (ssd.jpl.nasa.gov).

General bvar values:
  * "Number" - SPK number
  * "Name" - name
  * "CB" - central body
  * "GM" - gravitational parameter (km^3/s^2)
  * "radius" - mean radius (km)
  * "J2" - gravitational harmonic due to oblateness
  * "rotation_period" - sidereal rotation period (solar days)

Ephemeris bvar values:
  * "R" - position vector relative to SSB or second body (km)
  * "V" - velocity vector relative to SSB or second body (km/s)
  * "X" - state vector relative to SSB or second body (km,km/s)
  * "A" - acceleration vector relative to SSB or second body (km,km/s,km/s^2,km/s^3)
  * "RTN" - rotating frame direction cosine matrix unit[R ; RxVxR ; RxV]
  * "reference" - ephemeris reference file
  * "date" - creation date of ephem file (days since J2000)

Small body specific bvar values:
  * "mag" - absolute magnitude
  * "occ" - orbital condition code
  * "type" - spectral type, Tholen/SMASII taxonomic classification
  * "close_approach" - table of close approach times, bodies, and distances

Orientation bvar values:
  * "Pole" - spin axis unit vector
  * "PM" - prime meridian unit vector
  * "Eqx" - equinox unit vector (Pole x H)
  * "Poledcm" - pole direction cosine matrix unit[Node;Q;Pole] where Node = [0;0;1] x Pole, Q = Pole x Node
  * "PMdcm" - prime meridian direction cosine matrix unit[Node;Q;Pole] where Node = prime meridian, Q = Pole x Node
  * "Eqxdcm" - equinox direction cosine matrix unit[Node;Q;Pole] where Node = equinox, Q = Pole x Node

For orientation `bvar` values the letter "q" can be appended on the front of each
command for "quick" orientation versions, which approximate values based off of
known values at J2000. This can be highly inaccurate for satellite bodies.

All vectors are output in ecliptic and mean equinox of J2000 reference frame.

Note that data retrieval may take a up to a minute on the first run, and is much
faster once data has been cached. Cache files are saved at the current project
working directory, which can be determined with the command pwd().

"""
function boddat(bvar::AbstractString,bi::AbstractArray{Int64, 1},
    dict::Dict{String,Any}=Dict{String,Any}(),
    t::AbstractArray{P}=Float64[],
    opti::AbstractArray{Bool}=Bool[]) where P

    #SOURCES
    #small body data https://ssd.jpl.nasa.gov/sbdb.cgi
    #Horizons ephemeris and orientation https://ssd.jpl.nasa.gov/horizons_batch.cgi
    #list of major bodies
    #       https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=MB
    #Ephemeris file names and time spans https://ssd.jpl.nasa.gov/eph_spans.cgi
    #Inner planet, sun & barycenter GM, planet radius J2
    #       ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/
    #Outer planet & satellite GM, planet radius J2
    #       ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/
    #Satellite GM and mean radius https://ssd.jpl.nasa.gov/?sat_phys_par
    #analytic orientation and radii
    #       ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/

    varargout = Array{Any}(undef, 1)

    if isempty(dict)
        tryloaddict!(dict)
    end

    if length(t) == 1 # Reformat time as scalar
        t = t[1]
    end

    #if there's no time, the time is now
    if isempty(t)
        t = Dates.datetime2julian(now())-2451545.
    end

    if size(t,2)>size(t,1)
        t = transpose(t)
    end
    (eltype(t)==Int64) && (t = t*1.)

    #default options is [true false true], replace with any input and save in params
    opts=[true false true]
    opts[1:length(opti)]=opti
    param(dict,"save",opts[1])
    param(dict,"ssd",opts[2])
    param(dict,"keeps",opts[3])

    if eltype(bi) != String && eltype(bi) != Any
        be=copy(bi)
    else
        be = getspk(bi,dict)
    end
    b = zeros(Int64,size(be))

    for ii in 1:length(be)
        b[ii]=convert(Int64,real(be[ii]))
    end

    #find proper parameter and call function

    bvi = Int64[]

    if length(bvar)==1
        inBVA = findall(x->(occursin(bvar,x)),param(dict,"bva")[1:4])
    else
        inBVA = findall(x->(occursin(bvar,x)),param(dict,"bva"))
    end
        (!isempty(inBVA)) && (inBVA = inBVA[1])
    #find out which member of bva matches input
    if isempty(inBVA)
        bvar = uppercase(bvar)
        inBVA = findall(x->(occursin(bvar, x)),param(dict,"bva"))
        (!isempty(inBVA)) && (inBVA = inBVA[1])
        if isempty(inBVA)
            bvar = lowercase(bvar)
            inBVA = findall(x->(occursin(bvar,x)),param(dict,"bva"))
            (!isempty(inBVA)) && (inBVA = inBVA[1])
            if isempty(inBVA)
                bvar = uppercasefirst(bvar)
                inBVA = findall(x->(occursin(bvar, x)),param(dict,"bva"))
                (!isempty(inBVA)) && (inBVA = inBVA[1])
                if isempty(inBVA)
                    bvar = lowercasefirst(bvar)
                    bvar = replace(bvar,"pm" => "PM")
                    bvar = replace(bvar,"pole" => "Pole")
                    bvar = replace(bvar,"eq" => "Eq")
                    inBVA = findall(x->(occursin(bvar,x)),param(dict,"bva"))
                    (!isempty(inBVA)) && (inBVA = inBVA[1])
                end
            end
        end
    end
    (isempty(inBVA)) && (inBVA = 0)
    if inBVA > 4
        T1=uniquestr(param(dict,"bva"),bvar)
        bvi = T1[1]
    elseif inBVA <= 4
        bvi = inBVA
    end
    if iszero(bvi)
        @warn "no match found for $bvar"
        return NaN*ones(size(t))
    end
    #match bvi

    if length(b) == 1
        b = b[1]
    end

    if bvi > 5 && length(b) > 1
        @warn "Vector body input is only valid for ephemeris commands, returning data for first body"
    end

    if bvi < 5#ephemeris, {R V X A}
        n = 6
        (bvi == 4) && (n = 12)

        if (!param(dict,"ssd"))||(!isdefined(:Xeph))
            if typeof(t) == Float64
                Xeph = Array{eltype(t)}(undef, n)
            else
                Xeph = Array{eltype(t)}(undef, n, length(t))
            end
            ephem!(b,t,EphemType(bvi),dict,Xeph)
        end

        a = Xeph
        if bvi == 1
            a = a[1:3]
        elseif bvi == 2
            a = a[4:6]
        end
    elseif bvi == 5
        if typeof(t) == Float64
            Xeph = Array{eltype(t)}(undef, 6)
        else
            Xeph = Array{eltype(t)}(undef, 6, length(t))
        end
        ephem!(b,t,state,dict,Xeph)
        R=Xeph[1:3,:]
        V=Xeph[4:6,:]
        uH = zero(R)
        uR = zero(R)
        uHxR = zero(R)
        for ii in 1:size(R,2)
            H = cross(R[:,ii],V[:,ii])
            HxR = cross(H,R[:,ii])
            Hm = 0.; HxRm = 0.; Rm = 0.
            for kk in 1:3
                Hm = H[kk]^2 + Hm
                HxRm = HxR[kk]^2 + HxRm
                Rm = R[kk,ii]^2 + Rm
            end
                Rm = sqrt(Rm)
                HxRm = sqrt(HxRm)
                Hm = sqrt(Hm)
            for kk in 1:3
                uH[kk,ii] = H[kk]/Hm
                uHxR[kk,ii] = HxR[kk]/HxRm
                uR[kk,ii] = R[kk,ii]/Rm
            end
        end
        a = [uR; uHxR; uH]
    elseif bvi == 6
        a = getnumout(b,dict) #number
    elseif bvi == 7#central body, matches bi for name vs number output
		a = getcb(b)
    elseif bvi == 8
        #a = Array{String}(length(b))
        #for ii in 1:length(b)
        #    T2 = getname(b[ii],dict)
        #    a[ii]=T2[1] #name, works for NonBodyControlPoints
        #end
        a = zeros(3,length(b))
        for ii in 1:length(b)
            a[:,ii] = gettriaxial(b[1],dict)
        end
        #9--19 calls a function called getfn, which uses matching bva parameter to
        #      use the corresponding function  elseif bvi<20
    elseif bvi == 9
        for ii in 1:length(b)
            a = getname(bi[ii], dict)
        end
    elseif bvi < 20
        fn=param(dict,"bva")
        fn=fn[bvi]
        a = getfn(fn,b,dict)
    else #20--31 is orientation
        fn = param(dict,"bva")
        fn = fn[bvi]
        a = boddatorient(fn, b[1], dict, t, false)
    end#if bvi<5
    varargout = a #write output to varargout

    #save data to bodydata
    if param(dict,"save")
        savedict(dict)
    end

    return varargout
end

"""
	boddatephem!(X,bi,dict,t)

Quick ephemeris computation using input vector X for position, velocity, or
acceleration for a body number and filled dictionary.
"""
function boddatephem!(X::AbstractArray{P},bi::AbstractArray{Q},
                dict::Dict{String,Any},t::AbstractArray{P}=Float64[]) where {P, Q}
    # Faster boddat for ephemeris calls
    if typeof(t)==StepRangeLen{Float64,Base.TwicePrecision{Float64},
        Base.TwicePrecision{Float64}} || typeof(t)==StepRange{Int64,Int64}
        typ = true
        t = [t[1]; t[end]; t[2]-t[1]]
    else
        typ = false
    end

    if size(t,2)>size(t,1)
        t = transpose(t)
    end

    if iszero(imag(bi))
        b = Array{Int64}(undef, size(bi))
    else
        b = Array{Union{Complex{Int64},Int64}}(undef, size(bi))
    end

    for ii in 1:length(bi)
        if iszero(imag(bi[ii]))
            b[ii] = convert(Int64,real(bi[ii]))
        else
            b[ii] = convert(Complex{Int64},bi[ii])
        end
    end

    if size(X,1) == 12
        bvi = 4
    else
        bvi = 3
    end
    ephem!(b,t,typ,bvi,dict,X)
end

function boddatephem!(X::AbstractArray{P},bi::AbstractArray{Int64},
                dict::Dict{String,Any},t::AbstractArray{P}=Float64[]) where P
    # Faster boddat for ephemeris calls
    if typeof(t)==StepRangeLen{Float64,Base.TwicePrecision{Float64},
        Base.TwicePrecision{Float64}} || typeof(t)==StepRange{Int64,Int64}
        typ = true
        t = [t[1]; t[end]; t[2]-t[1]]
    else
        typ = false
    end

    if size(t,2)>size(t,1)
        t = transpose(t)
    end

    if size(X,1) == 12
        bvi = 4
    else
        bvi = 3
    end
    ephem!(bi,t,typ,bvi,dict,X)
end

"""
	orientationVec = boddatorient(parameter, body, dict, time, opts)

Retrieves the requested orientation parameter.
Recognized parameters are "Pole", "pm", "eqx", "Poledcm", "pmdcm", and "eqxdcm",
as well as their equivalent quick versions.
"""
function boddatorient(parameter::String, body::Int64, dict::Dict{String,Any},
     time::Float64 = -Inf64, save::Bool = true)::Array{Float64, 1}

     if isempty(dict)
         tryloaddict!(dict)
     end

     # Sort of a hack to determine if the user passed us a time, since the base
     # Float64 type can't be empty
     if time == -Inf64
         time = Dates.datetime2julian(now())-2451545.
     end

     parameter = lowercase(parameter)
     if(parameter[1] == 'q') # call quick versions
         func = parameter[2:end]
         if(func == "pole")
             result = qpole(body, dict)
         elseif(func == "pm")
             result = qpm(body, time, dict)
         elseif(func == "eqx")
	     @warn "May not be accurate for small bodies. Recommended to use eqx function for small bodies"		
             result = qeqx(body, dict)
         elseif(func == "poledcm")
             result = qpoledcm(body, dict)
         elseif(func == "pmdcm")
             result = qpmdcm(body, time, dict)
         elseif(func == "eqxdcm")
             result = qeqxdcm(body, dict)
         else
             @warn "Parameter not recognized"
             result = NaN*ones(3)
         end
     else # call normal versions
         func = parameter
         if(func == "pole")
             result = pole(body, time, dict)
         elseif(func == "pm")
             result = pm(body, time, dict)
         elseif(func == "eqx")
             result = eqx(body, time, dict)
         elseif(func == "poledcm")
             result = poledcm(body, time, dict)
         elseif(func == "pmdcm")
             result = pmdcm(body, time, dict)
         elseif(func == "eqxdcm")
             result = eqxdcm(body, time, dict)
         else
             @warn "Parameter not recognized"
             result = NaN*ones(3)
         end
     end

     if(save)
         savedict(dict)
     end

     return result
end

"""
	orientationArr = boddatorient(parameter, body, dict, time, opts)

Retrieves the requested orientation parameter over a time range.
Time range MUST be supplied for this function.
Recognized parameters are "Pole", "pm", "eqx", "Poledcm", "pmdcm", and "eqxdcm".
Quick versions are not supported by this function.
"""
function boddatorient(parameter::String, body::Int64, dict::Dict{String,Any},
     timeRange::AbstractArray{Float64}, save::Bool = true)::Array{Float64, 2}

     if isempty(dict)
         tryloaddict!(dict)
     end

     parameter = lowercase(parameter)
     if(parameter[1] == 'q') # Warn about quick versions
         @warn "Quick functions are not supported by timeRange boddatorient"
         result = NaN*ones(3, 1)
     else # call normal versions
         if(parameter == "pole")
             result = pole(body, timeRange, dict)
         elseif(parameter == "pm")
             result = pm(body, timeRange, dict)
         elseif(parameter == "eqx")
             result = eqx(body, timeRange, dict)
         elseif(parameter == "poledcm")
             result = poledcm(body, timeRange, dict)
         elseif(parameter == "pmdcm")
             result = pmdcm(body, timeRange, dict)
         elseif(parameter == "eqxdcm")
             result = eqxdcm(body, timeRange, dict)
         else
             @warn "Parameter not recognized"
             result = NaN*ones(3, 1)
         end
     end

     if(save)
         savedict(dict)
     end

     return result
end

"""
	a = getfn(fn, b, dict)

Relates desired input to designated function resulting in desired output.
"""
function getfn(fn::String,b::Int64,dict::Dict{String,Any})
    a = []
    if (fn == "gm") || (fn=="GM")
        a = getgm(b,dict)
    elseif fn == "radius"
        a = getradius(b,dict)
    elseif (fn=="j2") || (fn=="J2")
        a = getj2(b,dict)
    elseif (fn=="rot") || (fn=="rotation_period") || (fn=="ROT")
        a = getrotperiod(b,dict)
    elseif (fn=="occ")
        a = getocc(b,dict)
    elseif (fn=="typ") || (fn=="type") || (fn=="types")
        a = gettype(b,dict)
    elseif (fn=="ref") || (fn=="reference")
        a = getreference(b,dict)
    elseif (fn=="clo") || (fn=="close") || (fn=="close_approach")
        a = getclose(b,dict)
    elseif (fn=="dat") || (fn=="date")
        a = getdate(b,dict)
    elseif (fn=="mag") || (fn=="magnitude")
        a = getmagnitude(b,dict)
    end
    return a
end
