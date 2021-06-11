# Copyright (c) 2020 California Institute of Technology (“Caltech”) and Embry Riddle Aeronautical University.
# U.S. Government sponsorship acknowledged.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
#  • Redistributions of source code must retain the above copyright notice, this 
#    list of conditions and the following disclaimer.
#  • Redistributions in binary form must reproduce the above copyright notice, 
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





mutable struct MB
    Names::AbstractArray{String}
    numbers::AbstractArray{Int64}
    gm::AbstractArray{Float64}
    rad::AbstractArray{Float64}
    j2::AbstractArray{Float64}
    rot::AbstractArray{Float64}
    ephref::AbstractArray{String}
    ephdate::AbstractArray{Float64}
    pp::AbstractArray{AbstractArray{Float64, 1}, 1}
    tri::AbstractArray{AbstractArray{Float64, 1}, 1}
end

mutable struct SB
    Names::AbstractArray{String}
    numbers::AbstractArray{Int64}
    gm::AbstractArray{Float64}
    rad::AbstractArray{Float64}
    rot::AbstractArray{Float64}
    h::AbstractArray{Float64}
    occ::AbstractArray{Float64}
    types::AbstractArray{String}
    ephref::AbstractArray{String}
    ephdate::AbstractArray{Float64}
    pp::AbstractArray{AbstractArray{Float64, 1}, 1}
    tri::AbstractArray{AbstractArray{Float64, 1}, 1}
end

"""
	gendict!(dict, opts)

Populates a dictionary for use with boddat commands.
"""
function gendict!(dict::Dict{String,Any}, opts::AbstractArray{Bool}=Bool[])
    if isempty(dict)
        sb = SB([],[],[],[],[],[],[],[],[],[],[],[])
        mb = MB([],[],[],[],[],[],[],[],[],[])
        tephf = TEPHF([],[],[],[])
        dict["sb"] = sb
        dict["mb"] = mb
        dict["tephf"] = tephf
    end

    #default options is [true false true], replace with any input and save in params
    optsInt=[true false true]
    optsInt[1:length(opts)] = opts
    param(dict,"save",optsInt[1])
    param(dict,"ssd",optsInt[2])
    param(dict,"keeps",optsInt[3])
    #boddat.m directory, tells where to save data
    if isempty(param(dict,"bdir"))
        bdir = pwd()
        param(dict,"bdir",bdir)
    end
    #parameter outputs, save in params as character array
    if isempty(param(dict,"bva"))
        ov = ["Pole","PM","Eqx"]
        bva = ["R","V","X","A","RTN","number","CB","triaxial", "Name","GM","radius",
        "j2","rotation_period","reference","date","mag","occ","types",
        "close_approach"]
        bva = [bva;ov]
        for jj in 1:3
            for ii in 1:length(ov)
                (jj == 1) && (ov[ii] = ov[ii]*"dcm")
                (jj == 2) && (ov[ii] = "q"*strip(ov[ii],['d','c','m']))
                (jj == 3) && (ov[ii] = ov[ii]*"dcm")
            end
            bva = [bva;ov]
        end
        param(dict,"bva",bva)
    end
    param(dict,"mbmod",false)
    param(dict,"sbmod",false)

    return
end

"""
	savedict(dict)

Saves the dictionary to a file. Overwrites existing file!
"""
function savedict(dict::Dict{String,Any})
    @static Sys.iswindows() ? (fn = string(param(dict,"bdir"), "\\bodydata.jld2")) : (fn=
    string(param(dict,"bdir"),"/bodydata.jld2"))
    if !isfile(fn)
        boddict = Dict()
        boddict["mb"] = MB([],[],[],[],[],[],[],[],[],[])
        boddict["sb"] = SB([],[],[],[],[],[],[],[],[],[],[],[])
        jldopen(fn,"w") do fid # Must use HDF5 and JLD package
            write(fid, "boddict", boddict)
        end
    end

    boddict = Dict()
    if param(dict,"mbmod")
        mb,_=getsbmb("mb",dict)
        boddict["mb"] = mb
        jldopen(fn,"r+") do fid
            if haskey(fid,"boddict")
                dT = read(fid, "boddict")
                boddict["sb"] = dT["sb"]
            end
        end
        jldopen(fn,"w") do fid
            write(fid, "boddict", boddict)
        end
        param(dict, "mbmod", false)
    end
    if param(dict,"sbmod")
        sb,_=getsbmb("sb",dict)
        boddict["sb"] = sb
        jldopen(fn,"r+") do fid
            if haskey(fid,"boddict")
                dT = read(fid, "boddict")
                boddict["mb"] = dT["mb"]
            end
        end
        jldopen(fn,"w") do fid
            write(fid, "boddict", boddict)
        end
        param(dict, "sbmod", false)
    end

    return
end

"""
	tryloaddict!(dict)

Attempts to load the previously saved dictionary from file.
"""
function tryloaddict!(dict::Dict{String,Any})
    if isempty(dict)
        gendict!(dict)
    end

    @static Sys.iswindows() ? (bdf = string(param(dict,"bdir"),"\\bodydata.jld2")) :
    (bdf = string(param(dict,"bdir"),"/bodydata.jld2"))
    if isfile(bdf)
        boddict = load(bdf,"boddict")
        dict["mb"] = boddict["mb"]
        dict["sb"] = boddict["sb"]
    else
        @warn "No dict file found"
    end

    return
end

"""
	x = getsymbol(body, symbol, dict)

Retrieves the requested symbol from the dictionary for the given body.
"""
function getsymbol(body::Int64,symbol::String,dict::Dict{String,Any})
    #get data from mb or sb structure
    x,_=getsbmb(body, dict)
    ii=findfirst(y -> y == body, x.numbers) #get datastructure and see if body is there

    if(ii == nothing)
        return [] # We've got no data
    end
    (length(ii)==1) && (ii = ii[1])

    o = []
    sym = Symbol(symbol)
    if !isempty(ii) && !param(dict,"ssd") && findfirst(y -> y == sym, fieldnames(typeof(x))) != nothing
        #see if should skip saved data and if datafield exists
        #see if entry exists in field for body

        field = getfield(x, sym)
        type = eltype(field)

        if isempty(field)
            return o
        elseif (type==Float64 || type==Int64) && length(field)>=ii
            if size(field,2) <= 1
                (!isnan(field[ii])) && (o=field[ii])
            else
                if any(.!isnan.(.!isnan.(field[ii,:])))
                    o=field[ii,:]
                end
            end
        elseif type == AbstractArray{Float64,1} && length(field)>=ii && !any(isnan, field[ii])# :tri and :pp
            o = field[ii]
        elseif eltype(getfield(x, sym)) == String
            if length(getfield(x, sym)) >= ii
                if !isempty(getfield(x, sym)[ii])
                    o = getfield(x, sym)[ii]
                end
            end
        end
    end
    return o
end

"""
	data = getsbmb(identifier, dict)

Gets small body or major body structure.
"""
function getsbmb(n::Integer, dict::Dict{String,Any})
    if n < 1e6 #convert number input to sb mb
        return getsbmb("mb", dict)
    else
        return getsbmb("sb", dict)
    end
end

function getsbmb(bi::String,dict::Dict{String,Any})
    data = param(dict,bi)
    if data == [] && bi == "mb" #read from saved data
        mb = MB([],[],[],[],[],[],[],[],[],[])
        param(dict,"mb",mb)
    elseif data == [] && bi == "sb"#read from saved data
        sb = SB([],[],[],[],[],[],[],[],[],[],[])
        param(dict,"sb",sb)
    elseif isempty(data.Names) || isempty(data.numbers)
        @static Sys.iswindows() ? (bdf = string(param(dict,"bdir"),"\\bodydata.jld2")) :
        (bdf = string(param(dict,"bdir"),"/bodydata.jld2"))
        if isfile(bdf)
            boddict = load(bdf,"boddict") # this will give "md" or "sd" elements
            data = boddict[bi]
            param(dict, bi, data)
        end
    end

    return data,bi
end

"""
	x = putsbmb(n, f, v, dict)

Writes data to structure.
"""
function putsbmb(n,f,v,dict::Dict{String,Any})
    #write data to structure
    x,bi=getsbmb(n[1],dict) #get structure
    bb=findall(z->(z==n[1]),x.numbers)
    if isempty(bb)&&length(n)>1
        n=n[2]
        bb=findall(z->(z==n),x.numbers)
    else
        n=n[1]
    end
    if isempty(bb) #see if body exists in field
        bb=length(x.numbers)+1
        f=["numbers" f]
        (size(v,2) > size(v,1)) ? (v = [n v]) : (v = [n; v])
    end

    if typeof(f) == String
        flds = 1
    else
        flds = length(f)
    end

    for ii in 1:flds
      #different fields
        Nm = Symbol("T")
        if flds == 1
            Nm = Symbol(f)
        else
            Nm = Symbol(f[ii])
        end
        if typeof(v[ii])==String || (typeof(v[ii])==SubString{String})
            #pad with spaces, then write
            if isempty(getfield(x, Nm))
                init = Array{String}(undef, length(x.numbers));
                for jj in 1:length(x.numbers)
                    init[jj] = "";
                end
                setfield!(x, Nm, init)
            end
            ct = 1
            nn=size(getfield(x,Nm),1)+1
            for jj in bb
                if iszero(bb[ct]- nn)
                    push!(getfield(x,Nm),v[ii])
                elseif bb[ct]>nn
                    for kk in 1:bb[ct]-nn
                        push!(getfield(x,Nm)[nn-1+kk],"_")
                    end
                    push!(getfield(x,Nm),v[ii])
                else
                    getfield(x,Nm)[jj]=v[ii]
                end
                ct = ct + 1
            end
        else
            #write, fill with nans
            if isempty(getfield(x, Nm))

                if Nm == :numbers
                    setfield!(x,Nm,zeros(Int64,length(v[ii])))
                elseif Nm == :tri || Nm == :pp
                    init = Array{Array{Float64, 1}}(undef, length(x.numbers))
                    for jj in 1:length(x.numbers)
                        init[jj] = NaN * zeros(length(v[ii]))
                    end
                    setproperty!(x, Nm, init)
                else
                    (isempty(x.numbers)) && (setfield!(x,Nm,NaN*zeros(length(v[ii]))))
                    (!isempty(x.numbers))&&(setfield!(x,Nm,NaN*zeros(length(x.numbers))))
                end

            end
            nn=size(getfield(x,Nm),1)+1
            ct = 1
            for jj in bb
                if bb[ct] >= nn
                    # Pad out dict if we're trying to store a value for a body we
                    # haven't populated up to yet. This usually happens with small
                    # bodies since the dict is of unknown length.
                    if Nm == :tri || Nm == :pp
                        for kk in 1:bb[ct]-nn
                            push!(getfield(x,Nm), NaN * zeros(length(v[ii])))
                        end
                    else
                        for kk in 1:bb[ct]-nn
                            push!(getfield(x,Nm)[nn-1+kk],NaN)
                        end
                    end

                    if typeof(getfield(x,Nm)[jj-1][1])==typeof(v[ii]) || Nm == :tri || Nm == :pp
                        push!(getfield(x,Nm),v[ii])
                    else
                        push!(getfield(x,Nm),convert(typeof(getfield(x,Nm)[jj-1][1]),v[ii]))
                    end
                else
                    if typeof(getfield(x,Nm)[jj][1])==typeof(v[ii])
                            getfield(x,Nm)[jj]=v[ii]
                    elseif Nm == :tri || Nm == :pp
                        getfield(x,Nm)[jj] = convert(typeof(getfield(x,Nm)[jj]), v[ii])
                    else
                        if flds == 1 && typeof(v) == String
                            getfield(x,Nm)[jj]=convert(typeof(getfield(x,Nm)[jj][1]),Meta.parse(v))
                        elseif eltype(v) == String
                            getfield(x,Nm)[jj]=convert(typeof(getfield(x,Nm)[jj][1]),Meta.parse(v[ii]))
                        else
                            getfield(x,Nm)[jj]=convert(typeof(getfield(x,Nm)[jj][1]),v[ii])
                        end
                    end
                end
                ct = ct + 1
            end
        end
    end #numel(f)
    (param(dict,bi*"mod")!= true) && (param(dict,bi*"mod",true))
    param(dict,bi,x)#save
    return x
end

"""
	o = param(vs, key, value)

Allows a variety of values to be stored in RAM under one name and sets value
to dictionary structure.
"""
function param(vs::Dict,key::String,value)
    # Allows a variety of values to be stored in RAM under one name
    # Sets value to dictionary structure
    vs[key] = value
end

"""
	o = param(vs, key)

Recalls value from dictionary structure made with param.
"""
function param(vs::Dict,key::String)
    # Recalls value from dictionary structure
    # If not in dictionary, returns empty matrix
    if !haskey(vs,key)
        o = []
    else
        o = vs[key]
    end
    return o
end
