"""
	h = getmagnitude(body, dict)

Retrieves the absolute magnitude of the given body.
"""
function getmagnitude(body::Integer,dict::Dict{String,Any})::Float64
    if body < 1e6 #only for small bodies
        @warn "Cannot retrieve magnitude of major body $body"
        return NaN
    end

    h = getsymbol(body,"h",dict) #check saved data (& no ssd flag)
    if isempty(h)
        bodydata,_ = getsb(body,dict)  #check ssd
        if !isnan(bodydata[1])
            h = bodydata[6]
        else
            @warn "Magnitude not found for body $body"
            return NaN
        end
    end

    return h
end

"""
	occ = getocc(body, dict)

Retrieves the orbital condition code of the given body.
"""
function getocc(body::Integer,dict::Dict{String,Any})::Integer
    if body < 1e6 #only for small bodies
        @warn "Cannot retrieve orbital condition code of major body $body"
        return -1
    end

    occ = getsymbol(body,"occ",dict) #check saved data (& no ssd flag)
    if isempty(occ)
        bodydata,_ = getsb(body,dict)#check ssd
        if !isnan(bodydata[1])
            occ = bodydata[7]
        else
            @warn "No orbital condition code found for body $body"
            return -1
        end
    end

    return occ
end

"""
	typ = gettype(body, dict)

Retrieves the spectral type of the given body.
"""
function gettype(body::Integer,dict::Dict{String,Any})::AbstractString
    if body < 1e6
        @warn "Cannot retrieve spectral type for major body $body"
        return ""
    end
    typ = getsymbol(body,"types",dict)#check saved data (& no ssd flag)
    if isempty(typ)
        bodydata,_=getsb(body,dict)#check ssd
        if !isnan(bodydata[1])
            typ = bodydata[8]
        else
            @warn "Spectral type not found for body $body"
            return ""
        end
    end#if
    typ = strip(typ,'_')#underscore denotes checked but nothing found

    return typ
end

"""
	approach_data = getclose(body, dict)

Retrieves close approach data for the given body. It is returned in the following
format: [Date of closest approach (J2000)] [Other body] [Nominal distance (AU)]
Each row is one data entry.
"""
function getclose(body::Integer,dict::Dict{String,Any})::AbstractArray{Any, 2}
    if body < 1e6 #only for small bodies
        @warn "Cannot retrieve close approach data for major body $body"
        return Array{Any, 2}(undef, 0, 0)
    end

    s = download("https://ssd.jpl.nasa.gov/sbdb.cgi?sstr="*string(body)*";cad=1")
    f = open(s); rd = readlines(f); close(f)
    rm(s)

    cls =[]
    clt = []
    for ii in 1:length(rd)
        (occursin("Close-Approach Data",rd[ii])) && (cls = [cls; ii])
        (occursin("<b>Close-Approach Data",rd[ii])) && (clt = [clt;ii])
    end

    if isempty(cls)
        @warn "No close approach data for body $body"
        return Array{Any, 2}(undef, 0, 0)
    end

    (length(clt)>1) && (@warn "Close-approach data might not be correct")
    clt = clt[1]
    ii = findall(x->(x==clt),cls)
    tst = cls[ii[1]]
    ten = cls[ii[1]+1]
    lne = []

    for ii in tst:ten
        (occursin("</tr>",rd[ii])) && (lne = [lne; ii])
    end

    cad = zeros(length(lne)-3,3)
    LastNum = ""
    for ii in 1:length(lne)-3
        lns = match(r"size=\"-2\">",rd[lne[ii]+1])
        ln = match(r"</font>",rd[lne[ii]+1])
        # Protect against faulty matches, such as those caused by the alternate
        # designations table
        if (lns == nothing || ln == nothing)
            cad = cad[1:size(cad, 1) - 1,:]
            continue
        end

        T1 = rd[lne[ii]+1][lns.offset+10:ln.offset-1]
        T1 = split(T1)
        (T1[1][6:8]=="Jan") && (T1[1]=T1[1][1:5]*"01"*T1[1][9:end])
        (T1[1][6:8]=="Feb") && (T1[1]=T1[1][1:5]*"02"*T1[1][9:end])
        (T1[1][6:8]=="Mar") && (T1[1]=T1[1][1:5]*"03"*T1[1][9:end])
        (T1[1][6:8]=="Apr") && (T1[1]=T1[1][1:5]*"04"*T1[1][9:end])
        (T1[1][6:8]=="May") && (T1[1]=T1[1][1:5]*"05"*T1[1][9:end])
        (T1[1][6:8]=="Jun") && (T1[1]=T1[1][1:5]*"06"*T1[1][9:end])
        (T1[1][6:8]=="Jul") && (T1[1]=T1[1][1:5]*"07"*T1[1][9:end])
        (T1[1][6:8]=="Aug") && (T1[1]=T1[1][1:5]*"08"*T1[1][9:end])
        (T1[1][6:8]=="Sep") && (T1[1]=T1[1][1:5]*"09"*T1[1][9:end])
        (T1[1][6:8]=="Oct") && (T1[1]=T1[1][1:5]*"10"*T1[1][9:end])
        (T1[1][6:8]=="Nov") && (T1[1]=T1[1][1:5]*"11"*T1[1][9:end])
        (T1[1][6:8]=="Dec") && (T1[1]=T1[1][1:5]*"12"*T1[1][9:end])
        T2 = split(T1[2],":")
        T2 = parse(Float64,T2[1])/24 + parse(Float64,T2[2])/60/24
        cad[ii,1] = Dates.datetime2julian(DateTime(T1[1]))-2451545. + T2

        lns = match(r"size=\"-2\">",rd[lne[ii]+3])
        ln = match(r"</font>",rd[lne[ii]+3])
        T1 = rd[lne[ii]+3][lns.offset+10:ln.offset-1]
        if (LastNum == T1)
            cad[ii,2] = cad[ii-1,2]
        else
            T2 = getspk(T1,dict)
            cad[ii,2] = T2[1]
        end
        LastNum = T1

        lns = match(r"size=\"-2\">",rd[lne[ii]+4])
        ln = match(r"</font>",rd[lne[ii]+4])
        cad[ii,3] = parse(Float64,rd[lne[ii]+4][lns.offset+10:ln.offset-1])
    end
    # find <tr> in beginning of line and </tr>\n at end of lines
    # Looking at first, third, and fourth column (date, body, nominal dist)
    return cad
end
