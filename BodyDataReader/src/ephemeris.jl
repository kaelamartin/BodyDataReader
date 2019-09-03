@enum EphemType begin
    position = 1
    velocity = 2
    state = 3
    acceleration = 4
end

@enum TimeType begin
    single = 1
    multiple = 2
end

mutable struct TEPHF
    numbers::AbstractArray{String}
    f::AbstractArray{String}
    tl::AbstractArray{Float64}
    sb::AbstractArray{Float64}
end

"""
	reference = getreference(body, dict)

Retrieves the ephemeris reference for the given body. If no cached ephemeris
reference is found, this function will check the SSD database.
"""
function getreference(body::Int64,dict::Dict{String,Any})::AbstractString
    ref = getsymbol(body,"ephref",dict)#check saved data (& no ssd flag)
    if isempty(ref) || param(dict,"ssd") #see what current file is
        if body < 1e6
            tephf = getephtime(dict)
            ii = findall(x -> x==string(body),tephf.numbers)
            ref = tephf.f[ii[1]] #get ephemeris from list
        else
            bodydata,e_ = getsb(body,dict)
            if !isnan(bodydata[1])
                ref = e_
            else
                ref = ""
            end
        end#get reference from ssd
    end
    if isempty(ref)
        @warn "File of last ephemeris for $body is unknown"
        return ""
    end

    return ref
end

"""
	date = getdate(body, dict)

Retrieves the date that ephemeris data was retrieved for the given body. Only
returns data after an ephemeris command has been run.
"""
function getdate(body::Integer,dict::Dict{String,Any})::Float64
    date = getsymbol(body,"ephdate",dict) #check saved data (& no ssd flag), no ssd analogue
    if isempty(date)
        date = NaN
        @warn "Date of last ephemeris for $body is unknown"
    end
    return date
end

"""
	ephem!(body, time, type, dict, X)

Retrieves the requested ephemeris data for the given body at the given time.
Output is stored in X.
"""
function ephem!(body::Integer, timeRange::AbstractArray{Float64, 1}, type::EphemType, dict::Dict{String,Any}, X::AbstractArray{Float64, 2})
    for k in 1:length(timeRange)
        time = timeRange[k]
        X[:,k] = ephem(body, time, type, dict)
    end
end

function ephem!(body::Integer, time::Float64, type::EphemType, dict::Dict{String,Any}, X::AbstractArray{Float64, 1})
    X[:] = ephem(body, time, type, dict)
end

"""
	ephem!(bodies, time, type, dict, X)

Retrieves the requested ephemeris data for bodies[1] relative to bodies[2] at the given time.
Output is stored in X, which must have size 6 for R, V, or X, and size 12 for A.
"""
function ephem!(bodies::AbstractArray{Int64, 1}, timeRange::AbstractArray{Float64, 1}, type::EphemType, dict::Dict{String,Any}, X::AbstractArray{Float64, 2})
    if length(bodies) != 2
        @error "Bodies array must have length 2 for relative ephemeris"
        X[:,:] = NaN*ones(length(X))
    end

    for k in 1:length(timeRange)
        time = timeRange[k]
        X[:,k] = ephem(bodies[1], time, type, dict)
        X[:,k] -= ephem(bodies[2], time, type, dict)
    end
end

function ephem!(bodies::AbstractArray{Int64, 1}, time::Float64, type::EphemType, dict::Dict{String,Any}, X::AbstractArray{Float64, 1})
    if length(bodies) != 2
        @error "Bodies array must have length 2 for relative ephemeris"
        X[:] = NaN*ones(length(X))
    end

    X[:] = ephem(bodies[1], time, type, dict)
    X[:] -= ephem(bodies[2], time, type, dict)
end

"""
	Xarr = ephem(body, timeRange, type, dict)

Retrieves the requested ephemeris data for the given body over the timeRange.
Returns X.
"""
function ephem(body::Int64, timeRange::AbstractArray{Float64, 1}, type::EphemType, dict::Dict{String,Any})::AbstractArray{Float64, 2}
    # Get the right size for output
    if type == acceleration
        n = 12
    else
        n = 6
    end

    Xarr = Array{Float64}(undef, n, length(timeRange)) # Preallocate
    for k in 1:length(timeRange)
        time = timeRange[k]
        Xarr[:,k] = ephem(body, time, type, dict)
    end
    return Xarr
end

"""
	Xvec = ephem(body, time, type, dict)

Retrieves the requested ephemeris data for the given body at the given time.
Returns X, which has size 6 for R, V, or X, and size 12 for A.
"""
function ephem(body::Int64, time::Float64, type::EphemType, dict::Dict{String,Any})::AbstractArray{Float64, 1}
    # Get the right size for output
    if type == acceleration
        n = 12
    else
        n = 6
    end

    #get saved data
    d = param(dict,"ephem_data")
    if isempty(d)
        d = Dict()
    end

    if param(dict,"ssd")
        X,_,err = makeephem(body, [time], point, dict) #go straight to horizons

        #flag if any data is out of bounds
        (length(err)==1) && (err = [err false])

        if err[1]
            @warn "Minimum Horizons ephemeris time for body $body has been exceeded"
        elseif err[2]
            @warn "Maximum Horizons ephemeris time for body $body has been exceeded"
        end

        if type == A
            @warn "Acceleration and jerk data cannot be retrieved from horizons"
        end

        return X #return state
    end

    cb = getcb(body) #central body, treat planets and satellites differently
    if !haskey(d, body) # if saved data doesn't exist
        @static Sys.iswindows() ?
        (efile = param(dict,"bdir")*"\\ephem\\"*string(body)*".jld2") :
        (efile = param(dict,"bdir")*"/ephem/"*string(body)*".jld2")

        if (!isfile(efile)) && (body > 1e6)
            sbdata,_ = getsb(body,dict)
            body = sbdata[1] # If body numbers change, we'll still be using the old one
            @static Sys.iswindows() ?
            (efile = param(dict,"bdir")*"\\ephem\\"*string(body)*".jld2") :
            (efile = param(dict,"bdir")*"/ephem/"*string(body)*".jld2")
            #see if number changed
        end

        if !isfile(efile) #te is saved time span, di is data to save, don't save nans
            te = [-1, 18261.5, 0.1 + 0.4*(cb==10) + 0.5*(body > 1e6)]
            di,td,_ = makeephem(body, te, multiple, dict)
            di = [td;di]
            if isempty(di)
                return NaN*ones(n)
            end#no good data
        else
            f = load(efile)
            di = f["d"]
        end#file exists, read data, 1st two entries is data size

        if cb != 10 #gives better accuracy for Moons, but takes ~twice time

            td = di[1,:]
            R = di[2:4,:]
            V = di[5:7,:]
            H = Array{Float64}(undef, 3, size(R, 2))
            for ii in 1:size(R,2)
                T = cross(R[:,ii],V[:,ii])
                H[:,ii] = T/sqrt(dot(T,T))
            end
            H = mean(H,dims=2)
            H = H/sqrt(dot(H,H))

            h12 = sqrt(H[1]^2+H[2]^2)
            N = [-H[2]/h12; H[1]/h12; 0]
            Q = [-H[3]*N[2]; H[3]*N[1]; h12]
            dcm = [N Q H]

            R = transpose(dcm)*R
            zr = R[3,:]
            R = R[1:2,:]
            r = zero(zr)
            for ii in 1:size(R,2)
                r[ii] = sqrt(dot(R[:,ii],R[:,ii]))
            end
            qr1 = atan.(R[2,:],R[1,:])

            dqr = diff(qr1, dims=1)
            T1 = falses(length(dqr))
            for ii in 1:length(dqr)
                T1[ii] = abs(dqr[ii])>pi
            end

            jj = findall(T1)
            for ii in 1:length(jj)
                if (ii == length(jj))
                    qr1[jj[ii]+1:end] = qr1[jj[ii]+1:end] .+ ii*2*pi
                else
                    qr1[jj[ii]+1:jj[ii+1]] = qr1[jj[ii]+1:jj[ii+1]] .+ ii*2*pi
                end
            end

            V = transpose(dcm)*V
            zv = V[3,:]
            V = V[1:2,:]
            v = zero(zr)
            qv = zero(zr)
            for ii in 1:size(V,2)
                v[ii] = dot(V[:,ii],R[:,ii])/r[ii]
                qv[ii] = (R[1,ii]*V[2,ii]-R[2,ii]*V[1,ii])/r[ii]^2
            end
            di[2:end,:] = [r';qr1';zr';v';qv';zv']
        else
            dcm = []
        end

        d = param(dict,"ephem_data")
        (isempty(d)) && (d = Dict{Int64,Dict}())
        d[body] = Dict{String,Array{Float64}}()

        if param(dict,"keeps") #save spline
            RVspline6!(di,d[body])
            d[body]["dcm1"] = dcm
            param(dict,"ephem_data",d)
        else
            RVAspline!(di,d[ii])
            d[body]["dcm1"] = dcm
            param(dict,"ephem_data",[])
        end

    end

    ist = false
    #see if data is out of bounds
    tl = d[body]["breaks"][1]
    if time < tl
        ist = true
        @warn "Minimum time for body $body ephem spline is $(Dates.julian2datetime(tl+2451545)), retrieving from Horizons instead"
    end
    tl = d[body]["breaks"][end]
    if time > tl
        ist = true
        @warn "Maximum time for body $body ephem spline is $(Dates.julian2datetime(tl+2451545)), retrieving from Horizons instead"
    end

    if !ist # We're on the spline

        X = zeros(n)
        if type != velocity || cb != 10 # Populate position
            pval!(d[body], [time], 6, view(X, 1:3))
        end
        if type != position # Populate velocity
            pval!(d[body], [time], 5, view(X, 4:6))
        end
        if type == acceleration # Populate acceleration and jerk
            pval!(d[body], [time], 4, view(X, 7:9))
            pval!(d[body], [time], 3, view(X, 10:12))
        end

        if cb != 10#gives better accuracy for Moons, but takes ~twice time
            dcm = d[body]["dcm1"]
            c = cos(X[2])
            s = sin(X[2])
            r = X[1]

            if type != velocity # Fix position
                X[1:3] = dcm*[r*c; r*s; X[3]]
            end

            if type != position # Fix velocity
                dr = X[4]
                dq = X[5]
                rdq = r*dq
                X[4:6] = dcm*[dr*c - rdq*s; dr*s + rdq*c; X[6]]
            end

            if type == acceleration # Fix acceleration and jerk
                d2r = X[7]
                d2q = X[8]
                r_ = d2r-rdq*dq
                q_ = r*d2q + 2*dr*dq
                X[7:9] = dcm*[r_*c - q_*s; r_*s + q_*c; X[9]]
                r_ = X[10] - 3*dr*dq^2 - 3*rdq*d2q
                q_ = r*X[11] + 3*d2r*dq + 3*dr*d2q - rdq*dq^2
                X[10:12] = dcm*[r_*c - q_*s; r_*s + q_*c; X[12]]
            end
        end
    else # We're off the spline
        save = param(dict,"save")
        param(dict,"save",false)
        X,_,_ = makeephem(body, [time], single, dict)
        param(dict,"save",sav)
    end

    return X
end

"""
	X, tt_transpose, error = makeephem(body, time, type, dict, latlong)

Reads ephemeris data from Horizons. Takes a body to retrieve for, and two params
to specify time: time and type. The param type is an enum that differentiates
between two distinct behaviors:
    If type is single, data will be retrieved for each entry in time
    If type is multiple, data will be retrieved over a range, and time must be
        formatted as time = [start, end, delta]
All data is retrieved relative to the solar system barycenter, except if latlong
is specified.
If latlong is specified, then data is retrieved for body as the central body with
an observer at the coordinates lat, long. If used, latlong must have length 2.
"""
function makeephem(body::Integer, time::AbstractArray{Float64, 1}, type::TimeType, dict::Dict{String,Any},
                latlong::AbstractArray{Float64, 1}=Float64[])
    #read from horizons
    err = [false, false]
    #get data wrt central body

    if body < 1e3 || body > 1e6
        tephf = getephtime(dict) #get ephemeris time span

        if body > 1e6 && body < 2e6 #comets are dumb
            CAP = "#3BCAP"
        else
            CAP = ""
        end

        if body > 1e6 #time limits for small bodies
            tlim=tephf.sb
            DES="DES="
            ieph=1
        else #time limits for major bodies
            ieph = findall(x -> x == string(body), tephf.numbers)
            tlim=tephf.tl[ieph[1],:]
            DES=""
        end

        if type != single && type != multiple
            @warn "Unrecognized time type, assuming single"
            type = single
        end

        #tt is expected time output from horizons, can only run 400 times in list,
        # or 90,000 times in {start,end,delta}
        if type == single
            ist = sortperm(collect(time)) #dev add collect to force column vector
            tt = time[ist]
            mnt = 400
        elseif type == multiple
            tt = collect(time[1]:time[3]:time[2])
            mnt = 90000
            if length(time) != 3
                @warn "Invalid {start time, end time, delta time} input to Horizons"
            end
        end

        if (isempty(ieph))
            @warn "No horizons ephemeris file found for body $body"
            return NaN*ones(6, length(tt)), tt, [true]
        end

        if minimum(tt) < tlim[1]
            ii = findall(x -> x < tlim[1], tt)
            nt1 = length(ii)
            err[1] = true
        elseif minimum(tt) >= tlim[1]
            nt1 = 0
            err[1] = false
        end

        if maximum(tt) > tlim[2]
            jj = findall(x-> x > tlim[2], tt)
            nt2 = length(jj)
            err[2] = true
        else
            nt2 = 0
            err[2] = false
        end

        n2t = length(tt) - nt2
        t12 = round.(Int64,range(nt1, stop=n2t, length=ceil(Int64,(n2t-nt1)/mnt)+1))
        #t12 breaks time span into different runs
        if body > 1e6
            tnow = argmin(abs.(Dates.datetime2julian(now()) - 2451545 .- tt))
            t12 = sort(unique([tnow-1; t12]))
        end

        X = vec(zeros(7,0))
        n12 = length(t12)
        println("Retrieving data for Body ",body," from Horizons")
        for ti in 1:n12 .- 1#only do mnt at a time
        #time input for url
            if type == multiple
                tstr =
                Printf.@sprintf("&START_TIME='JD%%20%.9f'&STOP_TIME='JD%%20%.9f'&STEP_SIZE='%d'",
                tt[t12[ti]+1]+2451545,tt[t12[ti+1]]+2451545,t12[ti+1]-t12[ti]-1)
            else
                tstr = "&TLIST='"
                for ii in t12[ti]+1:t12[ti+1]
                    tstr = tstr*Printf.@sprintf("%.9f%%0A",tt[ii]+51544.5)
                end
                tstr = tstr*"'"
            end
            #set target as center of body and observer as ll (lat,lon) for orientation
            if !isempty(latlong)
                centralBody = body
                coord="coord"
                llstr=Printf.@sprintf("&COORD_TYPE='GEODETIC'&SITE_COORD='%.8f,%.8f,0'",latlong[1],latlong[2])
            else
                centralBody = 0
                coord=""
                llstr=""
            end

            # horizons URL
            url=Printf.@sprintf(
            "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='%s%d%s'",
            DES,body,CAP)
            url = url*Printf.@sprintf("&CENTER='%s@%d%s'&MAKE_EPHEM='YES'",coord,centralBody,llstr)
            url = url*Printf.@sprintf("&TABLE_TYPE='VECTORS'%s&OUT_UNITS='KM-S'",tstr)
            url = url*"&VECT_TABLE='2'&REF_PLANE='ECLIPTIC'&REF_SYSTEM='J2000'"
            url = url*"&VECT_CORR='NONE'&VEC_LABELS='NO'&CSV_FORMAT='NO'"

            s = download(url) #run horizons url command and extract data
            f = open(s); rd = readlines(f); close(f)

            # Horrific file parsing nightmare
            mx = length(rd); i1 = 0; T1 = true
            while T1
                i1 = i1 + 1
                (occursin("\$\$SOE",rd[i1],)) && (T1 = false)
                (i1 == mx) && (T1 = false)
            end

            T1 = true; i2 = copy(i1)
            (i2==mx) && (T1=false)
            while T1
                i2 = i2 + 1
                (occursin("\$\$EOE",rd[i2])) && (T1 = false)
                (i2 == mx) && (T1 = false)
            end

            clw = length(rd[1])
            if i1 != mx
                rd = 0
                try
                    X1 = readdlm(s,skipstart=i1,use_mmap=false,skipblanks=true,
                    dims=(mx,clw))
                    X1 = X1[1:i2-i1-1,1:3]
                    XT = Array{Float64}(undef, 7)
                    for ii in 1:floor(Int64,size(X1,1)/3)
                        XT[1] = X1[3*ii-2,1]
                        XT[2] = X1[3*ii-1,1]
                        XT[3] = X1[3*ii-1,2]
                        XT[4] = X1[3*ii-1,3]
                        XT[5] = X1[3*ii,1]
                        XT[6] = X1[3*ii,2]
                        XT[7] = X1[3*ii,3]
                        append!(X,XT)
                    end
                catch
                    X1 = readlines(s)
                    X1 = X1[i1+1:i2-1]
                    nl = div(i2-i1-1,3)
                    XT = Array{Float64}(undef, 7)
                    for ii in 1:nl
                        XT[1] = Meta.parse(collect(m.match for m = eachmatch(r"[^ =]+",X1[3*ii-2],overlap = overlap))[1])
                        T1 = collect(m.match for m = eachmatch(r"[^ =]+",X1[3*ii-1],overlap=false))
                        XT[2] = Meta.parse(T1[2])
                        XT[3] = Meta.parse(T1[4])
                        XT[4] = Meta.parse(T1[6])
                        T1 = collect(m.match for m = eachmatch(r"[^ =]+",X1[3*ii],overlap=false))
                        XT[5] = Meta.parse(T1[2])
                        XT[6] = Meta.parse(T1[4])
                        XT[7] = Meta.parse(T1[6])
                        append!(X,XT)
                    end
                end
              # If we made it here, we've succeeded
            else
                @warn "Horizons error"
                return NaN*ones(6, length(tt)), tt, err
            end

            if ti != n12-1
                println("Horizons progress: ",floor(Int64,(t12[ti+1]-nt1)/(n2t-nt1)*100),"%")
            end

            rm(s) # Clean up temp files

        end#ti

        X = reshape(X,7,:)
        if iszero(nt1) && iszero(nt2)
            tt = X[1,:] .- 2451545.
        elseif nt1 == 0
            tt = [X[1,:] .- 2451545.; tt[end-nt2+1:end]]
        elseif iszero(nt2)
            tt = [tt[1:nt1]; X[1,:] .- 2451545.]
        else
            tt = [tt[1:nt1]; X[1,:] .- 2451545.; tt[end-nt2+1:end]]
        end

        X=[NaN*ones(6,nt1) (1-2*(!isempty(latlong)))*X[2:7,:] NaN*ones(6,nt2)]
        #pad out of range data with nans, if orient X=-X
        if type == single
            X[:,ist] = X
            tt[ist] = tt
        end#output times in same order as input

    else #Lagrange point
        if type == single
            tt = time
            ist = collect(1:length(time))
        else
            tt = collect(time[1]:time[3]:time[2])
        end
        X = lagrange(body, tt, dict)
    end

    @static Sys.iswindows() ? (edir = string(param(dict,"bdir"),"\\ephem")) :
    (edir = string(param(dict,"bdir"),"/ephem"))

    if !isdir(edir) #make "ephem" directory
        mkdir(edir)
    end

    if isempty(latlong) && param(dict,"save") #save ephemeris data
        @static Sys.iswindows() ? (f = edir*"\\"*string(body)*".jld2") :
        (f = edir*"/"*string(body)*".jld2")

        jldopen(f,"w") do fid
            tt_trans = transpose(tt)
            d = [tt_trans;X]
            if type == single
                d = d[:,ist]
            end
            if any(isnan.(X[1,:])) #save times in order, don't save nans
                rmc = Int64[]
                for ii in 1:size(X,2)
                    (!isnan(X[1,ii])) && (rmc = [rmc; ii])
                end
                d = copy(selectdim(d,2,rmc))
            end
            # write(fid,"sz",size(d)) # don't need size
            write(fid,"d",d)
        end

        if body > 1e6 #always update small bodies
            n,e_ = getsb(body,dict)
            putsbmb(body,["ephref" "ephdate"],
            [e_ Dates.datetime2julian(now())-2451545.],dict)
        elseif body < 1e3 #update reference file and time
            putsbmb(body,["ephref" "ephdate"],
            [tephf.f[ieph,:] Dates.datetime2julian(now())-2451545.],dict)
        else
            putsbmb(body,["ephref" "ephdate"],
            [getreference([floor(Int64,n/10)],dict) Dates.datetime2julian(now())-2451545.],
            dict)
        end
    end

    return X, transpose(tt), err
end

"""
	X = lagrange(point, time, dict)

Computes ephemeris for inputted Lagrange point.
"""
function lagrange(point::Int64,time::AbstractArray{Float64},dict::Dict{String,Any})
    i = mod(point,10) #Lagrange pt number
    body = floor(Int64,point/10) #central body
    cb = getcb(body)
    m1 = getgm([cb],dict)
    m2 = getgm([body],dict) #masses
    m1 = m1[1]; m2 = m2[1]
    #calculate distance for circular (a), then scale by distance (r)
    if i < 4 #L1,L2,L3
        #non-dim
        m=m1+m2
        m1=m1/m
        m2=m2/m
        m1i=m1
        m2i=m2

        #account for direction of pull
        (i==1||i==3) && (m2i=-m2)
        (i==3) && (m1i=-m1i)
        if i < 3 #initial guess
            L = exp(log(m2/m1/3)/3)
            (i==1) && (L=-L)
        elseif i == 3
            L = -2.
        end
        L=L+1e-42im #complex step
        e_ = 0.
        for cc=1:99
            e_=m1+L-m1i/(L+1.)^2-m2i/L^2#balance dynamics, rotation - gravity
            de=imag(e_)/1e-42
            e_=real(e_)
            (abs(e_)<9*eps()) && (break)
            #dL=e./de;ii=abs(dL)>xL;if any(ii);dL(ii)=sign(dL(ii))*xL;end
            L=L-e_/de
        end
        L=real(L)
        (abs(e_)>9*eps()) && (@warn "Lagrange not converged")
        #e=m1+L-m1i./(L+1).^2-m2i./L.^2;dL=-imag(e)/1e-42 ./de
        #L=L+dL*1e-42i;X=L*ephem(b,t+1e-42i,1);X=[real(X);imag(X)/1e-42/86400];
        X,_ = ephem(body,t,state,dict)
        X=L*X #point is along position of secondary wrt primary

    else#L4,L5
        #triangular points
        X,_ = ephem(body,t+1e-42im,state,dict)
        R=X[1:3,:]
        V=X[4:6,:]

        for ii in 1:size(R,2)
            y=1. #equilateral
            H=cross(R[:,ii],V[:,ii])
            HR=cross(H,R[:,ii])
            r = sqrt(sum(R.*R))
            hr = sqrt(sum(HR.*HR))
            #y=exp(-log(sum(H.^2)./r/(m1+m2))*2/3);#adjust height of triangle?
            Y = HR/hr
            y=r*sqrt(y-.25) #distance from axis connecting bodies
            (i==5) && (y=-y) #L5 trails
            T1=-R[:,ii]/2. + Y*y
            X[:,ii]= [real(T1); imag(T1)/1e-42/86400] #complex step
        end
        X = real(X)
    end#i
    return X
end

"""
	ed = getephdata(body, dict)

Retrieves ephemeris header data from JPL's Solar System Dynamics database (SSD).
"""
function getephdata(body::Int64,dict::Dict{String,Any})
    #ephemeris header data
    ef1=getephtime(dict)
    ii = findall(x->(x==string(body)),ef1.numbers)
    ef = lowercase(match(r"^[^-._ ]*",ef1.f[ii[1]]).match)
    #only need first file if merged
    ef = convert(String,ef)

    if isempty(ef)
        @warn "No ephemeris header file found for $body"
        ed=[]
        return ed
    end

    ed=param(dict,ef) #see if file data already saved
    if (isempty(ed)) && (body < 400) #inner planets and barycenters
        #read in DE### header data
        ef=match(r"de\d+",ef).match #deXXX format
        ef = convert(String,ef)
        efdir="ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/"*ef*"/" #directory

        client = FTP(efdir)
        dir = readdir(client)
        header = dir[end-1] # Header is always second to last file
        s = download(client, header)
        rd = read(s, String); close(s)
        close(client)

        ii= something(findfirst("GROUP   1040",rd), 0:-1)
        jj = something(findfirst("GROUP   1041",rd), 0:-1)
        kk = something(findfirst("GROUP   1050",rd), 0:-1)

        T1 = split(rd[ii[end]+1:jj[1]-1])  #read in variable names (each 6 chars)
        vars = Array{AbstractString}(undef, length(T1)-1)
        for ii in 1:length(T1)-1
            vars[ii] = T1[ii+1]*""
        end

        T1 = split(rd[jj[end]+1:kk[1]-1])  #read in variable names (each 6 chars)
        vals = Array{Float64}(undef, length(T1)-1)
        for ii in 1:length(T1)-1
            TD = something(findfirst("D",T1[ii+1]), 0:-1)
            T11 = parse(Float64,T1[ii+1][1:TD[1]-1])
            T1E = parse(Float64,T1[ii+1][TD[1]+1:end])
            vals[ii] = T11*10 .^T1E
        end

        ed = Dict()
        ed["vars"]=vars
        ed["vals"]=vals
        #write gm values
        au=ed["vals"][uniquestr(ed["vars"],"AU",1)]; au = au[1]
        ed["gm10"]=ed["vals"][uniquestr(ed["vars"],"GMS")]*au^3/86400 .^2 #sun
        ed["gm3"]=ed["vals"][uniquestr(ed["vars"],"GMB")]*au^3/86400 .^2 #E-M bary
        emr=ed["vals"][uniquestr(ed["vars"],"EMRAT")]; emr = emr[1]#E/M ratio
        ed["gm301"]=ed["gm3"]/(1 .+emr);ed["gm399"]=ed["gm3"]*emr/(1 .+emr);
        #Moon and Earth barycenters
        ed["gm0"]=ed["gm10"]+ed["gm3"]

        for ii in [1:2;4:9]
            GM="gm"*string(ii)
            gm=ed["vals"][uniquestr(ed["vars"],"GM"*string(ii))]*au^3/86400^2
            ed[GM] = gm
            ed["gm0"]=ed["gm0"]+gm
        end

        ed["gm199"]=ed["gm1"];ed["gm299"]=ed["gm2"];#Mercury and Venus are Trouble
        #radii and J2
        ed["rad10"]=ed["vals"][uniquestr(ed["vars"],"ASUN")]
        ed["j2_10"]=ed["vals"][uniquestr(ed["vars"],"J2SUN")]
        ed["rad399"]=ed["vals"][uniquestr(ed["vars"],"RE")]
        ed["j2_399"]=ed["vals"][uniquestr(ed["vars"],"J2E",1)]
        ed["rad301"]=ed["vals"][uniquestr(ed["vars"],"AM",1)]
        ed["j2_301"]=ed["vals"][uniquestr(ed["vars"],"J2M")]
        ed["rad199"]=ed["vals"][uniquestr(ed["vars"],"RAD1")]
        ed["rad299"]=ed["vals"][uniquestr(ed["vars"],"RAD2")]
        param(dict,ef,ed)#save data

    elseif isempty(ed) #outer planets and satellites
        #check default directory and eph filename

        client = FTP("ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/")
        dir = readdir(client)
        rd = join(dir, " ")

        if rd[1]!='<'
            fil = collect(m.match for m = eachmatch(r"\S+(?=.txt)",rd[1:end],overlap=false))
            fils = union(fil)
        else
            fil = collect(m.match for m = eachmatch(r">\S+(?=.txt)",rd[1:end],overlap=false))
            fils = union(fil)
            for ii in 1:length(fils)
                fils[ii] = strip(fils[ii],'>')
            end
        end
        ii = uniquestr(fils,ef,0)
        if isempty(ii)
            jj = findall(x->(occursin(ef,x)),fils)
            (isempty(jj)) && (ii = uniquestr(fils,ef))
            (isempty(ii)) ? (err = false) : (err = true)
            if isempty(ii) # Look for the first 3 leters of the highest # file
                fils2 = Array{eltype(fils)}(size(fils))
                nf = length(fils2)
                for kk in 1:nf
                    fils2[kk] = fils[nf-kk+1]
                end
                    ii = uniquestr(fils2,ef[1:3])
                    ii = nf - ii + 1 #
                    (isempty(ii)) ? (err = false) : (err = true)
            end
            #could also check
            #ftp://ssd.jpl.nasa.gov/pub/eph/satellites/rckin/rckin.*.log,
            # but GM likely 0
        else
            err = true
        end

        if err

            s =  download(client, fils[ii[1]]*".txt")
            rd = read(s, String); close(s)

            #s is ephemeris header file
            #get GM from "Bodies on the File" table: skip space,keep 3 #s,
            # skip some space, keep some #s, decimal, & non-space,skip space
            T1 = true; fils = [];
            fil = collect(m.match for m = eachmatch(r"\s(\d{3})\s+((\d+\.\S*))\s",rd,overlap=false))
            fils = union(fil)
            gm = Array{AbstractString}(undef, length(fils),2)
            for ii in 1:length(fils)
                gm[ii,:] = split(fils[ii])
            end

            ed = Dict()
            for ii in 1:size(gm,1)
                ed["gm"*gm[ii,1]] = parse(Float64,gm[ii,2])
            end
            #get J2 and radius of planet
            T1 = something(findfirst("J"*string(floor(Int64,body/100))*"02", rd), 0:-1)
            T2 = something(findfirst(r" [PJ]",rd[T1[end]:end]), 0:-1)
            ed["j2_"*string(convert(Int64,body))] =
            parse(Float64,strip(rd[T1[end]+1:T2[1]+T1[end]-2]))
            T1 = something(findfirst("RADIUS",rd), 0:-1)
            T2 = something(findfirst("J",rd[T1[end]:end]), 0:-1)
            ed["rad"*string(convert(Int64,body))] =
            parse(Float64,strip(rd[T1[end]+1:T2[1]+T1[end]-2]))

        else
            @warn "ephemeris header file $ef for $body not found"
        end#err
        param(dict,ef,ed)
		close(client)
    end
    return ed
end

"""
	tephf = getephtime(dict)

Reads in time span data from Horizons and converts it to a Julian date.
"""
function getephtime(dict::Dict{String,Any})
    #Ephemeris file and time spans used by Horizons
    tephf=param(dict,"tephf")
    if isempty(tephf.numbers)
        s = download("https://ssd.jpl.nasa.gov/eph_spans.cgi?id=A") #Planets
        f = open(s); rd = read(f, String); close(f)
        rm(s)

        #Read in bodynumber, begin time, " not " or " to " flag, end time, and file
        ssT = Array{AbstractString}(undef, 1, 5)
        off = 1; T1 = true; sss = Array{AbstractString}(undef, 0)
        while T1
            ss=match(
            r"<td.*?>(\d+)</td>.*?<\w\w?>(.*?) (not|to) (.*?)</\w\w?>&nbsp;.*?&nbsp;(.*?)&nbsp;",
            rd[off:end])
            if ss == nothing
                T1 = false
            else
                off = ss.offsets[end]+off
                for ii in 1:5
                    ssT[ii] = ss.captures[ii]
                end
                if isempty(sss)
                    append!(sss,ssT)
                    sss=reshape(sss, 1, length(sss))
                else
                    sss = [sss; ssT]
                end
            end
        end
        # There are entries in the table without data that need to be removed
        # They read "not currently available", which is why we look for "not"
        jj = findall(x->(x=="not"),sss[:,3])
        rows = collect(1:size(sss, 1))
        deleteat!(rows, jj)

        sss = sss[rows,[1:2;4:5]]
        s = download("https://ssd.jpl.nasa.gov/eph_spans.cgi?id=B") #Satellites
        f = open(s); rd = read(f, String); close(f)
        rm(s)

        #Read in bodynumber, begin time, end time, and file (no need to flag " not "
        #or " to "
        ssT = Array{AbstractString}(undef, 1, 4)
        off = 1; T1 = true
        while T1
            ss=match(r"<td.*?>(\d+)</td>.*?<tt>(.*?) to (.*?)</tt>&nbsp;.*?&nbsp;(.*?)&nbsp;",
            rd[off:end])
            if ss == nothing
                T1 = false
            else
                off = ss.offsets[end]+off
                for ii in 1:4
                    ssT[ii] = ss.captures[ii]
                end
                sss = [sss; ssT]
            end
        end

        for ii in 1:size(sss,1) # Remove B.C. and A.D.
            (sss[ii,2][1:4]=="B.C.") && (sss[ii,2]="0000"*sss[ii,2][10:end])
            (sss[ii,3][1:4]=="A.D.") && (sss[ii,3]=sss[ii,3][6:end])
            (sss[ii,3][6:8]=="Jan") && (sss[ii,3]=sss[ii,3][1:5]*"01"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Feb") && (sss[ii,3]=sss[ii,3][1:5]*"02"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Mar") && (sss[ii,3]=sss[ii,3][1:5]*"03"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Apr") && (sss[ii,3]=sss[ii,3][1:5]*"04"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="May") && (sss[ii,3]=sss[ii,3][1:5]*"05"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Jun") && (sss[ii,3]=sss[ii,3][1:5]*"06"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Jul") && (sss[ii,3]=sss[ii,3][1:5]*"07"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Aug") && (sss[ii,3]=sss[ii,3][1:5]*"08"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Sep") && (sss[ii,3]=sss[ii,3][1:5]*"09"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Oct") && (sss[ii,3]=sss[ii,3][1:5]*"10"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Nov") && (sss[ii,3]=sss[ii,3][1:5]*"11"*sss[ii,3][9:end])
            (sss[ii,3][6:8]=="Dec") && (sss[ii,3]=sss[ii,3][1:5]*"12"*sss[ii,3][9:end])
            (sss[ii,2][6:8]=="Jan") && (sss[ii,2]=sss[ii,2][1:5]*"01"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Feb") && (sss[ii,2]=sss[ii,2][1:5]*"02"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Mar") && (sss[ii,2]=sss[ii,2][1:5]*"03"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Apr") && (sss[ii,2]=sss[ii,2][1:5]*"04"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="May") && (sss[ii,2]=sss[ii,2][1:5]*"05"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Jun") && (sss[ii,2]=sss[ii,2][1:5]*"06"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Jul") && (sss[ii,2]=sss[ii,2][1:5]*"07"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Aug") && (sss[ii,2]=sss[ii,2][1:5]*"08"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Sep") && (sss[ii,2]=sss[ii,2][1:5]*"09"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Oct") && (sss[ii,2]=sss[ii,2][1:5]*"10"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Nov") && (sss[ii,2]=sss[ii,2][1:5]*"11"*sss[ii,2][9:end])
            (sss[ii,2][6:8]=="Dec") && (sss[ii,2]=sss[ii,2][1:5]*"12"*sss[ii,2][9:end])
        end # convert to days from J2000
        #ss=strrep(ss,'--','-Aug-');#fix error on website
        tephf.numbers = sss[:,1]
        tephf.f = sss[:,4]
        dateformat = DateFormat("y-m-d")
        tephf.tl =
        [Dates.datetime2julian.(DateTime.(sss[:,2], dateformat)) Dates.datetime2julian.(DateTime.(sss[:,3], dateformat))] .- 2451545.

        s = download("https://ssd.jpl.nasa.gov/eph_spans.cgi?id=D") #Satellites
        f = open(s); rd = read(f, String); close(f)
        rm(s)

        ssT = Array{AbstractString}(undef, 1, 2)
        off = 1; T1 = true; ss1=Array{AbstractString}(undef, 0)
        while T1
            ss=match(r"&nbsp;(\S+) to (\S+)&nbsp;",rd[off:end])
            if ss == nothing
                T1 = false
            else
                off = ss.offsets[end]+off
                ssT[1] = ss.captures[1]
                ssT[2] = ss.captures[2]
                for ii in 1:2
                    (ssT[ii][6:8]=="Jan") && (ssT[ii]=ssT[ii][1:5]*"01"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Feb") && (ssT[ii]=ssT[ii][1:5]*"02"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Mar") && (ssT[ii]=ssT[ii][1:5]*"03"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Apr") && (ssT[ii]=ssT[ii][1:5]*"04"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="May") && (ssT[ii]=ssT[ii][1:5]*"05"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Jun") && (ssT[ii]=ssT[ii][1:5]*"06"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Jul") && (ssT[ii]=ssT[ii][1:5]*"07"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Aug") && (ssT[ii]=ssT[ii][1:5]*"08"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Sep") && (ssT[ii]=ssT[ii][1:5]*"09"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Oct") && (ssT[ii]=ssT[ii][1:5]*"10"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Nov") && (ssT[ii]=ssT[ii][1:5]*"11"*ssT[ii][9:end])
                    (ssT[ii][6:8]=="Dec") && (ssT[ii]=ssT[ii][1:5]*"12"*ssT[ii][9:end])
                end
                if isempty(ss1)
                    append!(ss1,ssT)
                    ss1=reshape(ss1, 1, length(ss1))
                else
                    ss1 = [ss1; ssT]
                end
            end
        end
        tephf.sb =
        [Dates.datetime2julian.(DateTime.(ss1[:,1], dateformat)) Dates.datetime2julian.(DateTime.(ss1[:,2], dateformat))] .- 2451545.
        param(dict,"tephf",tephf)
    end
    return tephf
end
