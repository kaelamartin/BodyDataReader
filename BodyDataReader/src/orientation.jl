"""
	poleArr = pole(body, timeRange, dict)

Retrieves the pole unit vector of the inputted body at each given time.
"""
function pole(body::Int64, timeRange::AbstractArray{Float64},dict::Dict{String,Any})::Array{Float64,2}
    poleArr = Array{Float64}(undef, 3, length(timeRange)) # Preallocate
    for k in 1:length(timeRange)
        time = timeRange[k]
        poleArr[:,k] = pole(body, time, dict)
        if isnan(poleArr[1,k]) # In case of undefined orientations
            return NaN*ones(3,1)
        end
    end
    return poleArr
end

"""
	poleVec = pole(body, time, dict)

Retrieves the pole unit vector of the inputted body at the given time.
"""
function pole(body::Int64,time::Float64,dict::Dict{String,Any})::Array{Float64,1}
    #small body pole not in Horizons, may be in pck file.
    if body > 1e3
        return qpole(body, dict)
    end

    #get saved spline data
    d=param(dict,"orient_data")
    if(isempty(d))
        d = Dict()
    end

    ssd=param(dict,"ssd")
    if !haskey(d, body) || !haskey(d[body], "pole_data") || ssd#read from file or horizons
        if ssd #te is saved time span
            te = [time]
            type = single
        else
            te = [-0.5; 18261.5; 0.1]
            type = multiple
        end

        err = 0
        @static Sys.iswindows() ?
        (edir = string(param(dict,"bdir"),"\\ephem\\pol",body,".jld2")) :
        (edir = string(param(dict,"bdir"),"/ephem/pol",body,".jld2"))
        rf = (!ssd && isfile(edir))
        if rf
            #read from file if rf, 1st two entries is data size
            di = load(edir,"di")
        else
            X,td,err = makeephem(body,te,type,dict,[0, 90.])#pole
            sp = falses(size(X,2))

            for kk in 1:size(X,2)
                if isnan(X[1, kk])
                    @warn "No pole data available for $body"
                    return NaN*ones(3)
                end
                (X[3,kk]<0) && (sp[kk]=true)
            end

            if any(sp) #z component is positive for sqrt
                X[:,sp]=-X[:,sp]
                if any(.!(sp))
                    @warn "pole crosses ecliptic"
                end
            end

            # Truncate X and normalize
            X = X[1:3,:]
            for kk in 1:size(X,2) #unit vector
                Xmag = sqrt(dot(X[1:3,kk],X[1:3,kk]))
                X[1:3,kk]=X[1:3,kk]/Xmag
            end

            di=[td;X[1:2,:]]
        end#rf

        if any(isnan, di[2,:]) #data to save, don't save nans
            jj = findall(isnan.(di[2,:]))
            for kk in 1:length(jj)
                di = [di[:,1:jj[kk]-1] di[:,1:jj[kk]+1]]
                jj = jj .- 1
            end
        end

        if isempty(di)
            return NaN*ones(3)#no data, output nans
        end

        #see if times are out of horizons range if ssd
        if ssd
            (length(err)==1) && (err = [err false])
            if err[1]
                @warn "Min Horizons ephem time for $body is EXCEEDED"
            elseif err[2]
                @warn "Max Horizons ephem time for $body is EXCEEDED"
            end
            #save data if no error and didn't read from existing file
        elseif (!rf) && (param(dict,"save")) && all(x -> x == 0, err)
            jldopen(edir,"w") do fid
                write(fid, "di", di)
            end
        end

        if !ssd#make spline if not ssd flag

            #see if body already exists in orient data
            if !haskey(d,body)
                d[body] = Dict()
            end

            d[body]["pole_data"] = Dict{String,Array{Float64}}()
            spline!(di[1,:], di[2:end,:], d[body]["pole_data"])
            if param(dict,"keeps")
                param(dict,"orient_data", d)
            end#save data
        end
    end#ssd

    ist=false
    if !ssd
        #see if data is out of bounds
        tl=d[body]["pole_data"]["breaks"][1]
        if time < tl
            ist = true
            @warn "Min time for $body ephem spline is $(Dates.julian2datetime(tl+2451545))"
        end
        tl=d[body]["pole_data"]["breaks"][end]
        if time > tl
            ist = true
            @warn "Max time for $body ephem spline is $(Dates.julian2datetime(tl+2451545))"
        end
    end#ssd

    if !ist #Only execute on the spline
        T1 = zeros(2,1)
        pvalorient!(d[body]["pole_data"],[time],4,T1)
        #get new pole data if not already in spline
        X = Array{Float64}(undef, 3) #pole vector
        X[1:2,:] = T1
        p12 = 0
        p12 = X[1]^2 + X[2]^2
        X[3] = sqrt(1-p12)
    else #get horizons data for out of range times
        sav=param(dict,"save")
        param(dict,"save",false); param(dict,"ssd",true)#reset flags
        X = pole(body, time, dict)
        param(dict,"save",sav)
        param(dict,"ssd",ssd)
    end

    return X
end

"""
	poledcmVec = poledcm(body, time, dict)

Retrieves the pole dcm vector of the inputted body at the given time.
"""
function poledcm(body::Int64, time::Float64, dict::Dict{String,Any})::Array{Float64, 1}
    #DCM is [Node;Q;Pole] where Node=cross([0;0;1],Pole), Q=cross(Pole,Node)
    Pole = pole(body, time, dict)
    Node = cross([0;0;1],Pole)
    Q = cross(Pole, Node)
    return [Node;Q;Pole]
end

"""
	pmArr = pm(body, timeRange, dict)

Retrieves the pm unit vector of the inputted body at each given time.
"""
function pm(body::Int64, timeRange::AbstractArray{Float64},dict::Dict{String,Any})::Array{Float64,2}
    pmArr = Array{Float64}(undef, 3, length(timeRange)) # Preallocate
    for k in 1:length(timeRange)
        time = timeRange[k]
        pmArr[:,k] = pm(body, time, dict)
        if isnan(pmArr[1,k])
            return NaN*ones(3,1)
        end
    end
    return pmArr
end

"""
	pmVec = pm(body, time, dict)

Retrieves the prime meridian of the inputted body at the given time.
"""
function pm(body::Int64,time::Float64,dict::Dict{String,Any})::Array{Float64,1}
    #small body pole not in Horizons, may be in pck file. write xx 1-to-1 if tb
    if body > 1e3
        return qpm(body,time,dict)
    end

    #get saved spline data
    d=param(dict,"orient_data")
    if(isempty(d))
        d = Dict()
    end

    ssd=param(dict,"ssd")
    if !haskey(d, body) || !haskey(d[body], "pm_data") || ssd#read from file or horizons
        if ssd #te is saved time span
            te = [time]
            type = single
        else
            te = [-0.5; 18261.5; 0.1]
            type = multiple
        end

        err = 0
        @static Sys.iswindows() ?
        (edir = string(param(dict,"bdir"),"\\ephem\\pm",body,".jld2")) :
        (edir = string(param(dict,"bdir"),"/ephem/pm",body,".jld2"))
        rf = (!ssd && isfile(edir))
        if rf
            #read from file if rf, 1st two entries is data size
            di = load(edir,"di")
        else
            X,td,err = makeephem(body,te,type,dict,[0, 0.])
            if type == multiple
                P = pole(body, range(te[1], te[2]; step=te[3]), dict)
            else
                P = pole(body, time, dict)
            end

            if isnan(P[1,1]) # Not guaranteed protection, but doesn't impact performance much
                @warn "No prime meridian data available for $body"
                return NaN*ones(3)
            end

            X = X[1:3,:]
            for kk in 1:size(X,2) #unit vector
                if isnan(X[1, kk])
                    @warn "No prime meridian data available for $body"
                    return NaN*ones(3)
                end
                Xmag = sqrt(dot(X[1:3,kk],X[1:3,kk]))
                X[1:3,kk]=X[1:3,kk]/Xmag
            end

            N = Array{Float64}(undef, 3, size(P,2))
            Q = Array{Float64}(undef, 3, size(P,2))
            w = Array{Float64}(undef, size(P,2))
            for jj in 1:size(P,2)
                p12=sqrt(P[1,jj]^2+P[2,jj]^2)
                N[:,jj]=[P[2,jj]/p12; -P[1,jj]/p12; 0]
                Q[:,jj]=[-P[3,jj]*N[2,jj]; P[3,jj]*N[1,jj]; -p12]#Pole x [0;0;1]
                w[jj] = atan(dot(X[:,jj],Q[:,jj]),dot(X[:,jj],N[:,jj]))
            end
            jj = []
            #clock angle (period has to be longer than 2X sample for unwrap)
            for kk in 1:length(w)-1
                if (abs(w[kk]-w[kk+1])>1.9*pi)
                    (!isempty(jj)) && (jj = [jj; kk])
                    (isempty(jj)) && (jj = kk)
                end
            end
            for kk in jj
                w[kk:end] = w[kk:end] .+ 2*pi
            end
            di=[td;w']
        end#rf

        if any(isnan.(di[2,:])) #data to save, don't save nans
            jj = findall(isnan.(di[2,:]))
            for kk in 1:length(jj)
                di = [di[:,1:jj[kk]-1] di[:,1:jj[kk]+1]]
                jj = jj .- 1
            end
        end

        if isempty(di)
            return NaN*ones(3,1) #no data, output nans
        end

        #see if times are out of horizons range if ssd
        if ssd
            (length(err)==1) && (err = [err false])
            if err[1]
                @warn "Min Horizons ephem time for $body is EXCEEDED"
            elseif err[2]
                @warn "Max Horizons ephem time for $body is EXCEEDED"
            end
            #save data if no error and didn't read from existing file
        elseif (!rf) && (param(dict,"save")) && all(x -> x == 0, err)
            jldopen(edir,"w") do fid
                write(fid, "di", di)
            end
        end

        if !ssd#make spline if not ssd flag

            #see if body already exists in orient data
            if !haskey(d,body)
                d[body] = Dict()
            end

            d[body]["pm_data"] = Dict{String,Array{Float64}}()
            spline!(di[1,:], di[2:end,:], d[body]["pm_data"])
            if param(dict,"keeps")
                param(dict,"orient_data", d)
            end#save data
        end
    end#ii|ssd

    ist = false
    if !ssd
    #see if data is out of bounds
        tl=d[body]["pm_data"]["breaks"][1]
        if time < tl
            ist = true
            @warn "Min time for $body ephem spline is $(Dates.julian2datetime(tl+2451545))"
        end
        tl=d[body]["pm_data"]["breaks"][end]
        if time > tl
            ist = true
            @warn "Max time for $body ephem spline is $(Dates.julian2datetime(tl+2451545))"
        end
    end#ssd

    if !ist #Only execute on the spline
        if (!haskey(d[body], "pole_data")) #get new pole data if not already in spline
            P = pole(body,time,dict)
        else
            T1 = zeros(2,1)
            pvalorient!(d[body]["pole_data"],[time],4,T1)
            P = Array{Float64}(undef, 3); P[1:2] = T1
        end

        p12 = 0
        p12 = P[1]^2 + P[2]^2
        P[3] = sqrt(1-p12)

        #X is pole node: Pole x [0;0;1]
        p12 = sqrt(p12)
        X = Array{Float64}(undef, 3, 1)
        X[1] = P[2]/p12
        X[2] = -P[1]/p12
        X[3] = 0

        w = [0.]
        pvalorient!(d[body]["pm_data"],[time],4,w)
        XT = X[:]
        X[1] = cos(w[1])*XT[1]+sin(w[1])*-P[3]*XT[2]
        X[2] = cos(w[1])*XT[2]+sin(w[1])*P[3]*XT[1]
        X[3] = cos(w[1])*XT[3]+sin(w[1])*-p12[1]
    else #get horizons data for out of range times
        sav=param(dict,"save")
        param(dict,"save",false); param(dict,"ssd",true)#reset flags
        X=pm(body,time,dict)
        param(dict,"save",sav)
        param(dict,"ssd",ssd)
    end

    return X[:]
end

"""
	pmdcmVec = pmdcm(body, time, dict)

Retrieves the prime meridian dcm vector of the inputted body at the given time.
"""
function pmdcm(body::Int64, time::Float64, dict::Dict{String,Any})::Array{Float64, 1}
    #DCM is [Node;Q;Pole] where Node=Prime Meridian, Q=cross(Pole,Node)
    Pole = pole(body, time, dict)
    Node = pm(body, time, dict)
    Q = cross(Pole, Node)
    return [Node;Q;Pole]
end

"""
	eqxArr = eqx(body, timeRange, dict)

Retrieves the eqx unit vector of the inputted body at each given time.
"""
function eqx(body::Int64, timeRange::AbstractArray{Float64},dict::Dict{String,Any})::Array{Float64,2}
    eqxArr = Array{Float64}(undef, 3, length(timeRange)) # Preallocate
    for k in 1:length(timeRange)
        time = timeRange[k]
        eqxArr[:,k] = eqx(body, time, dict)
        if isnan(exqArr[1,k])
            return NaN*ones(3,1)
        end
    end
    return eqxArr
end

"""
	eqxVec = eqx(body, time, dict)

Retrieves the equinox of the inputted body at the given time.
"""
function eqx(body::Int64,time::Float64,dict::Dict{String,Any})::Array{Float64,1}
    #small body orientation not in Horizons, may be in pck file. write xx 1-to-1 if tb
    if body > 1e3
        return qeqx(body, dict)
    end

    #get saved spline data
    d=param(dict,"orient_data")
    if(isempty(d))
        d = Dict()
    end

    ssd=param(dict,"ssd")
    if !haskey(d, body) || !haskey(d[body], "eqx_data") || ssd#read from file or horizons
        if ssd #te is saved time span
            te = time
            typ = false
            X_ = Array{Float64}(undef, 6)
        else
            te = range(-0.5, 18261.5; step=0.1)
            typ = true
            X_ = Array{Float64}(undef, 6, length(te))
        end

        err = 0
        @static Sys.iswindows() ?
        (edir = string(param(dict,"bdir"),"\\ephem\\eqx",body,".jld2")) :
        (edir = string(param(dict,"bdir"),"/ephem/eqx",body,".jld2"))
        rf = (!ssd && isfile(edir))
        if rf
            #read from file if rf, 1st two entries is data size
            di = load(edir,"di")
        else
            cachedSave = param(dict, "save")
            cb = getcb(body)
            ephem!([body; cb], te, state, dict, X_)
            P = pole(body,te,dict)

            if isnan(P[1,1]) # Not guaranteed protection, but doesn't impact performance much
                @warn "No equinox data available for $body"
                return NaN*ones(3)
            end

            X = Array{Float64}(undef, 3, size(X_,2)) #equinox
            for kk in 1:size(X_,2)
                X[:,kk] = cross(P[:,kk],cross(X_[1:3,kk],X_[4:6,kk]))
            end

            X = X[1:3,:]
            for kk in 1:size(X,2) #unit vector
                Xmag = sqrt(dot(X[1:3,kk],X[1:3,kk]))
                X[1:3,kk]=X[1:3,kk]/Xmag
            end

            N = Array{Float64}(undef, 3, size(P,2));
            Q = Array{Float64}(undef, 3,size(P,2));
            w = Array{Float64}(undef, size(P,2))
            for jj in 1:size(P,2)
                p12=sqrt(P[1,jj]^2+P[2,jj]^2)
                N[:,jj]=[P[2,jj]/p12; -P[1,jj]/p12; 0]
                Q[:,jj]=[-P[3,jj]*N[2,jj]; P[3,jj]*N[1,jj]; -p12]#Pole x [0;0;1]
                w[jj] = atan(dot(X[:,jj],Q[:,jj]),dot(X[:,jj],N[:,jj]))
            end
            jj = []
            #clock angle (period has to be longer than 2X sample for unwrap)
            for kk in 1:length(w)-1
                if (abs(w[kk]-w[kk+1])>1.9*pi)
                    (!isempty(jj)) && (jj = [jj; kk])
                    (isempty(jj)) && (jj = kk)
                end
            end
            for kk in jj
                w[kk:end] = w[kk:end] .+ 2*pi
            end
            di=[collect(te)';w']

        end#rf

        if any(isnan.(di[2,:])) #data to save, don't save nans
            jj = findall(isnan.(di[2,:]))
            for kk in 1:length(jj)
                di = [di[:,1:jj[kk]-1] di[:,1:jj[kk]+1]]
                jj = jj .- 1
            end
        end

        if isempty(di)
            return NaN*ones(3,1)#no data, output nans
        end

        #see if times are out of horizons range if ssd
        if ssd
            (length(err)==1) && (err = [err false])
            if err[1]
                @warn "Min Horizons ephem time for $body is EXCEEDED"
            elseif err[2]
                @warn "Max Horizons ephem time for $body is EXCEEDED"
            end
            #save data if no error and didn't read from existing file
        elseif (!rf) && (param(dict,"save")) && all(x -> x == 0, err)
            jldopen(edir,"w") do fid
                write(fid, "di", di)
            end
        end

        if !ssd#make spline if not ssd flag

            #see if body already exists in orient data
            if !haskey(d,body)
                d[body] = Dict()
            end

            d[body]["eqx_data"] = Dict{String,Array{Float64}}()
            spline!(di[1,:], di[2:end,:], d[body]["eqx_data"])
            if param(dict,"keeps")
                param(dict,"orient_data", d)
            end#save data
        end
    end#ii|ssd

    ist = false
    if !ssd
        #see if data is out of bounds
        tl=d[body]["eqx_data"]["breaks"][1]
        if time < tl
            ist = true
            @warn "Min time for $body ephem spline is $(Dates.julian2datetime(tl+2451545))"
        end
        tl=d[body]["eqx_data"]["breaks"][end]
        if time > tl
            ist = true
            @warn "Max time for $body ephem spline is $(Dates.julian2datetime(tl+2451545))"
        end
    end#ssd

    if !ist
        if (!haskey(d[body], "pole_data"))
            P = pole(body,time,dict)
        else
            T1 = zeros(2,1)
            pvalorient!(d[body]["pole_data"],[time],4,T1)
            P = Array{Float64}(undef, 3); P[1:2] = T1
        end
        #get new pole data if not already in spline
        p12 = 0
        p12 = P[1]^2 + P[2]^2
        P[3] = sqrt(1-p12)

        p12 = sqrt(p12)
        X = Array{Float64}(undef, 3) #pole vector
        X[1] = P[2]/p12
        X[2] = -P[1]/p12
        X[3] = 0

        w = [0.]
        pvalorient!(d[body]["eqx_data"],[time],4,w)
        XT = X[:]
        X[1] = cos(w[1])*XT[1]+sin(w[1])*-P[3,1]*XT[2]
        X[2]=cos(w[1])*XT[2]+sin(w[1])*P[3,1]*XT[1]
        X[3]=cos(w[1])*XT[3]+sin(w[1])*-p12[1]
    else #get horizons data for out of range times
        sav=param(dict,"save")
        param(dict,"save",false); param(dict,"ssd",true)#reset flags
        X = eqx(body,time,dict)
        param(dict,"save",sav)
        param(dict,"ssd",ssd)
    end

    return X
end

"""
	eqxdcmVec = eqxdcm(body, time, dict)

Retrieves the Equinox dcm vector of the inputted body at the given time.
"""
function eqxdcm(body::Int64, time::Float64, dict::Dict{String,Any})::Array{Float64, 1}
    #DCM is [Node;Q;Pole] where Node=Equinox, Q=cross(Pole,Node)
    Pole = pole(body, time, dict)
    Node = eqx(body, time, dict)
    Q = cross(Pole, Node)
    return [Node;Q;Pole]
end

"""
	poleVec = qpole(body, dict)

Quickly approximates the pole vector of the inputted body.
"""
function qpole(body::Int64,dict::Dict{String,Any})::Array{Float64,1}
    pp=getsymbol(body,"pp",dict)
    if isempty(pp)#check saved data
        pp = populatepp(body, dict)
    end

    return pp[1:3]
end

"""
	poledcmVec = qpoledcm(body, dict)

Quickly approximates the pole dcm vector of the inputted body at the given time.
"""
function qpoledcm(body::Int64, dict::Dict{String,Any})::Array{Float64, 1}
    #DCM is [Node;Q;Pole] where Node=cross([0;0;1],Pole), Q=cross(Pole,Node)
    Pole = qpole(body, dict)
    Node = cross([0;0;1],Pole)
    Q = cross(Pole, Node)
    return [Node;Q;Pole]
end

"""
	pmVec = qpm(body, time, dict)

Quickly approximates the prime meridian of the inputted body at the given time.
"""
function qpm(body::Int64, time::Float64, dict::Dict{String,Any})::Array{Float64,1}
    #quick orientation calculation

    pp=getsymbol(body,"pp",dict)
    if isempty(pp)#check saved data
        pp = populatepp(body, dict)
    end

    X=pp[1:3] #X = pole

    #need X = pole node
    P=X
    p12=sqrt(P[1]^2+P[2]^2)
    X=[P[2]/p12;-P[1]/p12; 0*p12]
    T1 = getrotperiod(body,dict)
    w = pp[4]+2*pi/T1[1]*time
    #rotate along pole
    T1 = X; X = Array{Float64}(undef, 3, 1); T2 = P;
    X = T1*cos(w)+cross(T2,T1)*sin(w)

    return X
end

"""
	pmdcmVec = qpmdcm(body, time, dict)

Quickly approximates the prime meridian dcm vector of the inputted body at the given time.
"""
function qpmdcm(body::Int64, time::Float64, dict::Dict{String,Any})::Array{Float64, 1}
    #DCM is [Node;Q;Pole] where Node=Prime Meridian, Q=cross(Pole,Node)
    Pole = qpole(body, dict)
    Node = qpm(body, time, dict)
    Q = cross(Pole, Node)
    return [Node;Q;Pole]
end

"""
	eqxVec = qeqx(body, dict)

Quickly approximates the equinox of the inputted body.
"""
function qeqx(body::Int64, dict::Dict{String,Any})::Array{Float64,1}
    #quick orientation calculation

    pp=getsymbol(body,"pp",dict)
    if isempty(pp)#check saved data
        pp = populatepp(body, dict)
    end

    X=pp[1:3] #X = pole

    #need X = pole node
    P=X
    p12=sqrt(P[1]^2+P[2]^2)
    X=[P[2]/p12;-P[1]/p12; 0*p12]
    w=pp[5]

    #rotate along pole
    T1 = X; X = Array{Float64}(undef, 3, 1); T2 = P
    X = T1*cos(w)+cross(T2,T1)*sin(w)

    return X
end

"""
	eqxdcmVec = qeqxdcm(body, dict)

Quickly approximates the Equinox dcm vector of the inputted body.
"""
function qeqxdcm(body::Int64, dict::Dict{String,Any})::Array{Float64, 1}
    #DCM is [Node;Q;Pole] where Node=Equinox, Q=cross(Pole,Node)
    Pole = qpole(body, dict)
    Node = qeqx(body, dict)
    Q = cross(Pole, Node)
    return [Node;Q;Pole]
end

"""
	pp = populatepp(body, dict)

Populates the pp field in dict and returns it.
"""
function populatepp(body::Int64,dict::Dict{String,Any})::Array{Float64,1}
    adp = [getpck(body,"POLE_RA",dict)';getpck(body,"POLE_DEC",dict)';getpck(body,"pm",dict)']*pi/180 #get data from pck file
    if length(adp) == 9
        a = adp[1]; d = adp[2]; p = adp[3]; o = 84381.448/3600*pi/180
        #obliquity @J2000: Lieske, J., "Precession Matrix Based on IAU (1976)
        # System of Astronomical Constants"
        #rotate from EME to EMO
        #basis orthogonal to pole, N is pole node: pole x [0;0;1]
        P = [cos(d)*cos(a); sin(d)*sin(o)+cos(d)*cos(o)*sin(a);
        cos(o)*sin(d)-cos(d)*sin(o)*sin(a)] #pole
        N = [P[2];-P[1];0]
        N = N/sqrt(dot(N,N))
        Q = cross(P,N)
        #basis orthogonal to pole, N is pole node: Earth pole x pole

        #prime meridian
        X = [-cos(p)*sin(a)-cos(a)*sin(d)*sin(p);cos(a)*cos(p)-sin(a)*sin(d)*sin(p);
        cos(d)*sin(p)]
        #rotate to EMO # get angle @ J2000
        X = [X[1];cos(o)*X[2]+sin(o)*X[3];cos(o)*X[3]-sin(o)*X[2]]
        p = atan(dot(X,Q),dot(X,N))

        #get equinox
        X = ephem(body, 0., state, dict)
        X = cross(P,cross(X[1:3],X[4:6]))
        X = X/sqrt(dot(X,X))
        pe = atan(dot(X,Q),dot(X,N))

        #pp = [pole vector; pm angle; eqx angle] @ J2000
        pp = [P;p;pe]
        putsbmb(body,"pp",[pp],dict)
    else
        @warn "no analytic orientation data found for $body"
        return NaN*ones(5)
    end

    return pp
end
