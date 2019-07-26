"""
	b = getspk(body, dict)

Tries to match the given string to an SPK number and returns the body's number
if found.
For Lagrange points, format string as "Body1-Body2 L(1-5)" where Body1 is the
dominant body. Lagrange points may not exist in SPICE system.
"""
function getspk(bodies::AbstractArray{AbstractString,1}, dict::Dict{String,Any})::AbstractArray{Int64, 1}
    b = Array{Float64}(undef, length(bodies)) # Preallocate
    for k in 1:length(bodies)
        body = bodies[k]
        b[k] = getspk(body, dict)
    end
    return b
end

function getspk(body::AbstractString,dict::Dict{String,Any})::Int64
    if(isempty(dict))
        tryloaddict!(dict)
    end

    mb = param(dict, "mb")
    sb = param(dict, "sb")

    #get mb from saved data unless empty or ssd flagged
    if isempty(mb.Names) || param(dict,"ssd")
        mb = getmb(dict)
    end

    L = match(r"L\d", body) #Lagrange point

    if L == nothing # Looking for a body
        x = mb
        index = uniquestr(mb.Names,body,1) #check whole word in major body data

        if isempty(index)
            if !isempty(fieldnames(typeof(sb)))
                x = sb
                index = uniquestr(sb.Names,body,1) #check whole word in small body data
            else
                index = []
            end

            if (isempty(index))||(param(dict,"ssd"))
                id = uniquesb(body,dict)
                if id != NaN
                    return id[1]
                end
            end#check ssd webpages

            if isempty(index)
                @warn "No whole-word match found for $body"
                x = mb
                index = uniquestr(mb.Names,body,2) #check fragment in mb
                if isempty(index)
                    x = sb
                    index = uniquestr(sb.Names,body,2) #check fragment in sb
                end#sb any
                if !isempty(index)
                    @warn "Parsing $body as $(x.Names[index,:])"
                end
            end#sb word
        end#mb any

        if isempty(index)
            @warn "No match found for $body"
            return NaN
        else
            spkid = x.numbers[index,:][1]
        end
    else # Looking for a Lagrange point
        L = parse(Int64, string(L.match[end]))

        bodiesStrings = split(replace(body, r" L\d" => ""), "-")
        if(length(bodiesStrings) != 2)
            @error "Malformed Lagrange point string"
            return NaN
        end
        bodies = Array{Int}(undef, 2)
        bodies[1] = getspk(bodiesStrings[1], dict)
        bodies[2] = getspk(bodiesStrings[2], dict)
        if any(isnan.(bodies)) || L > 5
            @error "$body is an ambiguous Lagrange point"
            return NaN
        end

        if bodies[1] != getcb(bodies[2])
            @error "$body is an ambiguous Lagrange point"
            return NaN
        end

        # Langrange points are formatted like 3991 for S-E-M L1, 3992, 3993, ect.
        # However, there isn't data for many of these...
        spkid = bodies[1]*10 + L
        if any(spkid == mb.numbers)
            putsbmb(spkid,"numbers", spkid, dict)
            putsbmb(spkid, "names", body, dict)
        end
    end

    return spkid
end

"""
	bo = getnumout(body, dict)

Outputs the associated number of the specified body input.
"""
function getnumout(body::T,dict::Dict{String,Any}) where T
    if isnan(body)
        return NaN
    else
        n = getsymbol(body,"numbers",dict) #check saved data (& no ssd flag)
        if !isempty(n)#found it
            n = n[1]
        elseif body > 1e3 && body < 1e4 && mod(body,10) < 6
            n = body #Lagrange point
        elseif body < 1e6
            x = getmb(dict)
            n = findfirst(y-> y == body, x.numbers) #check list of major bodies
            (!iszero(n)) ? (n = x.numbers[n]) : (n = NaN)
        else
            n,_ = getsb(body,dict)
            n=n[1]#check ssd
        end
    end#isnan
    if isempty(n)
        @warn "No match found for $body"
        return NaN
    end
    return n
end

"""
	n = getname(body, dict)

Gets the associated name of the specified body input.
"""
function getname(body::Integer,dict::Dict{String,Any})::AbstractString
    n = getsymbol(body,"Names",dict) #check saved data (& no ssd flag)
    if !isempty(n)#found it
        return n;
    elseif body > 1e3 && body < 1e4 && mod(body,10) < 6
        L = mod(body,10) #Lagrange point
        body = floor(Int,body/10)
        T1 = getname(getcb(body),dict)
        T2 = getname(body,dict)
        n = [T1[1]*"-"*T2[1]*" "*string(L)]
        putsbmb(10*body+L, "Names", n, dict)
        putsbmb(10*body+L, "numbers", 10*body+L, dict)
    elseif body<1e3
        x=getmb(dict)
        n=findall(y->(y==body),x.numbers)#check list of major bodies
        if !isempty(n)
            n=x.Names[n][1]
        end
    else
        n,_ = getsb(body,dict)  #check ssd
        if !isnan(n[1])
            n=n[2]
        else
            n=String[]
        end
    end

    if isempty(n)
        n="NULL"
        @warn "No match found for $body"
    end#nomatch
    #format small body name, tokens in regexp are {number, name, (YYYY A1)}
    if body > 1e6
        n1 = something(findfirst("(",n), 0:-1)
        if n1 != 0:-1
            if n1[1] == 1
                nT2 = String[]
            else
                nT2 = n[1:n1[end]-1]
            end
            nT3 = n[n1[end]:end]
        else
            nT2 = String[]
            nT3 = String[]
        end
        n1=something(findfirst(r" ",n), 0:-1)
        if n1 != 0:-1
            n2 = n[n1[end]+1:end]
        end

        #if named, go with name, if numbered go with "# (provisional)",
        # otherwise output provisional designation
        if !isempty(nT2)
            n=nT2
        elseif isempty(nT3)
            n=n2
        elseif !isempty(nT3)
            n=nT3[2:end-1]
        end
    end
    return n
end

"""
	gm = getgm(body, dict)

Finds the inputted body's gravitational constant. Unless it is zero, it will estimate
based on density and size of body
"""
function getgm(body::Integer, dict::Dict{String,Any})::Float64
    # GM from ephemeris header constants unless = 0, then estimate from
    #   density and size
    gm = getsymbol(body,"gm",dict) #check saved data (& no ssd flag)

    if (!isempty(gm)) #found it

    elseif body > 1e3 && body < 1e4 && mod(body,10) < 6
        gm = NaN
        putsbmb(body, "gm", gm, dict) #Lagrange point
    elseif body < 1e6
        ed = getephdata(body, dict) #major body, check ephemeris header for GM
        if (!isempty(ed)) && (ed["gm"*string(floor(Int64,body))][1]>1e-9)
            gm = ed["gm"*string(floor(Int64,body))][1]
            putsbmb(body,"gm",gm,dict)
            #read general satellite page, get first # after >name < and ">"
        else
            ed = getsatdata(dict)
            n = getname(body,dict)
            T1 = Regex(join([">",n[1],"[\\s<]"],""))
            r1 = something(findfirst(T1,ed), 0:-1)
            if (r1 != 0:-1)
                r2 = something(findfirst(r".([\d.])+[\s&<]",ed[r1[1]:end]), 0:-1)
                gm = parse(Float64,ed[r2[1]+r1[1]:r2[end]+r1[1]-2])
                putsbmb(body,"gm",gm,dict)
            else
                @warn "No GM found for body $body"
                return NaN
            end ##ok<*WNTAG>
        end
    else
        gm,_ = getsb(body,dict)#small body
        if !isnan(gm[1])
            gm = gm[3]
        else
            @warn "No GM found for body $body"
            return NaN
        end
    end

    return gm
end

"""
	j2 = getj2(body, dict)

Gets the Julian date for specific time.
"""
function getj2(body::Integer,dict::Dict{String,Any})::Float64
    if body in [10; 399; 499; 599; 699; 799; 899; 999; 301]
    #J2 currently available for only a few bodies from ephemeris header
        j2 = getsymbol(body,"j2",dict)
        if isempty(j2)
            ed = getephdata(body,dict)
            j2 = ed["j2_"*string(body)]
            putsbmb(body,"j2",j2,dict)
        end
        return j2[1]
    else
        @warn "J2 is only available for the following bodies: $([10; 399; 499; 599; 699; 799; 899; 999; 301])"
        return NaN
    end
end

"""
	rx = getradius(body, dict)

Gets radius for inputted body or finds body with associated radius.
"""
function getradius(body::Integer,dict::Dict{String,Any})::Float64
    r = getsymbol(body,"rad",dict)
    if !isempty(r)#found it
    elseif body < 10 # System barycenters
        @warn "Cannot retrieve radius for barycenter $body"
        return NaN
    elseif body > 1e3 && body < 1e4 && mod(body,10)<6 #Lagrange point
        @warn "Cannot retrieve radius for lagrange point $body"
        putsbmb(body,"rad",r,dict)
        return NaN
    elseif (body < 1e6) && ((mod(body,100) == 99)||(body == 10))
        #planet or sun, ephemeris header
        ed = getephdata(floor(Int64,body),dict)
        r = ed["rad"*string(floor(Int64,body))]
        if isempty(r)
            r = gettriaxial(body,dict)
            r = r[1]
        else
            r = r[1]
        end
        putsbmb(body,"rad",r,dict)
    elseif body < 1e6
    #satellite, read general satellite page, get second # after >name < and ">"
    ed = getsatdata(dict)
    n = getname(body,dict)
    T1 = something(findfirst("eft>"*n[1],ed), 0:-1)
    if T1 != 0:-1
        T2 = something(findfirst(r"[\s<].*?>[\d.]",ed[T1[end]:end]), 0:-1) .+ T1[end]
        T3 = something(findfirst(r"[\s&<].*?>([\d.]+)",ed[T2[end]:end]), 0:-1) .+ T2[end]
        r = ed[T3[1]:T3[end]]
        T4 = something(findfirst(r"\d",r), 0:-1)
        r = r[T4[1]:end]
        (r[end]=='&') && (r = r[1:end-1])
        (r[end]=='<') && (r = r[1:end-1])
        putsbmb(body,"rad",r,dict)
        r = parse(Float64,r)
    else
        @warn "No radius found for body $body"
        return NaN
    end
    else #small body
        r,_ = getsb(body, dict)
        if !isnan(r[1])
            r=r[4]
        else
            @warn "No radius found for body $body"
            return NaN
        end
    end

    return r
end

"""
	tri = gettriaxial(body, dict)

Finds the triaxial radii for the given body.
"""
function gettriaxial(body::Integer,dict::Dict{String,Any})::AbstractArray{Float64, 1}
    r = getsymbol(body,"tri",dict)
    if isempty(r) || isnan(r[1]) # Check if we haven't retrieved any real data
        r = getpck(body,"RADII",dict) #check pck file
        if (isempty(r)) && (body > 1e6)#small body
            #extent is flag for triaxial radii,convert to #,
            #   make axisymmetric if only 2#s
            nn,ee = getsb(body,dict)
            s = download(string("https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=",string(body)))
            f = open(s); rd = read(f, String); close(f)
            rm(s)
            r1 = something(findfirst(r">extent<.*?>",rd), 0:-1)
            if r1 != 0:-1
                off = r1[end]
                r1 = something(findfirst(r">extent<.*?>",rd[off:end]), 0:-1)
                r2 = something(findfirst(r">\d",rd[off+r1[1]:end]), 0:-1)
                r3 = something(findfirst("<",rd[off+r1[1]+r2[1]:end]), 0:-1)
                rs = rd[off+r1[1]+r2[1]:off+r1[1]+r2[1]+r3[1]-2]
                T1 = split(rs,'x')
                r = zeros(size(T1))
                for ii in 1:length(r)
                    r[ii] = Meta.parse(T1[ii])
                end
                (length(r) == 2) && (r = [r[1] r])
            end
        end
        if isempty(r) #see if there's a single radius value
            T1 = getradius(body, dict)
            r = T1 .* [1.0; 1.0; 1.0]
        end
        (!isnan(r[1])) && (putsbmb(body,"tri",[r],dict))
    end#if

    return r
end

"""
	r = getrotperiod(body, dict)

Retrieves the sidereal rotation period of the requested body in solar days.
"""
function getrotperiod(body::Integer,dict::Dict{String,Any})::Float64
    r = getsymbol(body,"rot",dict)
    if isempty(r)
        r = getpck(body,"pm",dict) #check pck file
        if !isempty(r)
            r = 360/r[2]
            putsbmb(body,"rot",r,dict) #convert to days
        elseif body > 1e6 #small body
            r,_ = getsb(body,dict)
            if !isnan(r[1])
                r = r[5]
            else
                @warn "No rotation data found for body $body"
                return NaN
            end
        else
            @warn "No rotation data found for body $body"
            return NaN
        end
    end#if

    return r
end

"""
	sd = getsatdata(dict)

Retrieves satellite data (such as gm, mean radius, density, magnitude, and albedo)
from JPL's Solar System Dynamics database.
"""
function getsatdata(dict::Dict{String,Any})
    #Satellite data page, has GM and mean radius (& density, magnitude, albedo)
    sd=param(dict,"satdat")
    if isempty(sd)
        s=download("https://ssd.jpl.nasa.gov/?sat_phys_par")
        f = open(s); sd = read(f, String); close(f)
        rm(s)

        T1 = something(findfirst("Earth's Moon",sd), 0:-1)
        if (T1 != 0:-1)
            sd = sd[T1[1]:end]
            T2 = something(findfirst("References",sd), 0:-1)
            (T2 != 0:-1) && (sd = sd[1:T2[1]])
        end
        param(dict,"satdat",sd)
    end
    return sd
end

"""
	d = getpck(body, variable, dict)

Retrieves the requested variable for the given body from NAIF pack files
"""
function getpck(body::Int64,variable::String,dict::Dict{String,Any})
    #NAIF pck file, analytic orientation & contains some small body radii and
    #orientation data not on Horizons
    pck=param(dict,"pck")
    if isempty(pck)

        client = FTP("ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/")
        dir = readdir(client)
        close(client)
        pck = match(r"pck\d+.tpc",join(dir)).match

        #get latest file in directory
        s=download("ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/"*pck)
        #read contents to string
        f = open(s); rd = read(f, String); close(f)
        rm(s)

        rd = replace(rd,"\n" => ";"); rd = replace(rd,"(" => "["); rd = replace(rd,")" => "]")
        #turn newlines into ; and parens into brackets for matlabese
        rd = replace(rd,"2431010" => "2000243")
        pck = replace(rd,"9511010" => "2000951") #change # of Ida and Gaspra to SPK-ID
        param(dict,"pck",pck)
    end#save
    v1 ="BODY"*string(floor(Int64,body))*"_"*uppercase(variable) #variable name
    T1 = Regex(join([v1,".*?]"],""))
    if match(T1,pck) == nothing
        m1 = []
    else
        m1 = match(T1,pck).match#just single param
    end
    #d=regexp(pck,['\\begindata[;\s]+(' v{1} '.*?;)[;\s]+\\begintext'],'tokens');
    #everything in data block
    d = Float64[]
    if !isempty(m1)
        m1 = split(m1)
        for ii in 1:length(m1)
            (m1[ii]=="[") && (continue)
            (m1[ii]=="]") && (continue)
            (m1[ii]=="=") && (continue)
            if tryparse(Float64, m1[ii]) != nothing
                d = [d; parse(Float64,m1[ii])]
            else
                m1[ii] = replace(m1[ii],'D' => 'e')
                if tryparse(Float64, m1[ii]) != nothing
                    d = [d; parse(Float64,m1[ii])]
                end
            end
        end
    end
    return d
end

"""
	x = getmb(dict)

Retrieves major body list from Horizons.
"""
function getmb(dict::Dict{String,Any})
    #major body list
    nn=download("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=MB")
    f = open(nn); lines = readlines(f); close(f)
    rm(nn)

    NumNam = []
    for ii in 1:length(lines) #keep 1--3 numbers, skip spaces, keep stuff until two spaces
        split1 = split(lines[ii], "  "; keepempty=false) # skip spaces
        if isempty(split1)
            continue
        end
        if length(split1) > 1 # remove extra spaces
            split1[1] = strip(split1[1])
            split1[2] = strip(split1[2])
        end
        if all(isnumeric, split1[1])
            nums = parse(Int64,split1[1]); nams = split1[2]
            if (nums < 1e4) && nums >= 0 # keep 1--3 numbers and stuff until 2 spaces
                if isempty(NumNam)
                    NumNam = [nums nams]
                else
                    NumNam = vcat(NumNam, [nums nams]) #numbers and names
                end
            end
        end
    end

    jj = findall(x->(x>9),NumNam[:,1]) #match barycenters last
    nums = zeros(size(NumNam,1))
    nams = Array{AbstractString}(undef, size(NumNam,1))
    ct = 0
    for ii in jj
        nums[ct+1] = NumNam[ii,1]
        nams[ct+1] = NumNam[ii,2]
        ct = ct + 1
    end
    jj = findall(x->(x<=9),NumNam[:,1])
    for ii in jj
        nums[ii+ct] = NumNam[ii,1]
        nams[ii+ct] = NumNam[ii,2]
    end

    #b = vcat(nums[1:10],nums[15:20],nums[22:23],nums[24:50],nums[160:end])
    x,_=getsbmb("mb",dict)
    b=x.numbers
    jj = [] # see if body numbers already exist in data structure
    for ii in 1:length(b)
        T1 = findall(z->(z==b[ii]),nums)
        if isempty(jj)
            jj = T1
        else
            jj = vcat(jj,T1)
        end
    end

    T1 = Array{AbstractString}(undef, length(jj))
    ct = 0
    for ii in 1:length(jj) #put new bodies last
        splice!(nums,jj[ii]-ct)
        T1[ii] = nams[jj[ii]-ct]
        splice!(nams,jj[ii]-ct)
        ct = ct + 1
    end
    nums = vcat(b,nums)
    nams = vcat(T1,nams)

    x.numbers=nums;x.Names=nams
    param(dict,"mb",x) #save data
    return x
end

"""
	nnf, eph = getsb(n, dcit)

Retrieves small body list from Horizons.
"""
function getsb(n,dict::Dict{String,Any})
    #small body web pages
    if typeof(n) == String #name or number
        srch = escape_string(n)*"*" # Add wildcard to help match incomplete names
    else
        srch=string(convert(Int64,n))
    end
    s = download(string("https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=",srch))
    f = open(s); rd = read(f, String); close(f)
    rm(s)

    ss = something(findfirst(r"\+1\"><b>[^<]+", rd), 0:-1) #Name is bigger font "+1"
    if ss != 0:-1
        nam = rd[ss[8]:ss[end]]
        ss1 = something(findfirst(r">\d{7}<", rd), 0:-1)
        nums = rd[ss1[2]:ss1[end-1]]
        nums = parse(Int64, nums)
        if typeof(n) != String && n != nums
            println("SPK ID of ", n," has been changed to ", nums)
            n = [n nums]
        else
            n = nums
        end
        nnf = [nums nam] # [number name gm rad rot h occ type]
        #read in data off table
        ss1 = something(findfirst(r">GM<.*?>",rd), 0:-1)
        if ss1 != 0:-1
            ss2 = something(findfirst(r">([\d.e-]{1,})</",rd[ss1[end]:end]), 0:-1)
            a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
            (isempty(a)) && (a=0); (isempty(a)) || (a=parse(Float64,a))
        else
            a = 0
        end
        nnf = [nnf a]
        ss1 = something(findfirst(r">diameter<.*?>",rd), 0:-1)
        if ss1 != 0:-1
            ss2 = something(findfirst(r">([\d.]{1,})</",rd[ss1[end]:end]), 0:-1)
            a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
            (isempty(a)) && (a=0); (isempty(a)) || (a=parse(Float64,a)/2.)
        else
            a = 0
        end
        nnf = [nnf a]
        ss1 = something(findfirst(r">rot_per<.*?>",rd), 0:-1)
        if ss1 != 0:-1
            ss2 = something(findfirst(r">([\d.]{1,})</",rd[ss1[end]:end]), 0:-1)
            a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
            (isempty(a)) && (a=0); (isempty(a)) || (a=parse(Float64,a)/24.)
        else
            a = 0
        end
        nnf = [nnf a]
        ss1 = something(findfirst(r">H<.*?>",rd), 0:-1)
        if ss1 != 0:-1
            ss2 = something(findfirst(r">([\d.]{1,})</",rd[ss1[end]:end]), 0:-1)
            a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
            (isempty(a)) && (a=NaN); (isempty(a)) || (a=parse(Float64,a))
        else
            a=NaN
        end
        nnf = [nnf a]
        ss1 = something(findfirst(r">condition code<.*?sp;(\d)&nb",rd), 0:-1)
        if ss1 != 0:-1
            ct = 0
            for ii in 1:ss1[1]
                ct = ii-1
                (rd[ss1[end]-ii] == ';') && (break)
            end
            a = parse(Float64,rd[ss1[end]-ct:ss1[end]-3])
        else
            a=NaN
        end
        nnf = [nnf a]
        ss1 = something(findfirst(r">spec_T<.*?>",rd), 0:-1)
        if ss1 != 0:-1
            ss2 = something(findfirst(r">([\w])</",rd[ss1[end]:end]), 0:-1)
            if ss2!= 0:-1 #dev
                a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
            else
                a = []
            end
            (isempty(a)) && (a="")
        else
            a=""
        end
        ss1 = something(findfirst(r">spec_B<.*?>",rd), 0:-1)
        if ss1 != 0:-1
            ss2 = something(findfirst(r">(\w+)</",rd[ss1[end]:end]), 0:-1)
            (ss1 != 0:-1) ? (b = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]) : (b="_")
            (!isempty(a)) && (b="/"*b)
        else
            b="_"
        end
        nnf = [nnf a*b]
        ss1 = something(findfirst(r">Reference: <.*?>.*?>.*?>(.*?)<",rd), 0:-1)
        if ss1 != 0:-1
            ct = 0
            for ii in 1:ss1[1]
                ct = ii-1
                (rd[ss1[end]-ii] == '>') && (break)
            end
            eph = rd[ss1[end-ct]:ss1[end]-1] #Ephemeris ID
        end
        #save data
        putsbmb(n,["numbers" "Names" "gm" "rad" "rot" "h" "occ" "types"],nnf,dict)
    else
        #see if returned multiple matches
        @warn "multiple matches might be possible in boddat. possibly not coded correctly"
        ss2 = something(findfirst(r">([^<]+)</a></td>",rd), 0:-1)
        nnf = [NaN]; eph = ""
    end
    return nnf, eph
end

"""
	nn = uniquesb(b, dict)

Searches for the unique small body name.
"""
function uniquesb(bodyName::AbstractString, dict::Dict{String,Any})
    bs=[bodyName,string(bodyName,"*"),string("*",bodyName,"*")] #search exact, beginning, fragment
    nn = Float64[] #empty array that will be a float
    for b1 in bs
        nn,_ = getsb(b1,dict) #check ssd
        if (nn==NaN) && !isempty(nn)
            @warn "$bodyName returned multiple matches."
        end#if
        if !isnan(nn[1])
            return nn
        end#found one!
    end#for
    return nn
end

"""
	In = uniquestr(xni, bi, nb)

Looks for a unique match for a body input.

Example: if there are 5 asteroids that look for 'ces', it looks for a unique string
after this to differentiate between them.
"""
function uniquestr(xni,bi,nb::Int64=0)
    #find a unique match of bi in xni
    #1 for only whole word matches, 2 for only fragment matches, 0 for either
    if typeof(xni) == String
        T1 = false
        (occursin(bi,xni)) && (T1=true)
    else
        T1 = falses(length(xni))
        for ii in 1:length(xni)
            (occursin(bi,xni[ii])) && (T1[ii]=true)
        end
    end
    In = findall(T1)
    #in=strfind(xn,[' ' bi ' ']);#whole word

    if (length(In) > 1) && (nb == 1) #multiple entries with same word
      #In is strfind index, xm tracks current list of multiple matches
        xm = Array{AbstractString}(undef, length(In))
        for ii in 1:length(In)
            xm[ii] = xni[In[ii]]
        end
        In = In[findall(x->(x==bi),xm)] #exact match
        if length(In) > 1
            @warn "identical entries"
            In = In[1]
        end
        (isempty(In)) && (@warn "ambiguous string match")
    elseif (length(In) > 1)  && (nb!=1)
        if any(xni[In] .== bi)
        elseif any(x->(x==xni[In[1]]),xni[In[2:end]])
            @warn "ambiguous string match"
        elseif any(x->(x[1:3]==xni[In[1]][1:3]),xni[In[2:end]])
            @warn "ambiguous string match"
        end
        In=In[1]
    end
    return In
end

"""
	cb = getcb(body)

Retrieves the central body for the given body
"""
function getcb(body::Int64)::Int64
    #Calculates the central body
    cb = 10
    if body > 10 && body < 1000 && mod(body,100) != 99
        cb = floor(Int64,body/100)*100+99
    elseif body > 1e3 && body < 1e4
        cb = floor(Int64,body/10)
    end
    return cb
end
