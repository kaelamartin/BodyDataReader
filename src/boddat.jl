function boddat(bvars,bi::Array=[],dict::Dict=Dict(),t=[],opti=[])
# If have ssd flag, orient will not work and problems may occur elsewhere
# qorient and orient not fully tested and might throw errors


## Function needs following Julia packages: HDF5, JLD, URIParser
# Note: Install HDF5 before JLD
# Add packages by: Pkg.add("HDF5"); Pkg.add("JLD"); Pkg.add("URIParser")
# At start of file using boddat put:
#          using JLD,HDF5,URIParser
#           include("path/boddat.jld")
# bi should nominally be a row vector. [399, 499] will reslt in Earth's position
#     relative to Mars

#Body Data from Horizons and ssd.jpl.nasa.gov
#
#   varargout=boddat(parameters,bodies,dictionary,times,options)
#
#"parameters" is a cell array or comma, semicolon, or space delimited string
#   that designates which parameters to output in varargout. Case matters.
#The available parameters are:
#   Number (NAIF or SPK ID, e.g. Sun is 10, Earth is 399, Moon is 301,
#            Ceres is 2000001)
#   Name (body name)
#   CB (central body of bodies)
#   GM (gravitational paramter, km^3/s^2)
#   radius (mean radius, km)
#   triaxial (radii of a triaxial ellipsoid shape model, km)
#   J2 (gravitational harmonic due to oblateness)
#   rotation_period (period of rotation, days)
#Small-body-specific parameters
#   mag (absolute magnitude [ H ], magnitude at 1 AU from Sun and observer)
#   occ (orbital condition code 0-9, higher numbers indicate uncertain orbits)
#   type (spectral type, Tholen/SMASII taxonomic classification)
#   close_approach (table of close approach time, body, and AU distance)
#Ephemeris
#   R (position vector, km)
#   V (velocity vector, km/s)
#   X (position and velocity)
#   A (position, velocity, acceleration, & jerk)
#   RTN (rotating frame direction cosine matrix unit[R ; RxVxR ; RxV])
#   reference (orbit solution ID or ephemeris file, used to track changes)
#   date (Creation date of ephemeris datafile in "ephem" directory, days past
#         J2000. The latest generation date from Horizons is not available.)
#Orientation
#   Pole (spin axis of body, unit vector)
#   PM (prime meridian of body, unit vector)
#   Eqx (Equinox direction, Pole x H, unit vector)
#
#   If "dcm" follows Pole, PM, or Eqx then the direction cosine matrix is output
#   as a 9x1 array, with Pole as the Z direction, and PM or Eqx as X.
#   The X direction for the Pole dcm is Earth pole x body Pole, unit vector.
#   If "q" precedes Pole, PM or Eqx, then a constant pole and spin rate is
#   assumed for quicker computations.
#All vectors are output in ecliptic and mean equinox of J2000 reference frame.
#
#"bodies" is a numerical or cell row vector, or a comma or semicolon delimited
#   string of body numbers or names.
#   Small-body designations are designated with a character string of digits,
#   e.g '3' returns Juno, while [3] returns Earth-Moon Barycenter.
#   The second row is reserved for the "with respect to" body for ephemeris calls.
#   If no second row is input or is NaN, then the "wrt" body is the central body
#
#"dict" or "dictionary" is an array that saves data for use later similar
#   to Matlab's persistent variables. It can be passed in as an empty value
#   "Dict()" or a filled dictionary. It will update after a call and overwrite.
#
#"times" is days past noon 1/1/2000 (J2000) for ephemeris and orientation.
#   times can be a numerical array or 3-element cell of {start time, end time,
#   delta time}. If the number of time elements equals the number of body
#   entries then times and bodies are mapped 1-to-1, otherwise data is computed
#   at every time and concatenated by body along dim 3. The last element of t
#   input is logical, true/1 or false/0. True corresponds to creating an array
#   from [t[1],t[end],dt] while false corresponds to calculating ephemeris at
#   each point of input. If you only have three inputs, the code may think you
#   have an array.
#
#"options" is a 3-element logical array for [save ssd keep_spline] where
#save = true saves ephemeris and orientation data into the (boddat.m path)/ephem
#       /directory and other data in (boddat.m path)/bodydata.mat.
#ssd = true checks Horizons and SSD for data first, while ssd = 0 checks for data
#       in memory, /ephem/ or bodydata.jld before resorting to Horizons.
#   The default time span for saved ephemeris and orientation data is
#   1/1/2000-1/1/2060 (unless span available on Horizons is shorter). Planetary
#   ephemeris is saved in 0.5-day intervals, while moon and orientation data
#   are saved in 0.1-day intervals.
#keep_spline = true keeps the ephemeris and orientation splines in memory for
#   subsequent calls.
#   Note: creates bodydata.jld to save spline data.
#options is optional; the default is [true false true].
#
#Example: Ephemeris +/- one week of Apophis close approach in 0.1-day increments
#t = Dates.datetime2julian(DateTime(2029,04,13,22)) -
#      Dates.datetime2julian(DateTime(2000,01,01,12)) + collect(-7:.1:7)
# Pos,PM,St,n = boddat(["R","qpm","type","num"],["Apophis" "Earth"],Dict(),[t; false])
# currently errors out because of putsbmb needing to push variables to a list instead
# of adding to an already allocated array

#   (It can take up to a minute for Horizons to compute the ephemerides.)
#   Pos = the position of the asteroid Apophis and Earth with respect to the sun
#   (their central body). Pos[:,:,1] is Apophis and Pos[:,:,2] is Earth.
#   PM = the prime meridian of Apophis and Earth. Currently no orientation data
#   is available for Apophis, so PM[:,:,1]=NaN. The "q" in "qpm" flags a
#   constant pole and spin rate for decreased computation time (and accuracy).
#   St = the spectral type of Apophis. St[2] is empty because Earth doesn't
#       have a type.
#   n = the body numbers of Apophis and Earth
#
#Example (con'd): Using saved data
#RV = boddat("X",["Apophis", "Earth"],Dict(),[t;false]); RV = RV[1]
#
#   (The run time should be a few orders of magnitude faster.)
#   RV = the state of Apophis with respect to Earth computed from saved
#       polynomials.
#   Note that the "bodies" input has elements in the second row (transposed
#   from above example)indicating 'Earth' as the wrt body, so
#   RV[1:3,:] = Pos[:,:,1]-Pos[:,:,2]
#
#Example (con'd): Bypassing saved data
#RV_h = boddat("X",[2099942, 399],Dict(),[t[1],t[end],0.1,true],[false true])
#
#   RV_h = the state of Apophis with respect to Earth directly from Horizons.
#   The zero in the options input indicates to not overwrite any saved data, and
#   the one in the options input indicates to skip saved data and read directly
#   from Horizons. When using Horizons it is typically faster (but not
#   necessary) to input times as {start, end, delta}, rather than a list of
#   discrete times.

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
(typeof(bvars) == String) && (varargout = Array(Any,1))
(typeof(bvars) != String) && (varargout = Array(Any,length(bvars)))

if isempty(dict)
  sb = SB([],[],[],[],[],[],[],[],[],[],[])
  mb = MB([],[],[],[],[],[],[],[],[],[])
  tephf = TEPHF([],[],[],[])
  dict["sb"] = sb
  dict["mb"] = mb
  dict["tephf"] = tephf
end

#peek into params (may be useful for troubleshooting)
if isempty(bi)
  if (typeof(bvars) == String)
    varargout = param(dict,bvars)
  else
    for ii in 1:length(varargout)
      varargout[ii] = param(dict,bvars[ii])
      # varagout[ii] = [] # Return empty if don't want the name of the variables
    end
  end
  return varargout
end
#if there's no time, the time is now
if isempty(t)
  t = [Dates.datetime2julian(now())-2451545.]
  typ = false
else
  typ = t[end]
  if (typ==1)||(typ==0)||(typeof(typ)==Bool)
    typ = typ!=0
    t = t[1:end-1]
  else
    warn("Flag not included at end of time array to create interval")
    typ = false
  end
  if size(t,2)>size(t,1)
    t = transpose(t[1:end-1])
  end
  (typeof(t[1])==Int64) && (t = t*1.)
end

#default options is [true false true], replace with any input and save in params
opts=[true false true]
opts[1:length(opti)]=opti
param(dict,"save",opts[1])
param(dict,"ssd",opts[2])
param(dict,"keeps",opts[3])
#boddat.m directory, tells where to save data
if isempty(param(dict,"bdir"))
  bdir = pwd()
  param(dict,"bdir",bdir)
end
#parameter outputs, save in params as character array
if isempty(param(dict,"bva"))
  ov = ["Pole","PM","Eqx"]
  bva = ["R","V","X","A","RTN","number","CB","Name","GM","radius","triaxial",
        "j2","rotation_period","reference","date","mag","occ","types",
        "close_approach"]
  bva = [bva;ov]
  ov2 = identity(ov)
  for jj in 1:3
    for ii in 1:length(ov)
      (jj == 1) && (ov2[ii] = ov[ii]*"dcm")
      (jj == 2) && (ov2[ii] = "q"*ov[ii])
      (jj == 3) && (ov2[ii] = "q"*ov[ii]*"dcm")
    end
    bva = [bva;ov2]
  end
  param(dict,"bva",bva)
end
param(dict,"mbmod",false)
param(dict,"sbmod",false)

#get body numbers as unique identifiers, ephemeris (be) can be 2 x n,
# otherwise first row of b is used
ii = find(x->(typeof(x)==SubString{String}),bi)
for jj in ii
   bi[jj] = bi[jj]*""
end

if (isempty(find(x->(typeof(x)==String),bi)))
  be=identity(bi)
else
  be=getnum(bi,dict)
end
b = zeros(Int64,length(be))
for ii in 1:length(be)
  b[ii]=convert(Int64,be[ii])
end #if size(be,1)>2;be=be';end;b=be(1,:);

#input body parameter variables, if not a cell array then cellstr or split by
# comma, semicolon, space & or |
#if ~iscell(bvars);if size(bvars,1)>1;bvars=cellstr(bvars);else;
# bvars=regexp(bvars,'[,;\s&|]','split');
# bvars(cellfun('isempty',bvars))=[];end;end
#if strcmpi(bvars{1},'all');bvars=cellstr(params('bva'));end
#for each parameter input call the appropriate funciton and write
#result to varargout
if typeof(bvars) == String
  bvars = [identity(bvars)]
  nvars = 1
else
  nvars = length(bvars)
end
for bvr in 1:nvars
  bvi = []
  bvar=bvars[bvr]
  (typeof(bvar)==Char) && (bvar = string(bvar))
  if length(bvar)==1
    inBVA = find(x->(contains(x,bvar)),param(dict,"bva")[1:4])
  else
    inBVA = find(x->(contains(x,bvar)),param(dict,"bva"))
  end
  (!isempty(inBVA)) && (inBVA = inBVA[1])
  #find out which member of bva matches input
  if isempty(inBVA)
    bvar = uppercase(bvar)
    inBVA = find(x->(contains(x,bvar)),param(dict,"bva"))
    (!isempty(inBVA)) && (inBVA = inBVA[1])
    if isempty(inBVA)
      bvar = lowercase(bvar)
      inBVA = find(x->(contains(x,bvar)),param(dict,"bva"))
      (!isempty(inBVA)) && (inBVA = inBVA[1])
      if isempty(inBVA)
        bvar = ucfirst(bvar)
        inBVA = find(x->(contains(x,bvar)),param(dict,"bva"))
        (!isempty(inBVA)) && (inBVA = inBVA[1])
      end
    end
  end
  (isempty(inBVA)) && (inBVA = 0)
  if inBVA > 4
    T1=uniqstr(param(dict,"bva"),bvar)
    bvi = T1[1]
  elseif inBVA <= 4
    bvi = inBVA
  end
  if bvi == 0
    warn("no match found for ",bvar)
    varargout[bvr]=NaN*ones(size(t))
    continue
  end
  #match bvi
  if bvi<5#ephemeris, {R V X A}
    if (!param(dict,"ssd"))||(!isdefined(:Xeph)&&(!isdefined(:vars)))
      Xeph=ephem(b,t,typ,bvi,dict)
    end

    a=Xeph
    if bvi==1
      a=a[1:3,:,:]
    elseif bvi==2
      a=a[4:6,:,:]
    end
  elseif bvi==5
    Xeph=ephem(b,t,typ,3,dict)
    R=Xeph[1:3,:,:]
    V=Xeph[4:6,:,:]
    uH = zeros(R)
    uR = zeros(R)
    uHxR = zeros(R)
    for jj in 1:size(R,3)
      for ii in 1:size(R,2)
        H = cross(R[:,ii,jj],V[:,ii,jj])
        Hm = sqrt(dot(H,H))
        HxR = cross(H,R[:,ii,jj])
        HxRm = sqrt(dot(HxR,HxR))
        Rm = sqrt(dot(R[:,ii,jj],R[:,ii,jj]))
        for kk in 1:3
          uH[kk,ii,jj] = H[kk]/Hm
          uHxR[kk,ii,jj] = HxR[kk]/HxRm
          uR[kk,ii,jj] = R[kk,ii,jj]/Rm
        end
      end
    end
    a = [uR; uHxR; uH]
  elseif bvi==6
    a=getnumout(b,dict) #number
  elseif bvi==7#central body, matches bi for name vs number output
    cb=10*ones(Int64,length(b))
    for jj in 1:length(b)
      if (b[jj]>10)&&(b[jj]<1000)&&(mod(convert(Int64,real(b[jj])),100)!=99)
        cb[jj] = floor(b[jj]/100)*100 + 99
      end
      (imag(b[jj]) != 0) && (cb[jj] = real(b[jj]))
    end
    T1 = trues(length(bi))
    for ii in 1:length(bi)
      if (typeof(bi[ii])!=Float64)&&(typeof(bi[ii])!=Int64)
        T1[ii]=false
      end
    end
    jj = find(x->(x),T1)
    if isempty(jj)
      a = cb
    else
      a = identity(bi)
      for ib in 1:length(bi)
        if T1[ib]
          a[ib] = cb[ib]
        else
          T2 = getnam(cb[ib],dict)
          a[ib] = T2[1]
        end
      end
    end
  #8--19 calls a function called getxxx, where xxx is first three letters
  #       of matching bva parameter. input is xxx and b
  elseif bvi==8
    a=getnam(b,dict) #name, works for NonBodyControlPoints
  elseif bvi<20
    fn=param(dict,"bva")
    fn=fn[bvi]
    a = getfn(fn,b,dict)
  # ii=~imag(b);
  # if any(ii);fn=params('bva');fn=fn(bvi,1:min(3,end));
  # eval(['a(:,ii)=get' fn '(b(ii));']);end
  # if any(~ii);if iscell(a);a(~ii)={''};else;a(:,~ii)=nan;end;end
  #20--31 is orientation
  #ll=1: Pole
  #ll=2: [Node;Q;Pole] where Node=cross([0;0;1],Pole), Q=cross(Pole,Node)
  #ll=3: Prime Meridian
  #ll=4: [Node;Q;Pole] where Node=Prime Meridian
  #ll=5: Equinox, cross(Pole,H)
  #ll=6: [Node;Q;Pole] where Node=Equinox

  else
    fn=param(dict,"bva")
    fn=fn[bvi]
    if contains(fn,"Pole")
      ll = 2
    elseif contains(fn,"PM")
      ll = 4
    elseif contains(fn,"Eqx")
      ll = 6
    else
      ll = []
    end
    (!contains(fn,"dcm")) && (ll = ll - 1) #check for dcm
    if fn[1]=='q' #check for quick calculation flag
      a=qorient(b,t,typ,ll,dict)
    else
      a,_=orient(b,t,typ,ll,dict)
    end

  end#if bvi<5
  varargout[bvr]=a #write output to varargout
end
#save data to bodydata
if param(dict,"save")
  @static is_windows()? (fn= string(param(dict,"bdir"),"\\bodydata.jld")):(fn=
              string(param(dict,"bdir"),"/bodydata.jld"))
  if !isfile(fn)
    fid = jldopen(fn,"w") # Must use HDF5 and JLD package
    boddict = Dict()
    boddict["mb"] = MB([],[],[],[],[],[],[],[],[],[])
    boddict["sb"] = SB([],[],[],[],[],[],[],[],[],[],[])
    fid["boddict"] = boddict
    close(fid)
  end
  boddict = Dict()
  if param(dict,"mbmod")
    mb,_=getsbmb("mb",dict)
    boddict["mb"] = mb
    fid = jldopen(fn,"r+")
    if exists(fid,"boddict")
      dT = read(fid["boddict"])
      boddict["sb"] = dT["sb"]
    end
    close(fid)
    fid = jldopen(fn,"w")
    fid["boddict"] = boddict
    close(fid)
  end
  if param(dict,"sbmod")
    sb,_=getsbmb("sb",dict)
    fid = jldopen(fn,"r+")
    boddict["sb"] = sb
    if exists(fid,"boddict")
      dT = read(fid["boddict"])
      boddict["mb"] = dT["mb"]
    end
    close(fid)
    fid = jldopen(fn,"w")
    fid["boddict"] = boddict
    close(fid)
  end
end
return varargout
end

function getnum(bi::Array,dict::Dict)
#get body number, only called if bi isn't numerical
(typeof(param(dict,"mb"))!=MB) && (param(dict,"mb",MB([],[],[],[],[],[],[],[],[],[])))
mb = param(dict,"mb")
# sav = false

#change into string, bi is used in central body output
(typeof(bi)==String) && (b = Array(String,1))
(typeof(bi)!=String) && (b = Array(Any,length(bi)))
for ib in 1:length(bi)
  if typeof(bi[ib])!=String
    b[ib] = identity(convert(Int64,bi[ib]))
    continue
  else
    bib = identity(bi[ib])
  end
  if contains(bi[ib],"_CP")
    nbcp = true
    replace(bib,"_CP","")
  else
    nbcp = false
  end
  if typeof(bib)!=String
    b[ib] = bib #already a number
  else #match a string
    x,_=getsbmb("mb",dict)
    if isempty(mb.Names)&&(isempty(x.Names)||param(dict,"ssd"))
        x=getmb(dict)
    end
    mb=x #get mb from saved data unless empty or ssd flagged
    In=uniqstr(x.Names,bib,1) #check whole word in major body data
    if isempty(In)
      x,_=getsbmb("sb",dict)
      sb=x
      if !isempty(fieldnames(x))
        In=uniqstr(x.Names,bib,1) #check whole word in small body data
      else
        In = []
      end
      if (isempty(In))||(param(dict,"ssd"))
        uniqsb(bib,dict)
        x,_=getsbmb("sb",dict)
        sb=x
        In=uniqstr(x.Names,bib,1)
      end#check ssd webpages
      if isempty(In)
        warn("No whole-word match found for ",bib)
        x=mb;In=uniqstr(x.Names,bib,2);#check fragment in mb
        if isempty(In)
          x=sb;In=uniqstr(x.Names,bib,2);#check fragment in sb
          if isempty(In)
            x=getmb(dict);In=uniqstr(x.Names,bib)
          end#mb potentially skipped
        end#sb any
        if !isempty(In)
          println("> Found it! ",bib," set to "x.Names[In,:])
        end
      end#sb word
    end#mb any
    if isempty(In)
      warn("No match found for ",bib)
      b[ib]=NaN
    else
      b[ib]=x.numbers[In]
    end
  end#mb word
(nbcp) && (b[ib]=b[ib]+1im)
if (typeof(b[ib])!=Float64) || (typeof(b[ib])!=Float64)
  b[ib] = b[ib][1]
end
end#for ib
param(dict,"mb",mb)
return b
end

function getnumout{T}(bi::AbstractArray{T},dict::Dict)

bo = zeros(length(bi))
for ib in 1:length(bi)
  b1=identity(bi[ib])
  if imag(b1) != 0
    nbcp=true
    b1=real(b1)
  else
    nbcp=false
  end
  if isnan(b1)
    n=b1
  else
    n=getx(b1,"numbers",dict) #check saved data (& no ssd flag)
    if !isempty(n)#found it
    elseif b1<1e6
      x=getmb(dict)
      n=find(y->(y==b1),x.numbers) #check list of major bodies
      (!isempty(n)) && (n=x.numbers[n])
      (isempty(n)) && (n=NaN)
    else
      n=getsb(b1,dict)
      n=n[1]#check ssd
    end
  end#isnan
  (n==0) && (warn("No match found for ",b1))#nomatch
  (nbcp) && (n=n+1im)
  bo[ib]=n
end
return bo
end

function getnam(bs,dict::Dict)
#get body name
#write name that matches each number input to cell
nn=Array(AbstractString,length(bs))
for ib in 1:length(bs)
  b=bs[ib]
  if imag(b) != 0
    nbcp=true
    b=real(b)
  else
    nbcp=false
  end
  n=getx(b,"Names",dict) #check saved data (& no ssd flag)
  if !isempty(n)#found it
  elseif b<1e6
    x=getmb(dict)
    n=find(y->(y==b),x.numbers)#check list of major bodies
      if !isempty(n)
        n=x.Names(n)
      end
  else
    n,_,_=getsb(b,dict)  #check ssd
    (!isnan(n[1])) && (n=n[2])
    (isnan(n[1])) && (n=[])
  end

  if isempty(n)
    n="NULL"
    warn("No match found for ",b)
  end#nomatch
  #format small body name, tokens in regexp are {number, name, (YYYY A1)}
  if b>1e6;
    n1=search(n,r" ")
    n2 = n[n1[end]+1:end]
    nT1 = n[1:n1[end]-1]
    n1 = search(n2,"(")
    if n1[1] == 1
      nT2 = []
    else
      nT2 = n2[1:n1[end]-1]
    end
    nT3 = n2[n1[end]:end]
    #n=n{1};#n=deblank(n{1});
    #if named, go with name, if numbered go with "# (provisional)",
    # otherwise output provisional designation
    if !isempty(nT2)
      n=nT2
    elseif !isempty(nT1)
      n=n
    else
      n=nT3[2:end-1]
    end
  end
  (nbcp) && (n=string(n,"_CP"))
  nn[ib]=n
end#for ib
return nn
end

function getgm(bs::AbstractArray{Int64},dict::Dict)
# GM from ephemeris header constants unless = 0, then estimate from
#   density and size
gmx=zeros(length(bs))
for jj in 1:length(bs)
  b=identity(bs[jj])
  (imag(b)!=0) && (gm=NaN)
  (imag(b)==0) && (gm=getx(b,"gm",dict)) #check saved data (& no ssd flag)
  if (!isempty(gm)) && (gm !=NaN)#found it
  elseif b<1e6
    ed=getephdat(floor(Int64,b),dict) #major body, check ephemeris header for
    # GM>1e-9
      if (!isempty(ed)) && (ed["gm"*string(floor(Int64,b))][1]>1e-9)
        gm=ed["gm"*string(floor(Int64,b))][1]
        putsbmb(b,"gm",gm,dict)
      #read general satellite page, get first # after >name < and ">"
      else
        ed=getsatdat(dict)
        n=getnam(b,dict)
        T1 = Regex(join([">",n[1],"[\\s<]"],""))
        r1 = search(ed,T1)
        if (r1 != 0:-1)
          r2 = search(ed[r1[1]:end],r".([\d.])+[\s&<]")
          gm = parse(Float64,ed[r2[1]+r1[1]:r2[end]+r1[1]-2])
          putsbmb(b,"gm",gm,dict)
        else
          gm=NaN
          warn("no GM found for",n[1])
        end ##ok<*WNTAG>
      end
  else
    gm,_,_=getsb(b,dict)#small body
    if !isnan(gm[1])
      gm=gm[3]
    else
      gm=NaN
      warn("body ",b," not found for GM data")
    end
  end
  gmx[jj]=gm
end#jj
return gmx
end

function getj2(bs::AbstractArray{Int64},dict::Dict)
j2x=zeros(Float64,size(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if any(b in [10; collect(3:9)*100+99; 301])
  #J2 currently available for only a few bodies from ephemeris header
    if imag(b)!=0
      j2=NaN
    else
      j2=getx(b,"j2",dict)
    end
    if isempty(j2)
      ed=getephdat(b,dict)
      j2=ed["j2_"*string(b)]
      putsbmb(b,"j2",j2,dict)
    end
  elseif condition
    warn("J2 is only available for the following bodies: ",
          [10; collect(3:9)*100+99; 301])
    j2=NaN
  end
  j2x[jj]=j2[1]
end
return j2x
end

function getrad(bs::AbstractArray{Int64},dict::Dict)
rx=zeros(Float64,length(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if imag(b)!=0 #check saved data (& no ssd flag)
    r=NaN
  else
    r=getx(b,"rad",dict)
  end
  if !isempty(r)#found it
  elseif (b<1e6) && ((mod(b,100)==99)||(b==10))
    ed=getephdat(floor(Int64,b),dict)
    r=ed["rad"*string(floor(Int64,b))]
    if isempty(r)
      r=gettri(b,dict)
      r=r[1]
    end
    putsbmb(b,"rad",r,dict) #planet or sun, ephemeris header
    #satellite, read general satellite page, get second # after >name < and ">"
  elseif b<1e6
    ed=getsatdat(dict)
    n=getnam(b,dict)
    T1 = search(ed,"eft>"*n[1])
    if T1 != 0:-1
      T2 = search(ed[T1[end]:end],r"[\s<].*?>[\d.]") + T1[end]
      T3 = search(ed[T2[end]:end],r"[\s&<].*?>([\d.]+)") + T2[end]
      r = ed[T3[1]:T3[end]]
      T4 = search(r,r"\d")
      r = r[T4[1]:end]
      (r[end]=='&') && (r = r[1:end-1])
      (r[end]=='<') && (r = r[1:end-1])
      putsbmb(b,"rad",r,dict)
    else
      r=NaN
      warn("no radius found for ",n[1])
    end
  else
    r,_,_=getsb(b,dict)#small body
    if !isnan(r[1]);
      r=[4]
    else
      r=NaN
      warn("body ",b," not found for radius data")
    end
  end
  (length(r)>1) && (r = r[1]) # could be removed wiht other error checking
  rx[jj]=r[1] #dev: r->r[1]
end#for
return rx
end

function gettri(bs::Int64,dict::Dict)
trix=zeros(3)
for jj=1:length(bs)
  b=bs[jj]
  if imag(b)!=0
    r=NaN*[1 1 1]
  else
    r=getx(b,"tri",dict)
  end#check saved data (& no ssd flag)
  if isempty(r)
    r=getpck(b,"RADII",dict) #check pck file
    if (isempty(r)) && (b>1e6)#small body
      #extent is flag for triaxial radii,convert to #,
      #   make axisymmetric if only 2#s
      nn,ee,ss=getsb(b,dict)
      f = open(ss); rd = readstring(f); close(f)
      r1 = search(rd,r">extent<.*?>")
      if r1 != 0:-1
        off = r1[end]
        r1 = search(rd[off:end],r">extent<.*?>")
        r2 = search(rd[off+r1[1]:end],r">\d")
        r3 = search(rd[off+r1[1]+r2[1]:end],"<")
        rs = rd[off+r1[1]+r2[1]:off+r1[1]+r2[1]+r3[1]-2]
        r = split(rs,'x')
        for ii in 1:length(r)
          r[ii] =parse(r)
        end
        (length(r) == 2) && (r =[r[1] r])
      end
    end
    if isempty(r) #see if there's a single radius value
      T1 = getrad([b],dict)
      r = T1*[1 1 1]
    end
    (!isnan(r[1])) && (putsbmb(b,"tri",r,dict))
  end#if
  trix[:,jj] = r
end#for jj
return trix
end

function getrot(bs::AbstractArray{Int64},dict::Dict)
rx=zeros(Float64,size(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if imag(b)!=0
    r=NaN
  else #check saved data (& no ssd flag)
    r=getx(b,"rot",dict)
  end
  if isempty(r)
    r=getpck(b,"PM",dict) #check pck file
    if !isempty(r)
      r=360/r[2]
      putsbmb(b,"rot",r,dict) #convert to days
    elseif b>1e6 #small body
      r,_=getsb(b,dict)
      if !isnan(r[1])
        r=r[5]
      else
        r=NaN
      end
    else
      r=NaN
    end
    (isnan(r)) && (warn("body ",b," not found for rotation data"))
  end#if
  rx[jj]=r
end#for jj
return rx
end

function getmag(bs::AbstractArray{Int64},dict::Dict)
hx=zeros(Float64,size(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if (b<1e6)||(imag(b)!=0)
    hx[jj]=NaN
    continue
  end#only for small bodies
  h=getx(b,"h",dict) #check saved data (& no ssd flag)
  if isempty(h)
    h,_=getsb(b,dict)  #check ssd
    if !isnan(h[1])
      h=h[6]
    else
      h=NaN
      warn("body ",b," not found for absolute magnitude")
    end
  end#if
  hx[jj]=h
end#for jj
return hx
end

function getocc(bs::AbstractArray{Int64},dict::Dict)
hx=zeros(Float64,size(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if (b<1e6)||(imag(b)!=0) #only for small bodies
    hx[jj]=NaN
    continue
  end
  h=getx(b,"occ",dict) #check saved data (& no ssd flag)
  if isempty(h)
    h,_=getsb(b,dict)#check ssd
    if !isnan(h[1])
      h=h[7]
    else
      h=NaN
      warn("body ",b," not found for condition code")
    end
  end#if
  hx[jj]=h
end#for jj
return hx
end

function gettyp(bs::AbstractArray{Int64},dict::Dict)
hx=Array(String,size(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if (b<1e6)||(imag(b)!=0) #only for small bodies
    hx[jj]=""
    continue
  end
  h=getx(b,"type",dict)#check saved data (& no ssd flag)
  if isempty(h)
    h,_=getsb(b,dict)#check ssd
    if !isnan(h[1])
      h=h[8]
    else
      h=NN
      warn("body ",b," not found for spectral type")
    end
  end#if
  hx[jj]=strip(h,'_')#underscore denotes checked but nothing found
end#for jj
return hx
end

function getclo(bs::AbstractArray{Int64},dict::Dict)
cada=Array(Any,size(bs))
#ssd=params('ssd');params('ssd',true)#for getnum
for jj=1:length(bs)
  b=bs[jj]
  if imag(b)!=0
    cada[jj]=""
    continue
  end
  if b<1e6 #only for small bodies
    cada[jj]=[]
    continue
  end
  s=download("https://ssd.jpl.nasa.gov/sbdb.cgi?sstr="*string(b)*";cad=1")
  f = open(s); rd = readlines(f); close(f)
  if jj/100==round(Int64,jj/100)
    println("Close approach progress: ",floor(jj/length(bs)*100,0),"%")
  end
  cls =[]
  clt = []
  for ii in 1:length(rd)
    (contains(rd[ii],"Close-Approach Data")) && (cls = [cls; ii])
    (contains(rd[ii],"<b>Close-Approach Data")) && (clt = [clt;ii])
  end
  if isempty(cls)
    warn("no close approach data for body ",b)
    cada[jj] = []
    continue
  end
  (length(clt)>1) && (warn("close-approach data might not be correct"))
  clt = clt[1]
  ii = find(x->(x==clt),cls)
  tst = cls[ii[1]]
  ten = cls[ii[1]+1]
  lne = []
  for ii in tst:ten
    (contains(rd[ii],"</tr>")) && (lne = [lne; ii])
  end
  cad = zeros(length(lne)-3,3)
  LastNum = ""
  for ii in 1:length(lne)-3
    lns = match(r"size=\"-2\">",rd[lne[ii]+1])
    ln = match(r"</font>",rd[lne[ii]+1])
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
      T2 = getnum([T1],dict)
      cad[ii,2] = T2[1]
    end
    LastNum = identity(T1)

    lns = match(r"size=\"-2\">",rd[lne[ii]+4])
    ln = match(r"</font>",rd[lne[ii]+4])
    cad[ii,3] = parse(Float64,rd[lne[ii]+4][lns.offset+10:ln.offset-1])
  end
    # find <tr> in beginning of line and </tr>\n at end of lines
    # Looking at first, third, and fourth column (date, body, nominal dist)
  cada[jj]=cad
end#for jj
#params('ssd',ssd)
#if jj==1;cada=cada{:};end
return cada
end

function getref(bs::AbstractArray{Int64},dict::Dict)
hx=Array(String,size(bs))
for jj=1:length(bs)
  b=bs[jj]
  if imag(b)!=0
    hx[jj]=""
    continue
  end
  if b>1e6
    h=getx(b,"ephref",dict)#check saved data (& no ssd flag)
  else
    h = ""
  end
  if (isempty(h)) && (param(dict,"ssd")) #see what current file is
    if b<1e6
      tephf=getef(dict)
      ii = find(x->(x==string(b)),tephf.numbers)
      h=tephf.f[ii[1]] #get ephemeris from list
    else
      h,e=getsb(b,dict)
      if !isnan(h[1])
        h=e
      else
        h=""
      end
    end#get reference from ssd
  end#if, no ephemeris found
  (isempty(h)) && warn("file of last ephemeris for ",b," is unknown")
  hx[jj]=h
end#for jj
return hx
end

function getdat(bs::AbstractArray{Int64},dict::Dict)
hx=zeros(Float64,size(bs))
for jj in 1:length(bs)
  b=bs[jj]
  if imag[b]!=0
    hx[jj]=NaN
    continue
  end
  h=getx(b,"ephdate",dict) #check saved data (& no ssd flag), no ssd analogue
  if isempty(h)
    h=NaN
    warn("date of last ephemeris for ",b," is unknown")
  end
  hx[jj]=h
end#for jj
return hx
end

function qorient(bb::AbstractArray{Int64},tt::AbstractArray{Float64},
                                              typ::Bool,ll::Int64,dict::Dict)
#quick orientation calculation
#ll=1: Pole
#ll=2: [Node;Q;Pole] where Node=cross(Earth pole,Pole), Q=cross(Pole,Node)
#ll=3: Prime Meridian
#ll=4: [Node;Q;Pole] where Node=Prime Meridian
#ll=5: Equinox, cross(Pole,H)
#ll=6: [Node;Q;Pole] where Node=Equinox

sb=size(bb)
if (sb[1]>2) && (sb[2]>2)
  bb=bb[1,:]
elseif sb[1]>2
  bb=bb'
end

tb=false
if typ #match t and b 1-to-1 if numel(tt)==numel(bb)
  tt = collect(tt[1]:tt[3]:tt[2])
else
  tt=transpose(tt[:])'
  tb=length(tt)==length(bb)
end

# Find unique bodies and when they occur
it = collect(1:length(bb))
ct = collect(2:length(bb))
for ii in 1:length(bb)-1
  sp = 1
  T1 = identity(ct)
  for jj in ct
    if bb[ii] == bb[jj]
      splice!(T1,sp)
      it[jj] = ii
      sp = sp - 1
    end
    sp = sp + 1
  end
  ct = identity(T1)
  (isempty(ct)) && (continue)
end

cl2=ceil(Int,ll/2)
if ll == 1 || ll == 3
  n = (3,length(tt))
elseif ll == 4 || ll == 2
  n = (9,length(tt))
elseif ll == 5
  n = (3,1)
end

if tb
  xx =Array(Float64,n)
else
  xx = Array(Float64,(n[1],n[2],length(bb)))
end

for ibu in 1:length(bb)
  if it[ibu] != ibu #match times for each unique body
    if tb  && (tt[ibu]==tt[it[ibu]])## Change the input and copy for other values
      xx[:,ibu] = xx[:,itp[ibu]]
      continue
    elseif !tb
      xx[:,:,ibu] = xx[:,:,itp[ibu]]
      continue
    end
  end
  b=bb[ibu]
  tb? (t=tt[ibu]) : (t=tt)

  pp=getx(b,"pp",dict)
  if isempty(pp)#check saved data (& no ssd flag)
    adp=[getpck(b,"POLE_RA",dict)';getpck(b,"POLE_DEC",dict)';
        getpck(b,"PM",dict)']*pi/180 #get data from pck file
    if length(adp)==9
      a=adp[1]; d=adp[2]; p=adp[3]; o=84381.448/3600*pi/180
      #obliquity @J2000: Lieske, J., "Precession Matrix Based on IAU (1976)
      # System of Astronomical Constants"
      #rotate from EME to EMO
      #basis orthogonal to pole, N is pole node: pole x [0;0;1]
      P=[cos(d)*cos(a); sin(d)*sin(o)+cos(d)*cos(o)*sin(a);
          cos(o)*sin(d)-cos(d)*sin(o)*sin(a)] #pole
      N=[P[2];-P[1];0]
      N=N/sqrt(dot(N,N))
      Q=cross(P,N)
      #basis orthogonal to pole, N is pole node: Earth pole x pole
      #N=cross(qorient(399,0,1),P);N=N/sqrt(sum(N.*N));Q=cross(P,N);
      #prime meridian
      X=[-cos(p)*sin(a)-cos(a)*sin(d)*sin(p);cos(a)*cos(p)-sin(a)*sin(d)*sin(p);
          cos(d)*sin(p)]
      #rotate to EMO # get angle @ J2000
      X=[X[1];cos(o)*X[2]+sin(o)*X[3];cos(o)*X[3]-sin(o)*X[2]]
      p=atan2(dot(X,Q),dot(X,N))
      #get H
      sav=param(dict,"save")
      param(dict,"save",false)
      X,_=ephem1(b,[0.],typ,3,dict)
      param(dict,"save",sav)
      #equinox
      X=cross(P,cross(X[1:3],X[4:6]))
      X=X/sqrt(dot(X,X))
      pe=atan2(dot(X,Q),dot(X,N))
      #pp =[pole vector;PM angle;eqx angle] @ J2000
      pp=[P;p;pe]; putsbmb(b,"pp",pp,dict)
    else
      warn("no analytic orientation data found for",b)
      pp=NaN*ones(5,1)
    end
  end
  X=pp[1:3] #X = pole, cl2 is Pole|PM|eqx output
  if ll>1
    #need X = pole node
    P=identity(X)
    p12=sqrt(P[1]^2+P[2]^2)
    X=[P[2]/p12;-P[1]/p12; 0*p12]
    #if ll>1;P=X;X=cross(qorient(399,0,1),P);X=X/sqrt(sum(X.*X));#need X = pole node
    if cl2>1
      if cl2==2 #need angle
        T1 = getrot([b],dict)
        w=pp[4]+2*pi/T1[1]*t
      elseif cl2==3
        w=pp[5]
        t=0
      end
      #rotate along pole
      T1 = X; X = zeros(3,length(w)); T2 = P; P = zeros(3,length(w))
      for ii in 1:length(w)
        X[:,ii]=T1*cos(w[ii])+cross(T2,T1)*sin(w[ii])
        P[:,ii] = T2
      end
    end
  end
  if ll==2 && b!=399 #J2000 Node DCM
    X = zeros(3,size(P,2))
    T1 = qorient([399],[0.],typ,1,dict)
    for ii in size(P,2)
      X[:,ii]=cross(T1,P[:,ii])
      X[:,ii] = X[:,ii]/sqrt(dot(X[:,ii],X[:,ii]))
    end
  end
  if cl2==div(ll,2) #dcm
    T1 = X; X = zeros(9,size(X,2))
    for ii in 1:size(X,2)
      X[:,ii] = [T1[:,ii];cross(P[:,ii],T1[:,ii]);P[:,ii]]
    end
  end
  (cl2!=2) && (X=repmat(X,size(t))) #match size of t

  if tb
    xx[:,ibu]=X
  else
    xx[:,:,ibu]=X
  end
end#for ibu
return xx
end

function ephem(b::AbstractArray{Int64},t::AbstractArray{Float64},
                    typ::Bool,bvi::Int64,dict::Dict)
#ephemeris
# typ is if t is a linspace (true) or if t is vector (false)
sb=size(b)
if (sb[1]>2) && (sb[2]>2)
  warn("Body input for ephemeris should be 1 x n or 2 x n")
  b=b[1,:]
elseif sb[1]>2
  b=b'
end

#match t and b 1-to-1 if numel(t)==size(b,2)
tb=false
if !typ
  if size(t,1) > 1
    T1 = t[:,1]
    for ii in 2:size(t,2)
      T1 = [T1; t[:,ii]]
    end
    t = transpose(T1)
  end
  tb = (length(t)==size(b,2))
end

if size(b,1)==1
  b =[b; NaN*ones(b)]
end#b(2,:) is wrt body list

it = collect(1:1:size(b,2))
ct = collect(2:1:size(b,2))
for ii in 1:size(b,2)-1
  sp = 1
  T1 = identity(ct)
  for jj in ct
    if (isnan(b[2,ii]))
      if (b[1,ii] == b[1,jj]) && (isnan(b[2,jj]))
        splice!(T1,sp)
        it[jj] = ii
        sp = sp - 1
      end
    else
      if b[:,ii] == b[:,jj]
        splice!(T1,sp)
        it[jj] = ii
        sp = sp - 1
      end
    end
    sp = sp + 1
  end
  ct = identity(T1)
  (isempty(ct)) && (continue)
end

n = 6
(bvi == 4) && (n = 12)
if tb
  X = Array(Float64,(n,size(b,2)))
else
  X = Array(Float64,(n,length(t),size(b,2)))
end

for ibu in 1:size(b,2)
  if it[ibu] != ibu #match times for each unique body
    if tb && (t[ibu]==t[it[ibu]])
      X[:,ibu] = X[:,it[ibu]]
      continue
    elseif !tb
      X[:,:,ibu] = X[:,:,it[ibu]]
      continue
    end
  end
  bib=b[1,ibu]; b0b=b[2,ibu] #bib is target, b0b is wrt
  (typeof(bib)==Float64) && (bib = convert(Int64,bib))
  if isnan(b0b) #wrt cetral body, no need to "tree" body center
    #get state, X, from ephem1 function
    if tb #match 1-to-one, or all times per body
      X[:,ibu],_=ephem1(bib,[t[ibu]],false,bvi,dict)
    else
      X[:,:,ibu],_=ephem1(bib,t,typ,bvi,dict)
    end
  else#wrt some specified body, may need to tree
    if tb
      tt=[t[ibu]]
    else
      tt=identity(t)
      if isempty(sizeof(tt))
        tt = [tt]
      end
    end#match 1-to-one, or all times per body
    if bib==b0b#target=wrt, X=0
        if typ
          tt=collect(t[1]:t[3]:t[2])
        end
        Xb=zeros(6,length(tt))
        Xb0=0
    else#compute ephemeris
      cb=10
      if imag(bib) !=0
        cb=real(bib)
      elseif (bib>10) && (bib<1000) && (mod(real(bib),100)!=99)
        cb=floor(Int64,bib/100)*100+99
      end#central body of target
      cb0=10
      if (b0b>10) && (b0b<1000) && (mod(b0b,100)!=99)
        cb0=floor(Int64,b0b/100)*100+99
      end#central body of wrt body
      lcb=10
      if ((cb!=10)||(cb0!=10)) && ((cb==cb0)||(cb==b0b)||(bib==cb0))
        lcb=max(cb,cb0)
      end#lowest central body (planet or sun)
      Xb =0
      if lcb!=bib
        Xb,_=ephem1(bib,tt,typ,bvi,dict)
        if lcb!=cb
          XT,_=ephem1(cb,tt,typ,bvi,dict)
          Xb=Xb+XT
        end
      end#get target state wrt lcb
      Xb0=0
      if lcb!=b0b
        Xb0,_=ephem1(b0b,tt,typ,bvi,dict)
        if lcb!=cb0
          XT,_ = ephem1(cb0,tt,typ,bvi,dict)
          Xb0=Xb0+XT
        end
      end#get "wrt body" state wrt lcb
    end#if bib

  #X-Xb0 is target state wrt "wrt body" state
    if tb
      X[:,ibu]=Xb-Xb0
    elseif bib != b0b
      X[:,:,ibu]=Xb-Xb0
    elseif size(b,2) == 1
      X = zeros(6,length(t),1)
      X[:,:,ibu]=Xb-Xb0
    else
      warn("Resulting epemeris probably junk. Change ephem function in boddat")
    end#match 1-to-one, or all times per body
  end
end

return X
end

function ephem1(bib::Int64,t::AbstractArray{Float64},typ::Bool,
                  bvi::Int64,dict::Dict)
n=6
if bvi==4
  n=12
end
if imag(bib)!=0 #non-body control point
  if typ
    t=collect(t[1]:t[3]:t[2])
  end
  X=zeros(n,length(t))
  return
end

#get saved data, ideph is list of body numbers, deph is spline data
id=param(dict,"ideph")
d=param(dict,"deph") #clear double so can overwrite as struct

if (param(dict,"ssd"))
  X,t,err=mkeph(bib,t,typ,dict) #go straight to horizons
  #flag if any data is out of bounds
  (length(err)==1) && (err = [err false])
  if err[1]
    warn("Min Horizons ephem time for ",bib," is EXCEEDED")
    #,datestr(err[1]+730486.5))
  elseif err[2]
    warn("Max Horizons ephem time for ", bib," is EXCEEDED")
    #,datestr(err[2]+730486.5))
  end
  (bvi==4) && (warn("no accel or jerk from Horizons"))
  return X,t #return state
end

if (bib>10) && (bib<1e3) && (mod(bib,100)!=99)
  cb=floor(Int64,bib/100)*100+99
else
  cb=10
end#central body, treat planets and satellites differently
ii=find(x->(x==bib),id)#see if saved data exists
if isempty(ii)
  @static is_windows()? (efile = param(dict,"bdir")*"\\ephem\\"*string(bib)".jld"):(efile =
              param(dict,"bdir")*"/ephem/"*string(bib)*".jld")
  if (!isfile(efile)) && (bib>1e6)
    bib,_=getsb(bib,dict)
    bib=bib[1]
    @static is_windows()? (efile = param(dict,"bdir")*"\\ephem\\"*string(bib)".jld"):(efile =
                param(dict,"bdir")*"/ephem/"*string(bib)".jld")
                #see if number changed
  end
  if !isfile(efile) #te is saved time span, di is data to save, don't save nans
    te=[-1,60*365.25,.1+.4*(cb==10)+.5*(bib>1e6)]
    di,td,_=mkeph(bib,te,true,dict)
    di = [td;di]
    if isempty(di)
      (typ) && (t=collect(t[1]:t[3]:t[2]))
      X=NaN*ones(n,length(t))
      return X,t
    end#no good data
  else
    f=load(efile)
    di=identity(f["d"])
  end#file exists, read data, 1st two entries is data size

  if cb!=10#gives better accuracy for Moons, but takes ~twice time
    td=di[1,:]
    R=di[2:4,:]
    V=di[5:7,:]
    H = zeros(R)
    for ii in 1:size(R,2)
      T = cross(R[:,ii],V[:,ii])
      H[:,ii] = T/sqrt(dot(T,T))
    end
    #mH=repmat(mean(H,2),size(td));[v e]=eig((H-mH)*(H-mH)');
    #[~,mi]=min(diag(e));H=v(:,mi);#Rotate to Local Laplace plane
    H=mean(H,2)
    H = H/sqrt(dot(H,H))

    h12=sqrt(H[1]^2+H[2]^2)
    N=[-H[2]/h12; H[1]/h12; 0]
    Q=[-H[3]*N[2]; H[3]*N[1]; h12]
    dcm=[N Q H]
    R=transpose(dcm)*R
    zr=R[3,:]
    R=R[1:2,:]
    r = zeros(zr)
    for ii in 1:size(R,2)
      r[ii] = sqrt(dot(R[:,ii],R[:,ii]))
    end
    qr = atan2(R[2,:],R[1,:])

    dqr = diff(qr)
    T1 = falses(length(dqr))
    for ii in 1:length(dqr)
      T1[ii] = abs(dqr[ii])>pi
    end
    jj = find(T1)
    for ii in 1:length(jj)
      if (ii == length(jj)) && (jj != length(qr))
        qr[jj[ii]+1:end] = qr[jj[ii]+1:end] + ii*2*pi
      else
        qr[jj[ii]+1:jj[ii+1]] = qr[jj[ii]+1:jj[ii+1]] + ii*2*pi
      end
    end

    V=transpose(dcm)*V
    zv=V[3,:]
    V=V[1:2,:]
    v = zeros(zr)
    qv = zeros(zr)
    for ii in 1:size(V,2)
      v[ii] = dot(V[:,ii],R[:,ii])/r[ii]
      qv[ii] = (R[1,ii]*V[2,ii]-R[2,ii]*V[1,ii])/r[ii]^2
    end
    di[2:end,:]=[r';qr';zr';v';qv';zv']
  else
    dcm=[]
  end

  ii=size(id,1)+1
  if isempty(id)
    id = zeros(1,1)
    id[1,1] = bib
  else
    id=[id; bib]#add body to list
  end

  d = param(dict,"deph")
  (isempty(d)) && (d = Dict())
  d[ii] = Dict()
  if param(dict,"keeps") #save spline
    RVspline6!(di,d[ii])
    d[ii]["dcm1"]=dcm
    param(dict,"ideph",id)
    param(dict,"deph",d)
  else
    RVAspline!(di,d[ii])
    d[ii]["dcm1"]=dcm
    param(dict,"ideph",[])
    param(dict,"deph",[])
  end
else
  ii = ii[1]
end#isempty ii

#convert {start,end,delta} to list,
if typ
  t=collect(t[1]:t[3]:t[2])
end
ts=[]
ist=falses(length(t))
#see if data is out of bounds
tl=d[ii]["breaks"][1]
if minimum(t)<tl
  ii = find(x->(x<tl),t)
  for jj in ii
    ist[jj] = true
  end
  warn("Min time for ",bib," ephem spline is ",Dates.julian2datetime(tl+2451545))
end
tl=d[ii]["breaks"][end]
if maximum(t)>tl
  ii = find(x->(x>tl),t)
  for jj in ii
    ist[jj] = true
  end
  warn("Max time for ",bib," ephem spline is ",Dates.julian2datetime(tl+2451545))
end
if any(ist) #ist is indeces of data outside of spline range
  ts=t[ist]
  t[ist]=[]
end

X=zeros(6,length(t))
if (bvi!=2)||(cb!=10)
  (typeof(t) == Float64) && (t = [t])
  pval!(d[ii],t,1,view(X,1:3,:))
end
if bvi>1
  pval!(d[ii],t,2,view(X,4:6,:))
end
if bvi==4
  T1 = zeros(3,length(t))
  T2 = zeros(3,length(t))
  pval!(d[ii],t,3,T1)
  pval!(d[ii],t,4,T2)
  X=[X; T1; T2]
end

if cb!=10#gives better accuracy for Moons, but takes ~twice time
  dcm=d[ii]["dcm1"]
  c=cos(X[2,:])
  s=sin(X[2,:])
  r=X[1,:]
  if bvi!=2
    for ii in 1:size(X,2)
      X[1:3,ii] = dcm*[r[ii]*c[ii]; r[ii]*s[ii]; X[3,ii]]
    end
  end
  if bvi>1
    dr=X[4,:]
    dq=X[5,:]
    rdq = zeros(size(X,2)); dr = zeros(size(X,2)); dq = zeros(size(X,2))
    for ii in 1:size(X,2)
      dr[ii] = X[4,ii]
      dq[ii] = X[5,ii]
      rdq[ii] = r[ii]*dq[ii]
      X[4:6,ii] = dcm*[dr[ii]*c[ii]-rdq[ii]*s[ii]; dr[ii]*s[ii]+rdq[ii]*c[ii];
                  X[6,ii]]
    end
  end
  if bvi==4
    for ii in 1:size(X,2)
      d2r=X[7,ii]
      d2q=X[8,ii]
      r_=d2r-rdq[ii]*dq[ii]
      q_=r[ii]*d2q+2.*dr[ii]*dq[ii]
      X[7:9,ii] = dcm*[r_*c[ii]-q_*s[ii]; r_*s[ii]+q_*c[ii]; X[9,ii]]
      r_=X[10,ii]-3.*dr[ii]*dq[ii]^2-3.*rdq[ii]*d2q
      q_=r[ii]*X[11,ii]+3.*d2r*dq[ii]+3.*dr[ii]*d2q-rdq[ii]*dq[ii]^2
      X[10:12,ii] = dcm*[r_*c[ii]-q_*s[ii];r_*s[ii]+q_*c[ii]; X[12,ii]]
    end
  end
end
#get data outside of spline range from Horizons
if any(ist)
  sav=param(dict,"save")
  param(dict,"save",false)
  if bvi==4
    warn("no accel or jerk from Horizons")
  end
  Xs=mkeph(bib,t,typ,dict)
  X[:,!ist]=X
  X[1:6,ist]=Xs
  param(dict,"save",sav)
end
return X,t
end

function RVspline6!{T}(X::AbstractArray{T},pp::Dict)
#5th order spline fiting R,V at segment endpoints, continuous A,J
#nonuniform step ok
t=X[1,:]
R=X[2:4,:]
V=X[5:7,:]*86400
X = 0
n=length(t)-1
R_ = Array(Float64,(3,n)); V_ = Array(Float64,(3,n)); dt = Array(Float64,n);
b_ = Array(Float64,(3,n))
for ii in 1:n
  dt[ii] = t[ii+1]-t[ii]
  for jj in 1:3
    R_[jj,ii] = (R[jj,ii]+V[jj,ii]*dt[ii]-R[jj,ii+1])/dt[ii]^2
    V_[jj,ii]=(V[jj,ii]-V[jj,ii+1])/dt[ii]/2.
  end
end
#R_+c2+c3.*dt+c4.*dt2,V_+c2+3/2*c3.*dt+2*c4.*dt2,
#c5(i)*dt^3+6*R_-6*V_+c2(i)=c2(i+1) ,
#     (3*c5*dt^3+8*R_-6*V_+2*c2)/dt=(c5*dt^3-4*R_+2*V_-2*c2)/dt
#c5(i)=(6*(V_-R_)-c2(i)+c2(i+1))/dt^3 ,
#      (12*V_-10*R_-c2(i-1)+3*c2(i))/dt=(8*V_-10*R_-3*c2(i)+c2(i+1))/dt
#(-15*R_+16*V_-2*c2(i-1)+3*c2(i))/dt^2=(15*R_-14*V_+3*c2(i)-2*c2(i+1))/dt^2
b = Array(Float64,(3,n-1))
s = Array(Float64,(3*(n-1)))
for ii in 2:n
  for jj in 1:3
    b[jj,ii-1]=(12.*V_[jj,ii-1]-10.*R_[jj,ii-1])/dt[ii-1]-
                (8.*V_[jj,ii]-10.*R_[jj,ii])/dt[ii]
  end
  s[ii-1]= 1./dt[ii-1]
  s[ii+n-2] = -3./dt[ii-1]-3./dt[ii]
  s[ii+2*n-3] = 1./dt[ii]
end
#s*[c2(i-1);c2(i);c2(i+1)]=b, tridiagonal
ii = 2:n
jj=[ii-1; ii; ii+1] #sparsity pattern x,y,z
ii=[ii; ii; ii]
b_ = zeros(3) #continuous 4th der at 2
for kk in 1:3
  b_[kk]=(-15.*R_[kk,1]+16.*V_[kk,1])/dt[1]^2-
          (15.*R_[kk,2]-14.*V_[kk,2])/dt[2]^2
end
s_=[2./dt[1]^2; -3./dt[1]^2+3./dt[2]^2; -2./dt[2]^2]
b=[b_ b]
s=[s_; s]

ii=[1; 1; 1; ii]
jj=[1; 2; 3; jj]
b_ = zeros(3) #continuous 4th der at n
for kk in 1:3
  b_[kk]=(-15.*R_[kk,n-1]+16*V_[kk,n-1])/dt[n-1]^2-
          (15.*R_[kk,n]-14.*V_[kk,n])/dt[n]^2
end
s_=[2./dt[n-1]^2; -3./dt[n-1]^2+3./dt[n]^2; -2./dt[n]^2]
b=[b b_]
s=[s; s_]
ii=[ii; n+1; n+1; n+1]
jj=[jj; n-1; n; n+1]

s=sparse(jj,ii,s) #solve for c2
c2=b/s
c22=c2[:,2:n+1]
c2=c2[:,1:n]
#solve for coefficients
c5 = Array(Float64,(3,n)); c3 = Array(Float64,(3,n));
c4 = Array(Float64,(3,n))
c0 = Array(Float64,(3,n)); c1 = Array(Float64,(3,n));
for ii in 1:n
  for jj in 1:3
    c5[jj,ii] = (6.*(V_[jj,ii]-R_[jj,ii])-c2[jj,ii]+c22[jj,ii])/dt[ii]^3
    c5dt3 = c5[jj,ii]*dt[ii]^3
    c3[jj,ii]=(2.*V_[jj,ii]-4.*R_[jj,ii]-2.*c2[jj,ii]+c5dt3)/dt[ii]
    c4[jj,ii]=(3.*R_[jj,ii]-2.*V_[jj,ii]+c2[jj,ii]-2.*c5dt3)/dt[ii]^2
    c0[jj,ii] = R[jj,ii]
    c1[jj,ii] = V[jj,ii]
  end
end
pp["breaks"] = t
pp["coefs1"] = [c5[:] c4[:] c3[:] c2[:] c1[:] c0[:]]
pp["pieces"] = length(t)-1
pp["order1"] = 6
pp["dim"] = 3
pp["dcm1"] = []
pp["dcm2"] = []
pp["dcm3"] = []
pp["dcm4"] = []
#derivatives
pp["coefs2"] = [5.*c5[:] 4.*c4[:] 3.*c3[:] 2.*c2[:] c1[:]]/86400.
pp["order2"] = 5
pp["coefs3"] = [20.*c5[:] 12.*c4[:] 6.*c3[:] 2.*c2[:]]/86400.^2
pp["order3"] = 4
pp["coefs4"] = [60.*c5[:] 24.*c4[:] 6.*c3[:]]/86400.^3
pp["order4"] = 3
return pp
end

###############################################################################
# Following function RVspline5 is not currently used
###############################################################################
# function pp=RVspline5(X)
# t=X(1,:);R=X(2:4,:);V=X(5:7,:)*86400;
# dt=repmat(diff(t),3,1);dt2=dt.*dt;
# c0=R(:,1:end-1);c1=V(:,1:end-1);
# R_=(c0+c1.*dt-R(:,2:end))./dt2;V_=(c1-V(:,2:end))./dt/2;
# #R_+c2+c3.*dt+c4.*dt2,V_+c2+3/2*c3.*dt+2*c4.*dt2,6*(R_-V_)+c2(i)=c2(i+1),
# c2=((V_(:,2)-2*R_(:,2))./dt(:,2)-(10*R_(:,1)-9*V_(:,1))./dt(:,1)).*dt(:,1)/2;
# RV6=6*(R_-V_);c2=cumsum([c2 RV6],2);c2=c2(:,1:end-1);
# #c2=((V_(:,end)-2*R_(:,end))./dt(:,end)-
# #(3*V_(:,end-1)-2*R_(:,end-1))./dt(:,end-1)).*dt(:,end-1)/2;
# #c2=cumsum([c2 -RV6(:,end-1:-1:1)],2);c2=c2(:,end:-1:1);#/2+c2_/2;
# c3=(V_-2*R_-c2)*2./dt;c4=(3*R_-2*V_+c2)./dt2;
# pp.form='pp';pp.breaks=t;pp.coefs=[c4(:) c3(:) c2(:) c1(:) c0(:)];
# pp.pieces=numel(t)-1;pp.order=5;pp.dim=3;pp.dcm=[];
# pp(2)=pp;pp(2).coefs=[4*c4(:) 3*c3(:) 2*c2(:) c1(:)]/86400;pp(2).order=4;
# return
###############################################################################

function RVAspline!{T}(X::AbstractArray{T},pp::Dict)
#5th order spline using R,V,A at beginning and end of segments
#approximates A with central diff of R and V
#assumes uniform time step
t=X[1,:]
R=X[2:4,:]
V=X[5:7,:]*86400
n = length(t)
dt = Array(Float64,n-1) #uniform timestep, scale V
for ii in 1:n-1
  dt[ii] = t[ii+1] - t[ii]
end
ct = 1
while (ct<n-2)
  if abs(dt[ct+1]-dt[ct])>1e-9
    error("ERROR:nonuniform timestep for spline")
    ct = n
  end
  ct = ct + 1
end
dt=dt[1]
V=V*dt

c2 = Array(Float64,(3,n-2))#central diff for A (c2 = coeff on t^2 term)
for ii in 1:n-2
  for jj in 1:3
    c2[jj,ii]=R[jj,ii]-2.*R[jj,ii+1]+R[jj,ii+2]+(V[jj,ii]-V[jj,ii+2])/4
  end
end
c21 = Array(Float64,3)#initial A (uses first 3 R,V)
c2n = Array(Float64,3)#final A (uses last 3 R,V)
for jj in 1:3
  c21[jj]=-23./4.*R[jj,1]+4.*R[jj,2]+7./4.*R[jj,3]-3.*V[jj,1]-
            4.*V[jj,2]-V[jj,3]/2.
  c2n[jj]=7./4.*R[jj,n-2]+4.*R[jj,n-1]-23./4.*R[jj,n]+V[jj,n-2]/2.+
            4.*V[jj,n-1]+3.*V[jj,n]
end
# coeff from initial R,V,A and temp final R,V,A on segs
c0 = Array(Float64,(3,n-1)); c1 = Array(Float64,(3,n-1))
R2 = Array(Float64,(3,n-1)); V2 = Array(Float64,(3,n-1))
#R_+c2+c3.*dt+c4.*dt2,V_+c2+3/2*c3.*dt+2*c4.*dt2,
#satisfy R and V at end of segment
#c5(i)*dt^3+6*R_-6*V_+c2(i)=c2(i+1) ,
#(3*c5*dt^3+8*R_-6*V_+2*c2)/dt=(c5*dt^3-4*R_+2*V_-2*c2)/dt
#continuous 2nd & 3rd derivative
#c5(i)=(6*(V_-R_)-c2(i)+c2(i+1))/dt^3 ,
#(12*V_-10*R_-c2(i-1)+3*c2(i))/dt=(8*V_-10*R_-3*c2(i)+c2(i+1))/dt
#combine into single eqn in c2
#(-15*R_+16*V_-2*c2(i-1)+3*c2(i))/dt^2=
#(15*R_-14*V_+3*c2(i)-2*c2(i+1))/dt^2
#continuous 4th derivative at end points

#solve for R2,V2,A2 at seg endpts
for ii in 1:n-1
  for jj in 1:3
    c0[jj,ii] = R[jj,ii]
    R2[jj,ii] = c0[jj,ii]-R[jj,ii+1]
    c1[jj,ii] = V[jj,ii]
    V2[jj,ii] = V[jj,ii+1]
  end
end
A2=[c2 c2n]
c2=[c21 c2]
c3=-10.*R2-6.*c1-4.*V2-3.*c2+A2
c4= 15.*R2+8.*c1+7.*V2+3.*c2-2.*A2
c5=-6.*R2-3.*c1-3.*V2-c2+A2
c=zeros(3*n-3,6)
# unscale time
for jj in 1:3*n-3
  c[jj,6]=c0[jj]
  c[jj,5]=c1[jj]/dt
  c[jj,4]=c2[jj]/dt^2
  c[jj,3]=c3[jj]/dt^3
  c[jj,2]=c4[jj]/dt^4
  c[jj,1]=c5[jj]/dt^5
end
pp["breaks"] = t
pp["pieces"] = n-1
pp["coefs1"] = c
pp["order1"] = 6
pp["dcm1"] = []
pp["dim"] = 3
d=[5 0 0 0 0; 0 4 0 0 0; 0 0 3 0 0; 0 0 0 2 0; 0 0 0 0 1.; 0 0 0 0 0]
pp["coefs2"] = c*d/86400
pp["order2"] = 5
pp["dcm2"] = []
return pp
end

function orient(bb::AbstractArray{Int64},tt::AbstractArray{Float64},typ::Bool,
                ll::Int64,dict::Dict)
# currently assumes that tt is a 1 x n or n x vector
#ll=1: Pole
#ll=2: [Node;Q;Pole] where Node=cross([0;0;1],Pole), Q=cross(Pole,Node)
#ll=3: Prime Meridian
#ll=4: [Node;Q;Pole] where Node=Prime Meridian
#ll=5: Equinox, cross(Pole,H)
#ll=6: [Node;Q;Pole] where Node=Equinox

sb=size(bb)
if (sb[1]>2) && (sb[2]>2)
  bb=bb[1,:]
elseif sb[1]>2
  bb=bb'
end

#get saved data, idori is list of body numbers, dori is spline data
id=param(dict,"idori")
d=param(dict,"dori")
(isempty(d)) && (d= Dict()) # d is a dictionary
(isempty(id)) && (id = zeros(Int64,(length(bb),3)))

tb=false #match t and b 1-to-1 if numel(t)==size(b,2)
if !typ
  tt=transpose(tt[:])'
  tb=length(tt)==length(bb)
end

# Find unique bodies and when they occur
it = collect(1:length(bb))
ct = collect(2:length(bb))
for ii in 1:length(bb)-1
  sp = 1
  T1 = identity(ct)
  for jj in ct
    if bb[ii] == bb[jj]
      splice!(T1,sp)
      it[jj] = ii
      sp = sp - 1
    end
    sp = sp + 1
  end
  ct = identity(T1)
  (isempty(ct)) && (continue)
end

if ll == 1 || ll == 3
  (!typ) ? (nn = (3,length(tt))) : (nn = (3,length(collect(tt[1]:tt[3]:tt[2]))))
elseif ll == 4 || ll == 2
  (!typ) ? (nn = (9,length(tt))) : (nn = (9,length(collect(tt[1]:tt[3]:tt[2]))))
elseif ll == 5
  nn = (3,1)
end

if tb
  t = Array(Float64,length(tt))
  xx =Array(Float64,nn)
else
  t = Array(Float64,size(tt))
  xx = Array(Float64,(nn[1],nn[2],length(bb)))
end

for ibu in 1:length(bb) #match times for each unique body
  if it[ibu] != ibu #match times for each unique body
    if tb && (tt[ibu]==tt[it[ibu]])
      t = tt[it[ibu]]
      X[:,ibu] = X[:,it[ibu]]
      continue
    elseif !tb
      t = tt
      X[:,:,ibu] = X[:,:,it[ibu]]
      continue
    end
  end
  bib=bb[ibu]
  tb? (t=tt[ibu]) : (t=tt)

  #small body orientation not in Horizons, may be in pck file. write xx 1-to-1 if tb
  if bib>1e3
    X=qorient(bib,t,ll,dict)
    if tb
      xx[:,ibu]=X
    else
      xx[:,:,ibu]=X
    end
    continue
  end

  #cl2 is Pole|PM|eqx output, ii is cl2 data exists, ij is any data exists
  cl2=ceil(Int,ll/2)
  iN=find(x->(x==bib),id)
  if !isempty(iN) # possibly coded wrong since not sure what id is like
    T1 = size(id)
    jj = zeros(Int,length(iN)); ij = []
    for kk in 1:length(iN)
      T2 = iN[kk]
      while T2 > 0
        ij = T2
        T2 = T2 - T1[1]
        jj[kk] = jj[kk] + 1
      end
    end
  else
    jj = []
    ij = []
  end
  any(x->(x==cl2),jj) ? (ii = identity(ij)) : (ii = [])
  ssd=param(dict,"ssd")
  if isempty(ii)||ssd#read from file or horizons
    if ssd #te is saved time span
      te=t
      (length(te)==3) ? (typ2 = true) : (typ2 = false)
    else
      te=[-.5; 60*365.25;.1]
      typ2 = true
    end

    #file names
    fls=["pol" "pm" "eqx"]; err=0
    @static is_windows()? (edir= string(param(dict,"bdir"),"\\ephem\\",fls[cl2],bib,".jld")):(edir=
                string(param(dict,"bdir"),"/ephem/",fls[cl2],bib,".jld"))
    rf=!ssd&&isfile(edir)
    if rf
      #read from file if rf, 1st two entries is data size
      di = load(edir,"di")
      #n=fread(f,2,'double')';di=fread(f,n(1:2),'double')
    else
      if cl2 == 1#read from horizons
        X,td,err=mkeph(bib,te,typ2,dict,[0 90])#pole
        sp = falses(size(X,2))
        for kk in 1:size(X,2)
          (X[3,kk]<0) && (sp[kk]=true)
        end
        if any(sp) #z component is positive for sqrt
          X[:,sp]=-X[:,sp]
          if any(!sp)
            warn("pole crosses ecliptic")
          end
        end
      elseif cl2 == 2
        X,td,err=mkeph(bib,te,typ2,dict,[0 0])
        P,_=orient([bib],te,typ2,1,dict)#prime meridian
      elseif cl2 == 3
        sav=param(dict,"save")
        param(dict,"save",0)
        X_,td=ephem1([bib],te,typ2,3,dict)
        param(dict,"save",sav)#get H
        P,_=orient([bib],te,typ2,1,dict)
        X = zeros(3,size(X_,2)) #equinox
        for kk in size(X_,2)
          X[:,kk] = cross(P[:,kk],cross(X[1:3,kk],X[4:6,kk]))
        end
      end
      for kk in 1:size(X,2) #unit vector
        XT=X[1:3,kk]
        Xd = sqrt(dot(XT,XT))
        for jj in 1:3
          X[jj,kk]=X[jj,kk]/Xd
        end
      end
      X = X[1:3,:]
      if cl2==1 #X is pole node, Pole x [0;0;1]
        di=[td;X[1:2,:]]
        if ll==2
          P=identity(X)
          for jj in 1:size(X,2)
            X[1,jj] = P[2,jj]
            X[2,jj] = -P[1,jj]
            X[3,jj] = 0
            XT = sqrt(dot(X[:,jj],X[:,jj]))
            for kk in 1:3
              X[kk,jj] = X[kk,jj]/XT
            end
          end
        end
        #if cl2==1;di=[td;X(1:2,:)];if ll==2;P=X;X=cross(qorient(399,0*td,1),P);
        #X=X./([1;1;1]*sqrt(sum(X.*X)));end#X is pole node, Pole x [0;0;1]
      #set up basis orthogonal to pole, N is pole node (Pole x [0;0;1]) Q is P x N
      else
        N = zeros(3,size(P,2)); Q = zeros(3,size(P,2)); w = zeros(size(P,2))
        for jj in 1:size(P,2)
          p12=sqrt(P[1,jj]^2+P[2,jj]^2)
          N[:,jj]=[P[2,jj]/p12; -P[1,jj]/p12; 0]
          Q[:,jj]=[-P[3,jj]*N[2,jj]; P[3,jj]*N[1,jj]; -p12]#Pole x [0;0;1]
          #else;N=cross(qorient(399,0*td,1),P);N=N./repmat(sqrt(sum(N.*N)),3,1);
          #Q=cross(P,N);#Earth pole x pole
          w[jj] = atan2(dot(X[:,jj],Q[:,jj]),dot(X[:,jj],N[:,jj]))
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
          w[kk:end] = w[kk:end] + 2*pi
        end
        di=[td;w']
      end#cl2
    end#rf
    if any(isnan(di[2,:])) #data to save, don't save nans
      jj = find(isnan(di[2,:]))
      for kk in 1:length(jj)
        di = [di[:,1:jj[kk]-1] di[:,1:jj[kk]+1]]
        jj = jj -1
      end
    end
    if isempty(di)
      (typ) && (t=collect(t[1]:t[3]:t[2]))
      X=NaN*ones(3+6*(cl2==div(ll,2)),length(t))#no data, output nans
      (tb) ? (xx[:,ibu]=X) : (xx[:,:,ibu]=X)
      continue
    end

    #see if times are out of horizons range if ssd
    if ssd
      (length(err)==1) && (err = [err false])
      if err[1]
        warn("Min Horizons ephem time for ",bib," is EXCEEDED")
        #,datestr(err[1]+730486.5))
      elseif err[2]
        warn("Max Horizons ephem time for ", bib," is EXCEEDED")
        #,datestr(err[2]+730486.5))
      end
    #save data if no error and didn't read from existing file
    elseif (!rf) && (param(dict,"save")) && all(!err)
      fid = jldopen(edir,"w")
      fid["di"] = di
      close(fid)
    end

    if !ssd#make spline if not ssd flag
      if id==zeros(Int64,1,3)
        ii = 1
      elseif !isempty(ij)
        ii=ij[1]
      else
        ii=size(id,1)+1
      end

      if isempty(ii) || ii > size(id,1)
        T1 = zeros(Int,1,3); T1[cl2] = bib
        id = [id; T1]
      else
        id[ii,cl2]=bib
      end

      #see if body already exists in idori
      (!haskey(d,ii)) && (d[ii] = Dict())
      d[ii][cl2] = Dict()
      spline!(di[1,:],di[2:end,:],d[ii][cl2])
      if param(dict,"keeps")
        param(dict,"idori",id)
        param(dict,"dori",d)
      end#save data
    end
  end#ii|ssd

  #convert {start,end,delta} to list
  if typ
    t=collect(t[1]:t[3]:t[2])
  end
  ts=[]
  ist=falses(length(t))
  if !ssd
    #see if data is out of bounds
    tl=d[ii][cl2]["breaks"][1]
    if minimum(t)<tl
      kk = find(x->(x<tl),t)
      for jj in kk
        ist[jj] = true
      end
      warn("Min time for ",bib," ephem spline is ",Dates.julian2datetime(tl+2451545))
    end
    tl=d[ii][cl2]["breaks"][end]
    if maximum(t)>tl
      kk = find(x->(x>tl),t)
      for jj in kk
        ist[jj] = true
      end
      warn("Max time for ",bib," ephem spline is ",Dates.julian2datetime(tl+2451545))
    end
    if any(ist) #ist is indeces of data outside of spline range
      ts=t[ist]
      t[ist]=[]
    end
    if (id[ii,1]==0)
      P,_=orient([bib],t,typ,1,dict)
    else
      T1 = zeros(2,length(t))
      pval!(d[ii][1],t,3,T1)
      P = zeros(3,length(t)); P[1:2,:] = T1
    end
    #get new pole data if not already in spline
    X = zeros(3,length(t)) #pole vector
    p12 = zeros(length(t))
    for kk in 1:length(t)
      p12[kk] = P[1,kk]^2 + P[2,kk]^2
      P[3,kk] = sqrt(1-p12[kk])
      for jj in 1:3
        X[jj,kk] = P[jj,kk]
      end
    end

    if ll>1
      for jj in 1:length(t) #X is pole node: Pole x [0;0;1]
        p12[jj] = sqrt(p12[jj])
        X[1,jj] = P[2,jj]/p12[jj]
        X[2,jj] = -P[1,jj]/p12[jj]
        X[3,jj] = 0
      end

      if cl2>1 #get angle and rotate
        w = zeros(1,length(t))
        pval!(d[ii][cl2],t,3,w)
        for jj in 1:length(t)
          XT = X[:,jj]
          X[1,jj] = cos(w[jj])*XT[1]+sin(w[jj])*-P[3,jj]*XT[2]
          X[2,jj]=cos(w[jj])*XT[2]+sin(w[jj])*P[3,jj]*XT[1]
          X[3,jj]=cos(w[jj])*XT[3]+sin(w[jj])*-p12[jj]
        end
        #if ll>1;X=cross(qorient(399,0*t,1),P);
        #X=X./repmat(sqrt(sum(X.*X)),3,1);
        #X is pole node: Earth pole x Pole
        #if cl2>1;w=ppval(d(ii,cl2),t);
        #X=[1;1;1]*cos(w).*X+[1;1;1]*sin(w).*cross(P,X);
        #get angle and rotate
      end
    end
  end#ssd

  if ll==2 && bb!=399
    X = qorient([399],0*ones(t),typ,1,dict)
    for jj in 1:length(t)
      XL = sqrt(dot(X[:,jj],X[:,jj]))
      for kk in 1:3
        X[kk,jj]=X[kk,jj]/XL
      end
    end
  end
  if cl2==div(ll,2) #dcm
    XT = identity(X); X = zeros(9,size(X,2))
    for jj in 1:size(X,2)
      XP = cross(P[:,jj],XT[:,jj])
      for kk in 1:3
        X[kk,jj] = XT[kk,jj]
        X[kk+3,jj] = XP[kk]
        X[kk+6,jj] = P[kk,jj]
      end
    end
  end

  if any(ist) #get horizons data for out of range times
    sav=param(dict,"save")
    param(dict,"save",false); param(dict,"ssd",1);#reset flags
    Xs=orient([bib],ts,false,ll,dict)
    X[:,!ist]=X
    X[:,ist]=Xs
    param(dict,"save",sav)
    param(dict,"ssd",ssd)
  end

  if tb   #match 1-to-one, or all times per body
    xx[:,ibu]=X
  else
    xx[:,:,ibu]=X
  end

end#for ibu

return xx, t
end

function spline!{P}(x::AbstractArray{P},y::AbstractArray{P},pp::Dict)
# Cubic spline with x1 <= x2 <= .. <= xn
# Assumes dx is constant
# x is a 1 x n or n x 1 array
# y is a 2 X n, n X 2, 1 x n, n x 1 array

(size(x,2)>size(x,1)) && (x = x')
dx = x[2]-x[1]

# Find y difference
(size(y,2)>size(y,1)) && (y = y')
dy = diff(y)
n = length(x)

# Create sparse tridiagonal sparse matrix
iv = zeros(Int,3*n-2); jv = zeros(Int,3*n-2); av = zeros(Int,3*n-2);
# row coordinates
iv[1] = 1; iv[2] = 1; iv[3*n-2] = n; iv[3*n-3] = n
for ii in 1:n-2
  for jj in 1:3
    iv[3*ii+jj-1] = ii + 1
  end
end
# column coordinates
jv[1] = 1; jv[2] = 2; jv[3*n-2] = n; jv[3*n-3] = n-1
for ii in 1:n-2
  for jj in 1:3
    jv[3*ii+jj-1] = ii + jj - 1
  end
end
# values in tridiagonal
av[1] = 1; av[2] = 2; av[3*n-2] = 1; av[3*n-3] = 2
for ii in 1:n-2
  av[3*ii] = 1
  av[3*ii+1] = 4
  av[3*ii+2] = 1
end
A = sparse(iv,jv,av)*dx

# Create vector so that A D = C
if size(y,2) == 2
  C1 = zeros(n); C2 = zeros(n)
  C1[1] = (5*dy[1,1]+dy[2,1])/2; C1[end] = (5*dy[n-1,1]+dy[n-2,1])/2
  C2[1] = (5*dy[1,2]+dy[2,2])/2; C2[end] = (5*dy[n-1,2]+dy[n-2,2])/2
  for ii in 2:n-1
    C1[ii] = 3*(y[ii+1,1]-y[ii-1,1])
    C2[ii] = 3*(y[ii+1,2]-y[ii-1,2])
  end

  D1 = \(A,C1)
  D2 = \(A,C2)

  # Determine the coefficients of the third order polynomial
  coefs = zeros(2*(n-1),4)
  for jj in 1:2*(n-1)
    coefs[jj,4] = y[jj]
  end

  for jj in 1:n-1
    coefs[2*jj-1,3] = D1[jj]
    coefs[2*jj,3] = D2[jj]
    coefs[2*jj-1,2] = 3*dy[jj,1]/dx^2 - 2*D1[jj]/dx - D1[jj+1]/dx
    coefs[2*jj,2] = 3*dy[jj,2]/dx^2 - 2*D2[jj]/dx - D2[jj+1]/dx
    coefs[2*jj-1,1] = -2*dy[jj,1]/dx^3 + D1[jj]/dx^2 + D1[jj+1]/dx^2
    coefs[2*jj,1] = -2*dy[jj,2]/dx^3 + D2[jj]/dx^2 + D2[jj+1]/dx^2
  end
elseif size(y,2) == 1
  C1 = zeros(n)
  C1[1] = (5*dy[1,1]+dy[2,1])/2; C1[end] = (5*dy[n-1,1]+dy[n-2,1])/2
  for ii in 2:n-1
    C1[ii] = 3*(y[ii+1,1]-y[ii-1,1])
  end
  D1 = \(A,C1)

  # Determine the coefficients of the third order polynomial
  coefs = zeros(n-1,4)
  for jj in 1:n-1
    coefs[jj,4] = y[jj]
    coefs[jj,3] = D1[jj]
    coefs[jj,2] = 3*dy[jj,1]/dx^2 - 2*D1[jj]/dx - D1[jj+1]/dx
    coefs[jj,1] = -2*dy[jj,1]/dx^3 + D1[jj]/dx^2 + D1[jj+1]/dx^2
  end
end

# send out the results
pp["breaks"] = x
pp["pieces"] = n-1
pp["coefs3"] = coefs
pp["order3"] = 4
pp["dim"] = size(y,1)
end

function mkeph{P}(n::Int64,t::AbstractArray{Float64},typ::Bool,dict::Dict,
                ll::AbstractArray{P}=[])
#read from horizons
# typ is true if t is a range of times
err= [false; false]
tephf=getef(dict) #get ephemeris time span
if (n>10) && (n<1e3) && (mod(n,100)!=99) #get data wrt central body
  cb=floor(n/100)*100+99
else
  cb=10
end
if (n>1e6) && (n<2e6) #comets are dumb
  CAP="#3BCAP"
else
  CAP=""
end
if n>1e6 #time limits for small bodies
  tlim=tephf.sb
  DES="DES="
  ieph=1
else #time limits for bog bodies
  ieph=find(x->(x==string(n)),tephf.numbers)
  tlim=tephf.tl[ieph[1],:]
  DES=""
end

#tt is expected time output from horizons, can only run 400 times in list,
# or 90,000 times in {start,end,delta}
if !typ
  ist = sortperm(collect(t)) #dev add collect to force column vector
  tt = t[ist]
  mnt=400
else
  tt=collect(t[1]:t[3]:t[2])
  mnt=90000
  (length(tt)<2) && (warn("invalid {start time, end time, delta time} input to Horizons"))
end
if (isempty(ieph))
  warn("no horizons ephemeris file for ",n)
  X=NaN*ones(6,numel(tt))
  err = [true]
  return X,tt,err
end
if minimum(tt) < tlim[1]
  ii = find(x->(x<tlim[1]),tt)
  nt1=length(ii)
  err[1] = true
elseif minimum(tt) >= tlim[1]
  nt1 = 0
  err[1] = false
end
if maximum(tt) > tlim[2]
  jj = find(x->(x>tlim[2],tt))
  nt2=length(jj)
  err[2] = true
else
  nt2 = 0
  err[2] = false
end
n2t=length(tt)-nt2
t12=round(Int64,linspace(nt1,n2t,ceil(Int64,(n2t-nt1)/mnt)+1))
#t12 breaks time span into different runs
if n>1e6
  tnow=indmin(abs(Dates.datetime2julian(now())-2451545.-tt))
  t12=sort(unique([tnow-1; t12]))
end

X=zeros(7,0)
n12=length(t12)
println("Retrieving data for Body ",n," from Horizons")
for ti in 1:n12-1#only do mnt at a time
  #time input for url
  if typ
    tstr = @sprintf("&START_TIME='JD%%20%.9f'&STOP_TIME='JD%%20%.9f'&STEP_SIZE='%d'",
    tt[t12[ti]+1]+2451545,tt[t12[ti+1]]+2451545,t12[ti+1]-t12[ti]-1)
  else
    tstr = "&TLIST='"
    for ii in t12[ti]+1:t12[ti+1]
      tstr = tstr*@sprintf("%.9f%%0A",tt[ii]+51544.5)
    end
    tstr = tstr*"'"
  end
  #set target as center of body and observer as ll (lat,lon) for orientation
  if !isempty(ll)
    cb=n
    coord="coord"
    llstr=@sprintf("&COORD_TYPE='GEODETIC'&SITE_COORD='%.8f,%.8f,0'",ll[2],ll[2])
  else
    coord=""
    llstr=""
  end
  # horizons URL
  url=@sprintf("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='%s%d%s'",
      DES,n,CAP)
  url = url*@sprintf("&CENTER='%s@%d%s'&MAKE_EPHEM='YES'",coord,cb,llstr)
  url = url*@sprintf("&TABLE_TYPE='VECTORS'%s&OUT_UNITS='KM-S'",tstr)
  url = url*"&VECT_TABLE='2'&REF_PLANE='ECLIPTIC'&REF_SYSTEM='J2000'"
  #url = url*"&OBJ_DATA='NOPE'" # Originally included but errors out
  #url = url*"&VECT_CORR='NONE'&VEC_LABELS='NO'&CSV_FORMAT='NO'"
  s = download(url) #run horizons url command and extract goodies
  f = open(s); rd = readlines(f); close(f)
  mx = length(rd); i1 = 0; T1 = true
  while T1
    i1 = i1 + 1
    (contains(rd[i1],"\$\$SOE")) && (T1 = false)
    (i1 == mx) && (T1 = false)
  end
  T1 = true; i2 = identity(i1)
  (i2==mx) && (T1=false)
  while T1
    i2 = i2 + 1
    (contains(rd[i2],"\$\$EOE")) && (T1 = false)
    (i2 == mx) && (T1 = false)
  end
  clw = length(rd[1])
  (i1!=mx) && (rd = 0)
  if (i1!=mx)
    X1 = readdlm(s,skipstart=i1,use_mmap=true,skipblanks=true,
                    dims=(mx,clw))
    X1 = X1[1:i2-i1-1,1:3]
    XT = zeros(7)
    for ii in 1:floor(Int64,size(X1,1)/3)
      XT[1] = X1[3*ii-2,1]
      XT[2] = X1[3*ii-1,1]
      XT[3] = X1[3*ii-1,2]
      XT[4] = X1[3*ii-1,3]
      XT[5] = X1[3*ii,1]
      XT[6] = X1[3*ii,2]
      XT[7] = X1[3*ii,3]
      X = [X XT]
    end
  #successful run
  #tri=regexp(ss,'Center radii\s*: ([\d.]+) x ([\d.]+) x ([\d.]+) ','tokens');
  #putsbmb(cb,{'tri'},{str2num(cat(1,tri{1}{:}))});#triaxial also in horizons
  else
    println(rd)
    warn("Horizons error")
    X=NaN*ones(6,length(tt))
    return X, tt, err
  end
  if ti!=n12-1
    println("Horizons progress: ",floor(Int64,(t12[ti+1]-nt1)/(n2t-nt1)*100),"%")
  end
end#ti
if nt1 == 0 && nt2 == 0
  tt=transpose(X[1,:]-2451545.)
elseif nt1 == 0
  tt=[X[1,:]-2451545. transpose(tt[end-nt2+1:end])]
elseif nt2 ==0
  tt=[transpose(tt[1:nt1]) X[1,:]-2451545.]
else
  tt=[transpose(tt[1:nt1]) X[1,:]-2451545. transpose(tt[end-nt2+1:end])]
end


X=[NaN*ones(6,nt1) (1-2*(!isempty(ll)))*X[2:7,:] NaN*ones(6,nt2)]
#pad out of range data with nans, if orient X=-X
if !typ
  X[:,ist]=X
  tt[ist]=tt
end#output times in same order as input

@static is_windows()? (edir= string(param(dict,"bdir"),"\\ephem")):(edir=
            string(param(dict,"bdir"),"/ephem"))

if !isdir(edir)
  mkdir(edir)
end#make "ephem" directory
if (param(dict,"save")) && (isempty(ll)) #save ephemeris data
  @static is_windows()? (f= edir*"\\"*string(n)*".jld"):(f= edir*"/"*string(n)*".jld")
  fid = jldopen(f,"w")
  d=[tt;X]
  if !typ
    d=d[:,ist]
  end
  if sum(isnan(X[1,:]))!=0 #save times in order, don't save nans
    rmc = []
    for ii in 1:size(X,2)
      (!isnan(X[1,ii])) && (rmc = [rmc; ii])
    end
    d = slicedim(d,2,rmc)
  end
  # write(fid,"sz",size(d)) # don't need size
  write(fid,"d",d)
  close(fid)  #write data
  if n>1e6 #always update small bodies
    n,e=getsb(n,dict)
    putsbmb(n[1],["ephref" "ephdate"],
            [e Dates.datetime2julian(now())-2451545.],dict)
  else #update reference file and time
    putsbmb(n,["ephref" "ephdate"],
            [tephf.f[ieph,:] Dates.datetime2julian(now())-2451545.],dict)
  end
end # save
return X,tt,err
end

function getephdat(b::Int64,dict::Dict)
#ephemeris header data
ef1=getef(dict)
ii = find(x->(x==string(b)),ef1.numbers)
ef = lowercase(match(r"^[^-._ ]*",ef1.f[ii[1]]).match)
#only need first file if merged
ef = convert(String,ef)
if isempty(ef)
  warn("no ephemeris header file for ",b)
  ed=[]
  return ed
end
ed=param(dict,ef) #see if file data already saved
if (isempty(ed)) && (b<400) #inner planets and barycenters
  #read in DE### header data
  ef=match(r"de\d+",ef).match #deXXX format
  ef = convert(String,ef)
  efdir="ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/"*ef*"/" #directory
  s = download(efdir)
  f = open(s); rd = readstring(f); close(f)
  T1 = true; fils = []; off = 1; ct = 0; #sometimes in "header.XXX_YYY" format
  T2 = Regex(join(["header.",ef[3:end],"_?\\d*"],""))
  fil = matchall(T2,rd[off:end])
  fils = union(fil)

  if isempty(fils) #standard filename is "header.XX"
    T2 = Regex(join(["header.",ef[3:end]],""))
    fil = match(T2,rd[off:end])
    if fil == nothing
      warn("Ephemeris header file ",ef, " for ", b, "not found")
      return ed
    else
      fils = fil.match
      s = download(efdir*fils)
    end
  else
    s = download(efdir*fils[1])
  end
  f = open(s); rd = readstring(f); close(f)
  ii= search(rd,"GROUP   1040")
  jj = search(rd,"GROUP   1041")
  kk = search(rd,"GROUP   1050")

  T1 = split(rd[ii[end]+1:jj[1]-1])  #read in variable names (each 6 chars)
  vars = Array(AbstractString,length(T1)-1)
  for ii in 1:length(T1)-1
    vars[ii] = identity(T1[ii+1])
  end

  T1 = split(rd[jj[end]+1:kk[1]-1])  #read in variable names (each 6 chars)
  vals = zeros(length(T1)-1)
  for ii in 1:length(T1)-1
    TD = search(T1[ii+1],"D")
    T11 = parse(Float64,T1[ii+1][1:TD[1]-1])
    T1E = parse(Float64,T1[ii+1][TD[1]+1:end])
    vals[ii] = T11*10.^T1E
  end
  ed = Dict()
  ed["vars"]=vars
  ed["vals"]=vals
  #write gm values
  au=ed["vals"][uniqstr(ed["vars"],"AU",1)]; au = au[1]
  ed["gm10"]=ed["vals"][uniqstr(ed["vars"],"GMS")]*au^3/86400.^2 #sun
  ed["gm3"]=ed["vals"][uniqstr(ed["vars"],"GMB")]*au^3/86400.^2 #E-M bary
  emr=ed["vals"][uniqstr(ed["vars"],"EMRAT")]; emr = emr[1]#E/M ratio
  ed["gm301"]=ed["gm3"]/(1.+emr);ed["gm399"]=ed["gm3"]*emr/(1.+emr);
  #Moon and Earth barycenters
  ed["gm0"]=ed["gm10"]+ed["gm3"]
  for ii in [1:2;4:9]
    GM="gm"*string(ii)
    gm=ed["vals"][uniqstr(ed["vars"],"GM"*string(ii))]*au^3/86400^2
    ed[GM] = gm
    ed["gm0"]=ed["gm0"]+gm
  end
  ed["gm199"]=ed["gm1"];ed["gm299"]=ed["gm2"];#Mercury and Venus are Trouble
  #radii and J2
  ed["rad10"]=ed["vals"][uniqstr(ed["vars"],"ASUN")]
  ed["j2_10"]=ed["vals"][uniqstr(ed["vars"],"J2SUN")]
  ed["rad399"]=ed["vals"][uniqstr(ed["vars"],"RE")]
  ed["j2_399"]=ed["vals"][uniqstr(ed["vars"],"J2E",1)]
  ed["rad301"]=ed["vals"][uniqstr(ed["vars"],"AM",1)]
  ed["j2_301"]=ed["vals"][uniqstr(ed["vars"],"J2M")]
  ed["rad199"]=ed["vals"][uniqstr(ed["vars"],"RAD1")]
  ed["rad299"]=ed["vals"][uniqstr(ed["vars"],"RAD2")]
  param(dict,ef,ed)#save data
elseif isempty(ed) #outer planets and satellites
  #check default directory and eph filename
  s =  download("ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/")
  f = open(s); rd = readstring(f); close(f)
  if rd[1]!='<'
    fil = matchall(r"\S+(?=.txt)",rd[1:end])
    fils = union(fil)
  else
    fil = matchall(r">\S+(?=.txt)",rd[1:end])
    fils = union(fil)
    for ii in 1:length(fils)
      fils[ii] = strip(fils[ii],'>')
    end
  end
  ii = uniqstr(fils,ef,0)
  if isempty(ii)
    jj = find(x->(contains(x,ef)),fils)
    (isempty(jj)) && (ii = uniqstr(fils,ef))
    (isempty(ii))? (err = false) : (err = true)
    if isempty(ii) # Look for the first 3 leters of the highest # file
      fils2 = Array(typeof(fils[1]),size(fils))
      nf = length(fils2)
      for kk in 1:nf
        fils2[kk] = fils[nf-kk+1]
      end
      ii = uniqstr(fils2,ef[1:3])
      ii = nf - ii + 1 #
      (isempty(ii))? (err = false) : (err = true)
    end
  #could also check
  #ftp://ssd.jpl.nasa.gov/pub/eph/satellites/rckin/rckin.*.log,
  # but GM likely 0
  else
    err = true
  end

  if err
    s =  download("ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/"
          *fils[ii[1]]*".txt")
    f = open(s); rd = readstring(f); close(f)
    #s is ephemeris header file
    #get GM from "Bodies on the File" table: skip space,keep 3 #s,
    # skip some space, keep some #s, decimal, & non-space,skip space
    T1 = true; fils = [];
    fil = matchall(r"\s(\d{3})\s+((\d+\.\S*))\s",rd)
    fils = union(fil)
    gm = Array(AbstractString,(length(fils),2))
    for ii in 1:length(fils)
      gm[ii,:] = split(fils[ii])
    end

    ed = Dict()
    for ii in 1:size(gm,1)
      ed["gm"*gm[ii,1]] = parse(Float64,gm[ii,2])
    end
    #get J2 and radius of planet
    T1 = search(rd,"J"*string(floor(Int64,b/100))*"02")
    T2 = search(rd[T1[end]:end],"J")
    ed["j2_"*string(convert(Int64,b))] =
                      parse(Float64,strip(rd[T1[end]+1:T2[1]+T1[end]-2]))
    T1 = search(rd,"RADIUS")
    T2 = search(rd[T1[end]:end],"J")
    ed["rad"*string(convert(Int64,b))] =
                    parse(Float64,strip(rd[T1[end]+1:T2[1]+T1[end]-2]))
    #end#if
  else
    warn("ephemeris header file ",ef, " for ",b," not found")
  end#err
  param(dict,ef,ed)
end
return ed
end

function getef(dict::Dict)
#Ephemeris file and time spans used by Horizons
tephf=param(dict,"tephf")
if isempty(tephf.numbers)
  s = download("https://ssd.jpl.nasa.gov/eph_spans.cgi?id=A") #Planets
  f = open(s); rd = readstring(f); close(f)
  #Read in bodynumber, begin time, " not " or " to " flag, end time, and file
  ssT = Array(AbstractString,(1,5))
  off = 1; T1 = true; sss = []
  while T1
    ss=match(r"<td.*?>(\d+)</td>.*?<\w\w?>(.*?) (not|to) (.*?)</\w\w?>&nbsp;.*?&nbsp;(.*?)&nbsp;",rd[off:end])
    if ss == nothing
      T1 = false
    else
      off = ss.offsets[end]+off
      for ii in 1:5
        ssT[ii] = ss.captures[ii]
      end
      if isempty(sss)
        sss = identity(ssT)
      else
        sss = [sss; ssT]
      end
    end
  end
  #Mercury and Venus are Trouble (no ephemeris file, but still in table)
  jj = find(x->(x=="not"),sss[:,3])
  for ii in jj
    warn("Potential that Mercury and Venus information not correct")
    kk=find(x->(x==sss[ii,1]*"99"),sss[:,1])
    if !isempty(kk)
      sss[ii,:]=sss[kk[1],:]
      sss[ii,1]=string(sss[kk[1],1][1])
    end
  end
  sss = sss[:,[1:2;4:5]]
  s = download("https://ssd.jpl.nasa.gov/eph_spans.cgi?id=B") #Satellites
  f = open(s); rd = readstring(f); close(f)
  #Read in bodynumber, begin time, end time, and file (no need to flag " not "
  #or " to "
  ssT = Array(AbstractString,(1,4))
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
  tephf.tl =
    [Dates.datetime2julian(DateTime(sss[:,2])) Dates.datetime2julian(DateTime(sss[:,3]))]-2451545.

  s = download("https://ssd.jpl.nasa.gov/eph_spans.cgi?id=D") #Satellites
  f = open(s); rd = readstring(f); close(f)
  ssT = Array(AbstractString,(1,2))
  off = 1; T1 = true; ss1=[]
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
        ss1 = identity(ssT)
      else
        ss1 = [ss1; ssT]
      end
    end
  end
  tephf.sb =
    [Dates.datetime2julian(DateTime(ss1[:,1])) Dates.datetime2julian(DateTime(ss1[:,2]))]-2451545.
  param(dict,"tephf",tephf)
end
return tephf
end

function getsatdat(dict::Dict)
#Satellite data page, has GM and mean radius (& density, magnitude, albedo)
sd=param(dict,"satdat")
if isempty(sd)
  s=download("https://ssd.jpl.nasa.gov/?sat_phys_par")
  f = open(s); sd = readstring(f); close(f)
  T1 = search(sd,"Earth's Moon")
  if (T1 != 0:-1)
    sd = sd[T1[1]:end]
    T2 = search(sd,"References")
    (T2 != 0:-1) && (sd = sd[1:T2[1]])
  end
  param(dict,"satdat",sd)
end
return sd
end

function getpck(b::Int64,v::String,dict::Dict)
#NAIF pck file, analytic orientation & contains some small body radii and
#orientation data not on Horizons
pck=param(dict,"pck")
if isempty(pck)
  s = download("ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/")
  f = open(s); rd = readstring(f); close(f)
  pck=match(r"pck\d+.tpc",rd).match
  #get latest file in directory
  s=download("ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/"*pck)
  #read contents to string
  f = open(s); rd = readstring(f); close(f)
  rd = replace(rd,"\n",";"); rd = replace(rd,"(","["); rd = replace(rd,")","]")
  #turn newlines into ; and parens into brackets for matlabese
  rd = replace(rd,"2431010","2000243")
  pck = replace(rd,"9511010","2000951") #change # of Ida and Gaspra to SPK-ID
  param(dict,"pck",pck)
end#save
v1 ="BODY"*string(floor(Int64,b))*"_"*uppercase(v) #variable name
T1 = Regex(join([v1,".*?]"],""))
if match(T1,pck) == nothing
  m1 = []
else
  m1 = match(T1,pck).match#just single param
end
#d=regexp(pck,['\\begindata[;\s]+(' v{1} '.*?;)[;\s]+\\begintext'],'tokens');
#everything in data block
d = []
if !isempty(m1)
  m1 = split(m1)
  for ii in 1:length(m1)
    (m1[ii]=="[") && (continue)
    (m1[ii]=="]") && (continue)
    (m1[ii]=="=") && (continue)
    if isa(parse(m1[ii]),Number)
      (!isempty(d)) && (d = [d; parse(Float64,m1[ii])])
      (isempty(d)) && (d = parse(Float64,m1[ii]))
    end
  end
end
return d
end

function getx(b1,s::String,dict::Dict)
#get data from mb or sb structure
x,_=getsbmb(b1,dict)
ii=find(y->(y==b1),x.numbers) #get datastructure and see if body is there
(length(ii)==1) && (ii = ii[1])
if (!isempty(ii))&&(!param(dict,"ssd"))&&
    (!isempty(find(y->(y==Symbol(s)),fieldnames(x))))
  #see if should skip saved data and if datafield exists
  #see if entry exists in field for body
  if isempty(getfield(x,(Symbol(s))))
    o = []
  elseif (typeof(getfield(x,Symbol(s))[1])==Float64)||(typeof(getfield(x,Symbol(s))[1])==Int64)
    if length(getfield(x,Symbol(s)))>=ii
      (!isnan(getfield(x,Symbol(s))[ii])) && (o=getfield(x,Symbol(s))[ii])
      (isnan(getfield(x,Symbol(s))[ii])) && (o=[])
    end
  elseif typeof(getfield(x,Symbol(s)))==Array{AbstractString,1}
    if length(getfield(x,Symbol(s)))>=ii
      (!isempty(getfield(x,Symbol(s))[ii])) && (o = getfield(x,Symbol(s))[ii])
    end
  else
    o=[]
  end
else
  o=[]
end
return o
end

function getmb(dict::Dict)
#major body list
nn=download("https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=MB")
f = open(nn); lines = readlines(f); close(f)

NumNam = []
for ii in 1:length(lines) #keep 1--3 numbers, skip spaces, keep stuff until two spaces
  split1 = split(lines[ii],"  ",keep=false) # skip spaces
  if length(split1) > 1 # remove extra spaces
    split1[1] = strip(split1[1])
    split1[2] = strip(split1[2])
  end
  if isnumber(split1[1])
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

jj = find(x->(x>9),NumNam[:,1]) #match barycenters last
nums = zeros(size(NumNam,1))
nams = Array(AbstractString,size(NumNam,1))
ct = 0
for ii in jj
  nums[ct+1] = NumNam[ii,1]
  nams[ct+1] = NumNam[ii,2]
  ct = ct + 1
end
jj = find(x->(x<=9),NumNam[:,1])
for ii in jj
  nums[ii+ct] = NumNam[ii,1]
  nams[ii+ct] = NumNam[ii,2]
end

#b = vcat(nums[1:10],nums[15:20],nums[22:23],nums[24:50],nums[160:end])
x,_=getsbmb("mb",dict)
b=identity(x.numbers)
jj = [] # see if body numbers already exist in data structure
for ii in 1:length(b)
  T1 = find(x->(x==b[ii]),nums)
  if isempty(jj)
    jj = T1
  else
    jj = vcat(jj,T1)
  end
end

T1 = Array(AbstractString,length(jj))
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

function getsb(n,dict::Dict)
#small body web pages
if typeof(n) == String #name or number
  srch = escape(n)
else
  srch=string(convert(Int64,n))
end
s = download(string("https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=",srch))
f = open(s); rd = readstring(f); close(f)
ss = search(rd,r"\+1\"><b>[^<]+") #Name is bigger font "+1"
if ss != 0:-1
  nam = rd[ss[8]:ss[end]]
  ss1 = search(rd,r">\d{7}<")
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
  ss1 = search(rd,r">GM<.*?>")
  if ss1 != 0:-1
    ss2 = search(rd[ss1[end]:end],r">([\d.e-]{1,})</")
    a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
    (isempty(a)) && (a=0); (isempty(a)) || (a=parse(Float64,a))
  else
    a = 0
  end
  nnf = [nnf a]
  ss1 = search(rd,r">diameter<.*?>")
  if ss1 != 0:-1
    ss2 = search(rd[ss1[end]:end],r">([\d.]{1,})</")
    a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
    (isempty(a)) && (a=0); (isempty(a)) || (a=parse(Float64,a)/2.)
  else
    a = 0
  end
  nnf = [nnf a]
  ss1 = search(rd,r">rot_per<.*?>")
  if ss1 != 0:-1
    ss2 = search(rd[ss1[end]:end],r">([\d.]{1,})</")
    a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
    (isempty(a)) && (a=0); (isempty(a)) || (a=parse(Float64,a)/24.)
  else
    a = 0
  end
  nnf = [nnf a]
  ss1 = search(rd,r">H<.*?>")
  if ss1 != 0:-1
    ss2 = search(rd[ss1[end]:end],r">([\d.]{1,})</")
    a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
    (isempty(a)) && (a=NaN); (isempty(a)) || (a=parse(Float64,a))
  else
    a=NaN
  end
  nnf = [nnf a]
  ss1 = search(rd,r">condition code<.*?sp;(\d)&nb")
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
  ss1 = search(rd,r">spec_T<.*?>")
  if ss1 != 0:-1
    ss2 = search(rd[ss1[end]:end],r">([\w])</")
    if ss2!= 0:-1 #dev
        a = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]
    else
        a = [];
    end
    (isempty(a)) && (a="")
  else
    a=""
  end
  ss1 = search(rd,r">spec_B<.*?>")
  if ss1 != 0:-1
    ss2 = search(rd[ss1[end]:end],r">(\w+)</")
    (ss1 != 0:-1) ? (b = rd[ss1[end]+ss2[1]:ss1[end]+ss2[end]-3]) : (b="_")
    (!isempty(a)) && (b="/"*b)
  else
    b="_"
  end
  nnf = [nnf a*b]
  ss1 = search(rd,r">Reference: <.*?>.*?>.*?>(.*?)<")
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
  warn("multiple matches might be possible in boddat.",
  " possibly not coded correctly")
  ss2 = search(rd,r">([^<]+)</a></td>")
  nnf = [NaN]; eph = ""
end
return nnf, eph, s
end

function uniqsb(b,dict::Dict)
bs=[b,string(b,"*"),string("*",b,"*")] #search exact, beginning, fragment
for b1 in bs
  nn,_,_ = getsb(b1,dict) #check ssd
  if isnan(nn[1])&&!isempty(nn[2]);
    warn(b1,"returned multiple matches. not coded for multiple matches")
    #n is {number name provisional} of firt match, keep first of those
    #n=regexp(nn{2}(1,:),'(\d*\s?)([^\(]*)(.*)','tokens');
    #n=n{1}{find(~cellfun('isempty',n{1}),1)};
    #nn=getsb(deblank(n));if isnan(nn{1});warning('#s didn''t work',n);end
    #get data from webpage
  end#if
  if !isnan(nn[1])
    return nn
  end#found one!
end#for
return nn
end

function getsbmb(n,dict::Dict)
#get small or major body structure
if typeof(n)!=String #should be 'mb' or 'sb'
  if n<1e6; bi="mb"; else; bi="sb"; end#convert number input to sb mb
else
  bi=identity(n)
end
bo=param(dict,bi)
if (bo == []) && (bi == "mb")#read from saved data
  mb = MB([],[],[],[],[],[],[],[],[],[])
  param(dict,"mb",mb)
elseif (bo == []) && (bi == "sb")#read from saved data
  sb = SB([],[],[],[],[],[],[],[],[],[],[])
  param(dict,"sb",sb)
elseif (isempty(bo.Names)) || (isempty(bo.numbers))
  @static is_windows()? (bdf= string(param(dict,"bdir"),"\\bodydata.jld")):(bdf=
              string(param(dict,"bdir"),"/bodydata.jld"))
  if isfile(bdf)
    boddict = load(bdf,"boddict") # this will give "md" or "sd" elements
    bo = boddict[bi]
    param(dict,bi,bo)
  end
end#save
return bo,bi
end

function putsbmb(n,f,v,dict::Dict)
#write data to structure
x,bi=getsbmb(n[1],dict) #get structure
bb=find(x->(x==n[1]),x.numbers)
if isempty(bb)&&length(n)>1
  n=n[2]
  bb=find(x->(x==n),x.numbers)
else
  n=n[1]
end
if isempty(bb) #see if body exists in field
  bb=length(x.numbers)+1
  f=["numbers" f]
  (size(v,2) > size(v,1))? (v = [n v]) : (v = [n; v])
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
      setfield!(x, Nm, Array(AbstractString,length(x.numbers)))
    end
    for jj in bb
     getfield(x, Nm)[jj]= v[ii]
    end
  else
    #write, fill with nans
    if isempty(getfield(x, Nm))
      if Nm == :numbers
        setfield!(x, Nm,zeros(Int64,length(v[ii])))
      else
        (isempty(x.numbers)) && (setfield!(x, Nm, NaN*zeros(length(v[ii]))))
        (!isempty(x.numbers))&&(setfield!(x, Nm, NaN*zeros(length(v[ii]),length(x.numbers))))
      end
    end
    nn=size(getfield(x,Nm),2)+1
    ct = 1
    for jj in bb
      if typeof(getfield(x,Nm)[:,jj][1])==typeof(v[ii])
        getfield(x,Nm)[:,jj]=v[ii]
      else
        getfield(x,Nm)[:,jj]=convert(typeof(getfield(x, Nm)[:,jj][1]),v[ii])
      end
      if bb[ct]>nn
        getfield(x,Nm)[:,nn:bb[ct]-1]=NaN
      end
      ct = ct + 1
    end
  end
end #numel(f)
 (param(dict,bi*"mod")!= true) && (param(dict,bi*"mod",true))
param(dict,bi,x)#save
return x
end

function uniqstr(xni,bi,nb=[])
#find a unique match of bi in xni
(nb==[]) && (nb=0)
#1 for only whole word matches, 2 for only fragment matches, 0 for either
if typeof(xni) == String
  T1 = false
  (contains(xni,bi)) && (T1=true)
else
  T1 = falses(length(xni))
  for ii in 1:length(xni)
    (contains(xni[ii],bi)) && (T1[ii]=true)
  end
end
In = find(x->(x==true),T1)
#in=strfind(xn,[' ' bi ' ']);#whole word

if (length(In) > 1) && (nb == 1) #multiple entries with same word
  #In is strfind index, xm tracks current list of multiple matches
  xm = Array(AbstractString,length(In))
  for ii in 1:length(In)
    xm[ii] = xni[In[ii]]
  end
  In = In[find(x->(x==bi),xm)] #exact match
  if length(In) > 1
    warn("identical entries")
    In = In[1]
  end
  (isempty(In)) && (warn("ambiguous string match"))
elseif (length(In) > 1)  && (nb!=1)
    if any(x->(x==xni[In[1]]),xni[In[2:end]])
      warn("ambiguous string match")
    elseif any(x->(x[1:3]==xni[In[1]][1:3]),xni[In[2:end]])
      warn("ambiguous string match")
    end
    In=In[1]
#if numel(in)>1
#  in1=in;xm=xni(ceil(in/nxn),:);in=in(in==floor(in/nxn)*nxn+1);
#beginning of entry
#  if isempty(in);in=in(regexp(xn(in-1),'\W'));#beginning of word
#      if numel(in)>1;xm=xni(ceil(in/nxn),:);elseif isempty(in);in=in1;end
#  elseif numel(in)>1;xm=xni(ceil(in/nxn),:); end
# in1=in;xm=xni(ceil(in/nxn),:);in=in(isspace(xn(in-1)));#beginning of word
# if numel(in)>1
#     in1=in;xm=xni(ceil(in/nxn),:);in=in(in==floor(in/nxn)*nxn+2);
#beginning of entry
#     if numel(in)>1;xm=xni(ceil(in/nxn),:);elseif isempty(in);in=in1;end
# elseif isempty(in);in=in1; end
      #if numel(in)>1;matches=xm,warning('ambiguous string match');in=in(1);end
end
return In
end

# function param(p::AbstractString)
# #variables global to function, works with only single input
# #saves p in local space and write v to p, then retrieve p from any function via o
# fid = jldopen("vs.jld", "r+")
# if exists(fid,p) # recalls v from saved space
#   o=read(fid,p)
# else # Output empty matrix
#   o=[]
# end
# close(fid)
# return o
# end
#
# function param{T}(p::AbstractString,v::T)
# #variables global to function, works with only single input
# #saves p in local space and write v to p,
# #  then retrieve p from any function via o
# fid = jldopen("vs.jld", "r+")
# # sets v to saved space and overwrites value if already there
# if exists(fid, p)
#   nt = length(names(fid))
#   jj = find(x->(x==p),names(fid))
#   fid2 = jldopen("vs2.jld","w")
#   for ii in 1:nt
#     if ii != jj[1]
#       fid2[names(fid)[ii]] = read(fid[names(fid)[ii]])
#     end
#   end
#   close(fid); close(fid2)
#   rm("vs.jld"); mv("vs2.jld","vs.jld")
#   fid = jldopen("vs.jld","r+")
# end
# fid[p] = v
# o = []
# close(fid)
# return o
# end

function param(vs::Dict,key::String,value)
# Allows a variety of values to be stored in RAM under one name
# Sets value to dictionary structure
vs[key] = value
end

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

function idxhist!{T}(xs::AbstractArray{T},rg::AbstractArray{Float64},
        index::AbstractArray{Int64},lx::Int64,l::Int64)
sr = 1
for jj in 1:lx
  while (index[jj] == 0)
    if (xs[jj] >= rg[sr]) && (xs[jj]<rg[sr+1])
      index[jj] = sr
      (sr > l) && (index[jj] = l)
    end
    sr = sr + 1
    (sr > l ) && (continue)
  end
  (sr!=1) && (sr = sr - 1)
  (sr > l ) && (continue)
end

end

function pval!{T}(pp::Dict,xx::AbstractArray{T},order::Int64,v::AbstractArray{T})
# Assumes that xx is sorted from smallest to largest
# pp is dictionary
# xx is a 1 x n or n x 1
# output dimension is always 3
xs = vec(xx)
lx = length(xs)
b = pp["breaks"]
c = pp[string("coefs",order)]
l = pp["pieces"]
k = pp[string("order",order)]
(size(v,2)==lx)? (vs = size(v,1)) : (vs = size(v,2))

index = zeros(Int64,lx)
rg = [-Inf;b[2:end];Inf]

idxhist!(xs,rg,index,lx,l)

ii = find(xs==Inf)
(!isempty(ii)) && (index[ii] = l)
badin = find(index==0) # Remove NaNs
(!isempty(badin)) && (index[badin] = 1)

xs = xs-b[index]

index = vs*index
T1 = zeros(vs*lx); T2 = zeros(Int64,vs*lx)
for ii in 1:lx
  for jj in -(vs-1):0
    T1[vs*ii+jj] = xs[ii]
    T2[vs*ii+jj] = index[ii]+jj
  end
end
xs = T1
index = T2

T1 = zeros(vs*lx)
for ii in 1:vs*lx
  T1[ii] = c[index[ii],1]
end

for jj in 2:k
  for ii in 1:vs*lx
    T1[ii] = xs[ii]*T1[ii] + c[index[ii],jj]
  end
end

for ii in 1:lx
  for jj in -(vs-1):0
    v[jj+vs,ii] = T1[vs*ii+jj]
  end
end

end


function p3val!{T}(pp::Dict,xx::AbstractArray{T},order::Int64,v::AbstractArray{T})
# Assumes that xx is sorted from smallest to largest
# pp is dictionary
# output dimension is always 3
xs = vec(xx)
lx = length(xs)
b = pp["breaks"]
c = pp[string("coefs",order)]
l = pp["pieces"]
k = pp[string("order",order)]

index = zeros(Int64,lx)
rg = [-Inf;b[2:end];Inf]

idxhist!(xs,rg,index,lx,l)

ii = find(xs==Inf)
(!isempty(ii)) && (index[ii] = l)
badin = find(index==0) # Remove NaNs
(!isempty(badin)) && (index[badin] = 1)

xs = xs-b[index]

index = 3*index
T1 = zeros(3*lx); T2 = zeros(Int64,3*lx)
for ii in 1:lx
  T1[3*ii-2] = xs[ii]
  T1[3*ii-1] = xs[ii]
  T1[3*ii] = xs[ii]
  T2[3*ii-2] = index[ii]-2
  T2[3*ii-1] = index[ii]-1
  T2[3*ii] = index[ii]
end
xs = T1
index = T2

T1 = zeros(3*lx)
for ii in 1:3*lx
  T1[ii] = c[index[ii],1]
end

for jj in 2:k
  for ii in 1:3*lx
    T1[ii] = xs[ii]*T1[ii] + c[index[ii],jj]
  end
end

for ii in 1:lx
  v[1,ii] = T1[3*ii-2]
  v[2,ii] = T1[3*ii-1]
  v[3,ii] = T1[3*ii]
end

end

function getfn(fn::String,b::AbstractArray{Int64},dict::Dict)
  a = []
  if (fn == "gm") || (fn=="GM")
    a = getgm(b,dict)
  end
  if fn == "radius"
    a = getrad(b,dict)
  end
  if (fn=="j2") || (fn=="J2")
    a = getj2(b,dict)
  end
  if (fn=="rot") || (fn=="rotation_period") || (fn=="ROT")
    a = getrot(b,dict)
  end
  if (fn=="occ")
    a = getocc(b,dict)
  end
  if (fn=="typ") || (fn=="types")
    a = gettyp(b,dict)
  end
  if (fn=="ref") || (fn=="reference")
    a = getref(b,dict)
  end
  if (fn=="clo") || (fn=="close") || (fn=="close_approach")
    a = getclo(b,dict)
  end
  if (fn=="dat") || (fn=="date")
    a = getdat(b,dict)
  end
  if (fn=="mag") || (fn=="magnitude")
    a = getmag(b,dict)
  end
return a
end

type MB
  Names::Array{AbstractString}
  numbers::AbstractArray{Int64}
  gm::AbstractArray{Float64}
  rad::AbstractArray{Float64}
  j2::AbstractArray{Float64}
  rot::AbstractArray{Float64}
  ephref::Array{AbstractString}
  ephdate::AbstractArray{Float64}
  pp::AbstractArray{Float64}
  tri::AbstractArray{Float64}
end

type SB
  Names::Array{AbstractString}
  numbers::AbstractArray{Int64}
  gm::AbstractArray{Float64}
  rad::AbstractArray{Float64}
  rot::AbstractArray{Float64}
  h::AbstractArray{Float64}
  occ::AbstractArray{Float64}
  types::Array{AbstractString}
  ephref::Array{AbstractString}
  ephdate::AbstractArray{Float64}
  pp::AbstractArray{Float64}
end

type TEPHF
  numbers::Array{AbstractString}
  f::Array{AbstractString}
  tl::AbstractArray{Float64}
  sb::AbstractArray{Float64}
end
