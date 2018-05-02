function varargout=boddat(bvars,bi,t,opti)
%Body Data from Horizons and ssd.jpl.nasa.gov
%
%   varargout=boddat(parameters,bodies,times,options)
%
%"parameters" is a cell array or comma, semicolon, or space delimited string that
%   designates which parameters to output in varargout.
%The available parameters are:
%   Number (NAIF or SPK ID, e.g. Sun is 10, Earth is 399, Moon is 301, Ceres is 2000001)
%   Name (body name)
%   CB (central body of bodies)
%   GM (gravitational paramter, km^3/s^2)
%   radius (mean radius, km)
%   triaxial (radii of a triaxial ellipsoid shape model, km)
%   J2 (gravitational harmonic due to oblateness)
%   rotation_period (period of rotation, days)
%Small-body-specific parameters
%   mag (absolute magnitude [ H ], magnitude at 1 AU from Sun and observer)
%   occ (orbital condition code 0-9, higher numbers indicate uncertain orbits)
%   type (spectral type, Tholen/SMASII taxonomic classification)
%   close_approach (table of close approach time, body, and AU distance)
%Ephemeris
%   R (position vector, km)
%   V (velocity vector, km/s)
%   X (position and velocity)
%   A (position, velocity, acceleration, & jerk)
%   RTN (rotating frame direction cosine matrix unit[R ; RxVxR ; RxV])
%   reference (orbit solution ID or ephemeris file, used to track changes) 
%   date (Creation date of ephemeris datafile in "ephem" directory, days past J2000.
%         The latest generation date from Horizons is not available.) 
%Orientation
%   Pole (spin axis of body, unit vector)
%   PM (prime meridian of body, unit vector)
%   Eqx (Equinox direction, Pole x H, unit vector)
%
%   If "dcm" follows Pole, PM, or Eqx then the direction cosine matrix is output 
%   as a 9x1 array, with Pole as the Z direction, and PM or Eqx as X.
%   The X direction for the Pole dcm is Earth pole x body Pole, unit vector.
%   If "q" precedes Pole, PM or Eqx, then a constant pole and spin rate is assumed
%   for quicker computations.
%All vectors are output in ecliptic and mean equinox of J2000 reference frame.
%
%"bodies" is a numerical or cell row vector, or a comma or semicolon delimited string
%   of body numbers or names. 
%   Small-body designations are designated with a character string of digits, 
%   e.g '3' returns Juno, while [3] returns Earth-Moon Barycenter.
%   The second row is reserved for the "with respect to" body for ephemeris calls.
%   If no second row is input or is NaN, then the "wrt" body is the central body
%
%"times" is days past high noon 1/1/2000 TDB (J2000) for ephemeris and orientation.
%   times can be a numerical array or 3-element cell of {start time, end time, delta time}.
%   If the number of time elements equals the number of body entries
%   then times and bodies are mapped 1-to-1, otherwise data is computed at every time
%   and concatenated by body along dim 3.
%
%"options" is a 3-element logical array for [save ssd keep_spline] where
%save = 1 saves ephemeris and orientation data into the (boddat.m path)/ephem/ directory
%   and other data in (boddat.m path)/bodydata.mat.
%ssd = 1 checks Horizons and SSD for data first, while ssd = 0 checks for data in 
%   memory, /ephem/ or bodydata.mat before resorting to Horizons. The default time span for
%   saved ephemeris and orientation data is 1/1/2000-1/1/2060 (unless span available on
%   Horizons is shorter). Planetary ephemeris is saved in 0.5-day intervals, while moon and 
%   orientation data are saved in 0.1-day intervals.
%keep_spline = 1 keeps the ephemeris and orientation splines in memory for subsequent calls.
%   Note: "clear boddat" will clear all the persistent variables in boddat.
%options is optional; the default is [1 0 1].
%
%Example: Ephemeris +/- one week of Apophis close approach in 0.1-day increments
%t = datenum('4/13/2029 22:00:00')-datenum('1/1/2000 12:00:00')+[-7 : 0.1 : 7];
%[Pos PM St n] = boddat('R,qpm;type num',{'Apophis' , 'Earth'},t);
%
%   (It can take up to a minute for Horizons to compute the ephemerides.)
%   Pos = the position of the asteroid Apophis and Earth with respect to the sun
%   (their central body). Pos(:,:,1) is Apophis and Pos(:,:,2) is Earth.
%   PM = the prime meridian of Apophis and Earth. Currently no orientation data is available
%   for Apophis, so PM(:,:,1)=NaN. The "q" in "qpm" flags a constant pole and spin rate for 
%   decreased computation time (and accuracy).
%   St = the spectral type of Apophis. St{2} is empty because Earth doesn't have a type.
%   n = the body numbers of Apophis and Earth
%
%Example (con'd): Using saved data
%RV = boddat('X',{'Apophis' ; 'Earth'},t);
%
%   (The run time should be a few orders of magnitude faster.)
%   RV = the state of Apophis with respect to Earth computed from saved polynomials.
%   Note that the "bodies" input has elements in the second row (transposed from above example)
%   indicating 'Earth' as the wrt body, so RV(1:3,:) = Pos(:,:,1)-Pos(:,:,2).
%
%Example (con'd): Bypassing saved data
%RV_h = boddat('X',[2099942 ; 399],{t(1),t(end),0.1},[0 1]);
%
%   RV_h = the state of Apophis with respect to Earth directly from Horizons.
%   The zero in the options input indicates to not overwrite any saved data, and
%   the one in the options input indicates to skip saved data and read directly from Horizons.
%   When using Horizons it is typically faster (but not necessary) to input times as 
%   {start, end, delta}, rather than a list of discrete times.

%SOURCES
%small body data https://ssd.jpl.nasa.gov/sbdb.cgi
%Horizons ephemeris and orientation https://ssd.jpl.nasa.gov/horizons_batch.cgi
%list of major bodies https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=MB
%Ephemeris file names and time spans https://ssd.jpl.nasa.gov/eph_spans.cgi
%Inner planet, sun & barycenter GM, planet radius J2 ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/
%Outer planet & satellite GM, planet radius J2 ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/
%Satellite GM and mean radius https://ssd.jpl.nasa.gov/?sat_phys_par
%analytic orientation and radii ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/


%peek into params (may be useful for troubleshooting)
if nargin==1;for ii=1:numel(bvars);varargout{ii}=params(bvars{ii});return;end;end


%if there's no time, the time is now
if nargin<3||isempty(t);t=now-730486.5;end
%default options is [1 0 1], replace with any input and save in params
opts=[1 0 1];if nargin<4;opti=[];end;opts(1:numel(opti))=opti;
params('save',opts(1));params('ssd',opts(2));params('keeps',opts(3));
%boddat.m directory, tells where to save data
if isempty(params('bdir'));bdir=mfilename('fullpath');ii=regexp(bdir,'/|\');params('bdir',bdir(1:ii(end)));end
%parameter outputs, save in params as character array
if isempty(params('bva'));ov={'Pole','PM','Eqx'};bva=char([
{'R','V','X','A','RTN'},...5
{'number','cb','triaxial','name','gm','radius','j2','rotation_period'},...13
{'reference','date'},...15
{'mag','occ','type','close_approach'},...19
ov,strcat(ov,'dcm'),strcat('q',[ov,strcat(ov,'dcm')])]);%31 (12 orientation variables)
%ov,strcat('dcm',ov),strcat('q',[ov,strcat('dcm',ov)])]);%31 (12 orientation variables)
params('bva',bva);
end
params('mbmod',false);params('sbmod',false);

%get body numbers as unique identifiers, ephemeris (be) can be 2 x n, otherwise first row of b is used
if ~isempty(bi);if isnumeric(bi);be=bi;else;[be bi]=getnum(bi);end;
    b=be;%if size(be,1)>2;be=be';end;b=be(1,:);
else varargout(1:nargout)={[]};return;end
%input body parameter variables, if not a cell array then cellstr or split by comma, semicolon, space & or |
if ~iscell(bvars);if size(bvars,1)>1;bvars=cellstr(bvars);else;bvars=regexp(bvars,'[,;\s&|]','split');bvars(cellfun('isempty',bvars))=[];end;end
%if strcmpi(bvars{1},'all');bvars=cellstr(params('bva'));end
%for each parameter input call the appropriate funciton and write result to varargout
for bvr=1:numel(bvars);bvar=bvars{bvr};
bvi=[];if numel(bvar)==1;bvi=strfind('RVXA',bvar);end%skip uniqstr for ephem
if isempty(bvi);bvi=uniqstr(params('bva'),bvar);end%find out which member of bva matches input
if isempty(bvi);warning('no match found for %s',bvar);varargout(bvr)={nan*ones(size(t))};continue;end
%match bvi
if bvi<5;%ephemeris, {R V X A}
    if ~params('ssd')|~exist('Xeph','var');Xeph=ephem(be,t,bvi);end
        a=Xeph;if bvi==1;a=a(1:3,:,:);elseif bvi==2;a=a(4:6,:,:);end
elseif bvi==5;Xeph=ephem(be,t,3);R=Xeph(1:3,:,:);H=cross(R,Xeph(4:6,:,:));a=cat(1,unit(R),unit(cross(H,R)),unit(H));
elseif bvi==6;a=getnumout(b);%number
elseif bvi==7;cb=cbfun(b);%central body, matches bi for name vs number output
    if isnumeric(bi);a=cb;else;a=bi;for ib=1:numel(bi);if isnumeric(bi{ib});a{ib}=cb(ib);else;a(ib)=getnam(cb(ib));end;end;end
%8--19 calls a function called getxxx, where xxx is first three letters of matching bva parameter. only input is b
elseif bvi==8;a=zeros(3,numel(b));for ub=unique(b(:))';iu=ub==b;a(:,iu)=repmat(gettri(ub),1,sum(iu(:)));end%triaxial 3 x n
elseif bvi<20;a=zeros(size(b));if any(bvi==[9 14 18 19]);a=cell(size(b));end%initialize, cell for name, ref, type, clo_app
        fn=params('bva');fn=fn(bvi,1:min(3,end));
        for ub=unique(b(:)');eval(['ua=get' fn '(ub);']);a(ub==b)=ua;end

% ii=~imag(b);
% if any(ii);fn=params('bva');fn=fn(bvi,1:min(3,end));eval(['a(:,ii)=get' fn '(b(ii));']);end
% if any(~ii);if iscell(a);a(~ii)={''};else;a(:,~ii)=nan;end;end
%20--31 is orientation
%ll=1: Pole
%ll=2: [Node;Q;Pole] where Node=cross([0;0;1],Pole), Q=cross(Pole,Node)
%ll=3: Prime Meridian
%ll=4: [Node;Q;Pole] where Node=Prime Meridian
%ll=5: Equinox, cross(Pole,H)
%ll=6: [Node;Q;Pole] where Node=Equinox
else fn=params('bva');fn=deblank(fn(bvi,:));
    ll=2*find(~cellfun('isempty',regexp(fn,{'Pole','PM','Eqx'})));
    if isempty(strfind(fn,'dcm'));ll=ll-1;end%check for dcm
    if fn(1)=='q';a=qorient(b,t,ll);else;a=orient(b,t,ll);end%check for quick calculation flag
end%if bvi<5
varargout(bvr)={a};%write output to varargout
end
%save data to bodydata
if params('save');
    fn=[params('bdir') 'bodydata.mat'];if ~exist(fn,'file');mb.names='';mb.numbers=[];sb=mb;save(fn,'sb','mb');end
    if params('mbmod');mb=getsbmb('mb');save(fn,'mb','-append');end
    if params('sbmod');sb=getsbmb('sb');save(fn,'sb','-append');end
end
return

function [b bi]=getnum(bi)
%get body number, only called if bi isn't numerical
mb=params('mb');%persistent mb
%change into string, bi is used in central body output
if ~iscell(bi);if size(bi,1)>1;bi=cellstr(bi);else;bi=regexp(bi,'[,&|;]','split');end;end;
b=zeros(size(bi));%initialize b to same size as bi
for ib=1:numel(bi)
if isnumeric(bi{ib});b(ib)=bi{ib};continue;else;bib=upper(strtrim(bi{ib}));end
if strfind(bib,'_CP');nbcp=1;bib=strrep(bib,'_CP','');else;nbcp=0;end
if isnumeric(bib);error('getnum b(ib)=bib;%already a number');end%delete
x=getsbmb('mb');if isempty(mb)&&(isempty(x.names)||params('ssd'));x=getmb;end;mb=x;%get mb from saved data unless empty or ssd flagged

L=regexp(bib,'L\d','match');%Lagrange point
if ~isempty(L);L=str2double(L{1}(2));
bib=regexp(regexprep(bib,'L\d',''),'[_-\s]','split');bib(cellfun('isempty',bib))=[];
n=getnum(bib);nb=numel(bib);if any(isnan(n))||nb<1||nb>2||L>5;error('%s is an ambiguous Lagrange point',bi{ib});end
if nb==1;n(2)=cbfun(n);end;if n(1)==cbfun(n(2));n=n([2 1]);elseif n(2)~=cbfun(n(1));error('%s is an ambiguous Lagrange point',bi{ib});end
b(ib)=n(1)*10+L;if ~any(b(ib)==x.numbers);putsbmb(b(ib),'numbers',b(ib));bib=getnam(n);putsbmb(b(ib),'names',sprintf('%s-%s L%d',bib{[2 1]},L));end

else%match a string
    in=uniqstr(x.names,bib,1);%check whole word in major body data
    if isempty(in);x=getsbmb('sb');sb=x;in=uniqstr(x.names,bib,1);%check whole word in small body data
    if isempty(in)||params('ssd');b_=uniqsb(bib);x=getsbmb('sb');sb=x;in=uniqstr(x.names,bib,1);%check ssd webpages
    if isempty(in);in=uniqstr(x.names,b_{2},1);end;end
        if isempty(in);fprintf('Warning: No whole-word match found for "%s"\n',bib);
            x=mb;in=uniqstr(x.names,bib,2);%check fragment in mb
            if isempty(in);x=sb;in=uniqstr(x.names,bib,2);%check fragment in sb
                if isempty(in);x=getmb;in=uniqstr(x.names,bib);end%mb potentially skipped
            end%sb any
        if ~isempty(in);fprintf('> Found it! "%s" set to "%s"\n\n',bib,x.names(in,:));end
        end%sb word
    end%mb any
if isempty(in);warning(sprintf('No match found for "%s"',bib));b(ib)=nan;else b(ib)=x.numbers(in);end
end%mb word
if nbcp;b(ib)=b(ib)+1i;end
end%for ib
return

function bo=getnumout(bi)
for ib=1:numel(bi);
b=bi(ib);if imag(b);nbcp=1;b=real(b);else;nbcp=0;end
if isnan(b);n=b;else
n=getx(b,'numbers');%check saved data (& no ssd flag)
if ~isempty(n)%found it
elseif b>1e3&&b<1e4&&mod(b,10)<6;n=b;%Lagrange point
elseif b<1e6;x=getmb;n=find(b==x.numbers);%check list of major bodies
    if ~isempty(n);n=x.numbers(n);else;n=nan;end
else;n=getsb(b);n=n{1};%check ssd
end
end%isnan
if isnan(n);warning('No match found for "%d"',b);end%nomatch
if nbcp;n=n+1i;end
bo(ib)=n;
end
return

function nn=getnam(bs)
%get body name
%write name that matches each number input to cell
nn=cell(size(bs));
for ib=1:numel(bs)
b=bs(ib);if imag(b);nbcp=1;b=real(b);else;nbcp=0;end
n=getx(b,'names');%check saved data (& no ssd flag)
if ~isempty(n)%found it
elseif b>1e3&&b<1e4&&mod(b,10)<6;L=mod(b,10);%Lagrange point
    b=floor(b/10);n={};n=[cbfun(b) b];n=getnam(n);
    n=sprintf('%s-%s L%d',n{:},L);putsbmb(10*b+L,'names',n);putsbmb(10*b+L,'numbers',10*b+L);
elseif b<1e6;x=getmb;n=find(b==x.numbers);%check list of major bodies
    if ~isempty(n);n=deblank(x.names(n,:));end
else;n=getsb(b);%check ssd
    if ~isnan(n{1});n=n{2};else;n=[];end
end
if isempty(n);n='NULL';warning('No match found for "%d"',b);end%nomatch
%format small body name, tokens in regexp are {number, name, (YYYY A1)}
if b>1e6;n=regexp(n,'(\d+\s)*([^\(]*)(.*)','tokens');n=n{1};%n=deblank(n{1});
    %if named, go with name, if numbered go with "# (provisional)", otherwise output provisional designation
    if ~isempty(n{2});n=deblank(n{2});elseif ~isempty(n{1});n=[n{:}];else;n=n{3}(2:end-1);end
end
if nbcp;n=[n '_CP'];end
nn{ib}=n;
end%for ib
return

function gmx=getgm(bs)
%GM from ephemeris header constants unless =0, then estimate from density and size
gmx=zeros(size(bs));
for jj=1:numel(bs)
b=bs(jj);if imag(b);gm=nan;else;gm=getx(b,'gm');end%check saved data (& no ssd flag)
if ~isempty(gm);%found it
elseif b>1e3&&b<1e4&&mod(b,10)<6;gm=nan;putsbmb(b,'gm',gm);%Lagrange point
elseif b<1e6;ed=getephdat(b);%major body, check ephemeris header for GM > 1e-9
    if ~isempty(ed)&&ed.(['gm' num2str(b)])>1e-9;gm=ed.(['gm' num2str(b)]);putsbmb(b,'gm',gm);
    %read general satellite page, get first # after >name < and ">"
    else;ed=getsatdat;n=getnam(b);gm=regexpi(ed,['>' n{1} '[\s<].*?>([\d.])+[\s&<]'],'tokens');
        if ~isempty(gm);gm=str2double(gm{1}{1});putsbmb(b,'gm',gm);
        else;gm=nan;warning('no GM found for %s',n{1});end %#ok<*WNTAG>
    end
else;gm=getsb(b);%small body
    if ~isnan(gm{1});gm=gm{3};else;gm=nan;warning('body %d not found for GM data',b);end
end
gmx(jj)=gm;
end%jj
return

function j2x=getj2(bs)
j2x=zeros(size(bs));
for jj=1:numel(bs);b=bs(jj);
if any(b==[10 [3:9]*100+99 301])%J2 currently available for only a few bodies from ephemeris header
if imag(b);j2=nan;else;j2=getx(b,'j2');end
if isempty(j2);ed=getephdat(b);j2=ed.(['j2_' num2str(b)]);putsbmb(b,'j2',j2);end
else;J2_bodies=[10 [3:9]*100+99 301],warning('J2 is only available for those bodies');j2=nan;
end
j2x(jj)=j2;end
return

function rx=getrad(bs)
rx=zeros(size(bs));
for jj=1:numel(bs)
b=bs(jj);if imag(b);r=nan;else;r=getx(b,'rad');end%check saved data (& no ssd flag)
if ~isempty(r)%found it
elseif b>1e3&&b<1e4&&mod(b,10)<6;r=nan;putsbmb(b,'rad',r);%Lagrange point
elseif b<1e6&mod(b,100)==99|b==10;ed=getephdat(b);r=ed.(['rad' num2str(b)]);
if isempty(r);r=gettri(b);r=r(1);end;putsbmb(b,'rad',r);%planet or sun, ephemeris header
            %satellite, read general satellite page, get second # after >name < and ">"
elseif b<1e6;ed=getsatdat;n=getnam(b);r=regexp(ed,['>' n{1} '[\s<].*?>[\d.]+[\s&<].*?>([\d.]+)[\s&<]'],'tokens');
    if ~isempty(r);r=str2double(r{1}{1});putsbmb(b,'rad',r);
    else;r=nan;warning('no radius found for %s',n{1});end
else;r=getsb(b);%small body
    if ~isnan(r{1});r=r{4};else;r=nan;warning('body %d not found for radius data');end
end
rx(jj)=r;
end%for 
return

function trix=gettri(bs)
trix=zeros([3 size(bs)]);
for jj=1:numel(bs)
b=bs(jj);if imag(b);r=nan*[1 1 1];else;r=getx(b,'tri');end%check saved data (& no ssd flag)
if isempty(r);r=getpck(b,'RADII');%check pck file
    if isempty(r)&&b>1e6%small body
    %extent is flag for triaxial radii,convert to #, make axisymmetric if only 2#s
    [nn ee ss]=getsb(b);r=regexp(ss,'>extent<.*?>(([\d.x\s]+)+[\d.]+)<','tokens');
    if ~isempty(r);r=r{1}{1};r=str2num(strrep(r,'x',' '));if numel(r)==2;r=[r(1) r];end;end
    end
    if isempty(r);r=getrad(b)*[1 1 1];end%see if there's a single radius value
    if ~isnan(r(1));putsbmb(b,'tri',r);end
end%if
trix(:,jj)=r;
end%for jj
return

function rx=getrot(bs)
rx=zeros(size(bs));
for jj=1:numel(bs)
b=bs(jj);if imag(b);r=nan;else;r=getx(b,'rot');end%check saved data (& no ssd flag)
if isempty(r);r=getpck(b,'PM');%check pck file
    if ~isempty(r);r=360/r(2);putsbmb(b,'rot',r);%convert to days
    elseif b>1e6;r=getsb(b);if ~isnan(r{1});r=r{5};else;r=nan;end%small body
    else;r=nan;end
    if isnan(r);warning('body %d not found for rotation data',b);end
end%if
rx(jj)=r;
end%for jj
return

function hx=getmag(bs)
hx=zeros(size(bs));
for jj=1:numel(bs);b=bs(jj);
if b<1e6||imag(b);hx(jj)=nan;continue;end%only for small bodies
h=getx(b,'h');%check saved data (& no ssd flag)
if isempty(h);h=getsb(b);%check ssd
if ~isnan(h{1});h=h{6};else;h=nan;warning('body %d not found for absolute magnitude',b);end
end%if
hx(jj)=h;
end%for jj
return

function hx=getocc(bs)
hx=zeros(size(bs));
for jj=1:numel(bs);b=bs(jj);
if b<1e6||imag(b);hx(jj)=nan;continue;end%only for small bodies
h=getx(b,'occ');%check saved data (& no ssd flag)
if isempty(h);h=getsb(b);%check ssd
if ~isnan(h{1});h=h{7};else;h=nan;warning('body %d not found for condition code',b);end
end%if
hx(jj)=h;
end%for jj
return

function hx=gettyp(bs)
hx=cell(size(bs));
for jj=1:numel(bs);b=bs(jj);
if b<1e6||imag(b);hx{jj}='';continue;end%only for small bodies
h=getx(b,'type');%check saved data (& no ssd flag)
if isempty(h);h=getsb(b);%check ssd
if ~isnan(h{1});h=h{8};else;h=nan;warning('body %d not found for spectral type',b);end
end%if
hx{jj}=strrep(h,'_','');%underscore denotes checked but nothing found
end%for jj
return


function cada=getclo(bs)
cada=cell(size(bs));
%ssd=params('ssd');params('ssd',true)%for getnum
for jj=1:numel(bs);b=bs(jj);if imag(b);cada{jj}='';continue;end
if b<1e6;cada{jj}=[];continue;end%only for small bodies
s=curl(sprintf('https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=%d;cad=1',b));
if jj/100==round(jj/100);fprintf('Close approach progress: %.0f%%\n',jj/numel(bs)*100);end
ca=regexp(s,'<tr>  (.+?)</tr>','tokens');cad=[];
for ii=1:numel(ca);
cc=regexp(ca{ii},'<.*?><.*?>(.+?)<','tokens');
cad(ii,1)=datenum(cc{1}{1}{1},'yyyy-mmm-dd HH:MM')-730486.5;
cad(ii,2)=getnum(cc{1}{3}{1});
cad(ii,3)=str2num(cc{1}{4}{1});
end
if isempty(cad);warning('no close approach data for body %d',b);end
cada{jj}=cad;
end%for jj
%params('ssd',ssd)
%if jj==1;cada=cada{:};end
return

function hx=getref(bs)
hx=cell(size(bs));
for jj=1:numel(bs);b=bs(jj);if imag(b);hx{jj}='';continue;end
h=getx(b,'ephref');%check saved data (& no ssd flag)
if isempty(h)&&params('ssd');%see what current file is
if b<1e6;tephf=getef;h=deblank(tephf.f(tephf.numbers==b,:));%get ephemeris from list
else;[h e]=getsb(b);if ~isnan(h{1});h=e;else;h='';end;end%get reference from ssd
end%if, no ephemeris found
if isempty(h);warning('file of last ephemeris for %d is unknown',b);end
hx{jj}=h;
end%for jj
return

function hx=getdat(bs)
hx=zeros(size(bs));
for jj=1:numel(bs);b=bs(jj);if imag(b);hx(jj)=nan;continue;end
h=getx(b,'ephdate');%check saved data (& no ssd flag), no ssd analogue
if isempty(h);h=nan;warning('date of last ephemeris for %d is unknown',b);end
hx(jj)=h;
end%for jj
return

function xx=qorient(bb,tt,ll)
%quick orientation calculation
%ll=1: Pole
%ll=2: [Node;Q;Pole] where Node=cross(Earth pole,Pole), Q=cross(Pole,Node)
%ll=3: Prime Meridian
%ll=4: [Node;Q;Pole] where Node=Prime Meridian
%ll=5: Equinox, cross(Pole,H)
%ll=6: [Node;Q;Pole] where Node=Equinox

sb=size(bb);if all(sb>2);bb=bb(1,:);elseif sb(1)>2;bb=bb';end

%match t and b 1-to-1 if numel(tt)==numel(bb)
tb=0;if iscell(tt);tt=[tt{1}:tt{3}:tt{2}];else;tt=tt(:)';tb=numel(tt)==numel(bb);end
[bb ~, it]=unique(bb);%do one body at time
for ibu=1:numel(bb)
%match times to body
i0=ibu==it;b=bb(ibu);if tb;t=tt(i0);else;t=tt;end

pp=getx(b,'pp');if isempty(pp);%check saved data (& no ssd flag)
adp=[getpck(b,'POLE_RA');getpck(b,'POLE_DEC');getpck(b,'PM')]*pi/180;%get data from pck file
if numel(adp)==9;a=adp(1);d=adp(2);p=adp(3);o=84381.448/3600*pi/180;
    %obliquity @J2000: Lieske, J., "Precession Matrix Based on IAU (1976) System of Astronomical Constants"
    %rotate from EME to EMO
    P=[cos(d)*cos(a);sin(d)*sin(o)+cos(d)*cos(o)*sin(a);cos(o)*sin(d)-cos(d)*sin(o)*sin(a)];%pole
    N=[P(2);-P(1);0];N=N/sqrt(sum(N.*N));Q=cross(P,N);%basis orthogonal to pole, N is pole node: pole x [0;0;1]
%    N=cross(qorient(399,0,1),P);N=N/sqrt(sum(N.*N));Q=cross(P,N);%basis orthogonal to pole, N is pole node: Earth pole x pole
    X=[-cos(p)*sin(a)-cos(a)*sin(d)*sin(p);cos(a)*cos(p)-sin(a)*sin(d)*sin(p);cos(d)*sin(p)];%prime meridian
    X=[X(1);cos(o)*X(2)+sin(o)*X(3);cos(o)*X(3)-sin(o)*X(2)];p=atan2(sum(X.*Q),sum(X.*N));%rotate to EMO % get angle @ J2000
    sav=params('save');params('save',0);X=ephem1(b,0,3);params('save',sav);%get H
    X=cross(P,cross(X(1:3,:),X(4:6,:)));X=X/sqrt(sum(X.*X));pe=atan2(sum(X.*Q),sum(X.*N));%equinox
    pp=[P;p;pe];putsbmb(b,'pp',pp);%pp =[pole vector;PM angle;eqx angle] @ J2000
else;warning('no analytic orientation data found for %d',b);pp=nan*ones(5,1);end
end
X=pp(1:3);cl2=ceil(ll/2);%X = pole, cl2 is Pole|PM|eqx output
if ll>1;P=X;p12=sqrt(P(1).^2+P(2).^2);X=[P(2,:)./p12;-P(1,:)./p12;0*p12];%need X = pole node
%if ll>1;P=X;X=cross(qorient(399,0,1),P);X=X/sqrt(sum(X.*X));%need X = pole node
if cl2>1;if cl2==2;w=pp(4)+2*pi/getrot(b)*t;elseif cl2==3;w=pp(5);t=0;end%need angle
X=X*cos(w)+cross(P,X)*sin(w);P=P*ones(size(t));%rotate along pole
end;end
if ll==2&b~=399;X=cross(qorient(399,0,1),P);X=X/sqrt(sum(X.*X));end%J2000 Node DCM
if cl2==ll/2;X=[X;cross(P,X);P];end%dcm
if cl2~=2;X=repmat(X,size(t));end%match size of t

if tb;xx(:,i0)=X;else;xx(:,:,i0)=repmat(X,[1 1 sum(i0)]);end%match 1-to-one, or all times per body
end%for ibu
return

function X=ephem(b,t,bvi)
%ephemeris
sb=size(b);
if all(sb>2);warning('Body input for ephemeris should be 1 x n or 2 x n');b=b(1,:);
elseif sb(1)>2;b=b';end

%match t and b 1-to-1 if numel(t)==size(b,2)
tb=0;if ~iscell(t);t=t(:).';tb=numel(t)==size(b,2);end

[b ii it]=unique(b','rows');b=b';%do for each unique body
if size(b,1)==1;b(2,:)=nan;end%b(2,:) is wrt body list
for ibu=1:size(b,2)
ii=ibu==it;%match times for each unique body
bib=b(1,ibu);b0b=b(2,ibu);%bib is target, b0b is wrt
if isnan(b0b)%wrt central body, no need to "tree" body center
%get state, X, from ephem1 function
if tb;X(:,ii)=ephem1(bib,t(ii),bvi);
else;X(:,:,ii)=repmat(ephem1(bib,t,bvi),[1 1 sum(ii)]);end%match 1-to-one, or all times per body

else%wrt some specified body, may need to tree
if tb;tt=t(ii);else;tt=t;end%match 1-to-one, or all times per body
if bib==b0b%target=wrt, X=0
    if iscell(t);tt=[t{1}:t{3}:t{2}];end
    Xb=zeros(6,numel(tt));Xb0=0;
else%compute ephemeris
cb=cbfun(bib);cb0=cbfun(b0b);%central body of target and wrt body
lcb=10;if (cb~=10|cb0~=10) & (cb==cb0|cb==b0b|bib==cb0);lcb=max(cb,cb0);end%lowest central body (planet or sun)

Xb =0;if lcb~=bib;Xb=ephem1(bib,tt,bvi);if lcb~=cb;Xb=Xb+ephem1(cb,tt,bvi);
        cbb=cbfun(cb);if lcb~=cbb;Xb=Xb+ephem1(cbb,tt,bvi);end
    end;end%get target state wrt lcb
Xb0=0;if lcb~=b0b;Xb0=ephem1(b0b,tt,bvi);if lcb~=cb0;Xb0=Xb0+ephem1(cb0,tt,bvi);
        cbb=cbfun(cb0);if lcb~=cbb;Xb0=Xb0+ephem1(cbb,tt,bvi);end
    end;end%get "wrt body" state wrt lcb 
end%if bib

%X-Xb0 is target state wrt "wrt body" state
if tb;X(:,ii)=Xb-Xb0;else;X(:,:,ii)=repmat(Xb-Xb0,[1 1 sum(ii)]);end%match 1-to-one, or all times per body
end;end

return

function [X t]=ephem1(bib,t,bvi)
n=6;if bvi==4;n=12;end
if imag(bib);if iscell(t);t=[t{1}:t{3}:t{2}];end;X=zeros(n,numel(t));return;end%non-body control point

%get saved data, ideph is list of body numbers, deph is spline data
id=params('ideph');d=params('deph');if isempty(d);clear d;end%clear double so can overwrite as struct

if params('ssd');[X t err]=mkeph(bib,t);%go straight to horizons
%flag if any data is out of bounds
if err(1);warning('Min Horizons ephem time for %d is %s\n',bib,datestr(err(1)+730486.5));
elseif err(2);warning('Max Horizons ephem time for %d is %s\n',bib,datestr(err(2)+730486.5));end
if bvi==4;warning('no accel or jerk from Horizons');end
return;end%return state

cb=cbfun(bib);%central body, treat planets and satellites differently
ii=find(id==bib);%see if saved data exists
if isempty(ii);
file=[params('bdir') 'ephem/' num2str(bib)];%see if file exists
if ~exist(file,'file')&bib>1e6;bib=getsb(bib);bib=bib{1};file=[params('bdir') 'ephem/' num2str(bib)];end%see if number changed
if ~exist(file,'file');
%te is saved time span, di is data to save, don't save nans
te={-1,60*365.25,.1+.4*(cb==10)+.5*(bib>1e6)};[di td]=mkeph(bib,te);di=[td;di];di(:,isnan(di(2,:)))=[];
if isempty(di);if iscell(t);t=[t{1}:t{3}:t{2}];end;X=nan*ones(n,numel(t));return;end%no good data
else;f=fopen(file);n=fread(f,2,'double')';di=fread(f,n(1:2),'double');fclose(f);end%file exists, read data, 1st two entries is data size

if cb~=10%gives better accuracy for Moons, but takes ~twice time
td=di(1,:);R=di(2:4,:);V=di(5:7,:);H=cross(R,V);H=H./([1;1;1]*sqrt(sum(H.*H)));
%mH=repmat(mean(H,2),size(td));[v e]=eig((H-mH)*(H-mH)');[~,mi]=min(diag(e));H=v(:,mi);%Rotate to Local Laplace plane
H=unit(mean(H,2));
h12=sqrt(H(1)^2+H(2)^2);N=[-H(2)./h12;H(1)./h12;0];Q=[-H(3).*N(2);H(3).*N(1);h12];dcm=[N Q H];
R=dcm.'*R;zr=R(3,:);R=R(1:2,:);r=sqrt(sum(R.*R));qr=unwrap(atan2(R(2,:),R(1,:)));
V=dcm.'*V;zv=V(3,:);V=V(1:2,:);v=sum(V.*R)./r;qv=(R(1,:).*V(2,:)-R(2,:).*V(1,:))./r.^2;
di(2:end,:)=[r;qr;zr;v;qv;zv];
else;dcm=[];end
ii=size(id,1)+1;id(ii,1)=bib;%add body to list

if params('keeps');d(ii,:)=RVspline6(di);d(ii,1).dcm=dcm;params('ideph',id);params('deph',d);%save spline
else;clear d;d(ii,:)=RVAspline(di);d(ii,1).dcm=dcm;params('ideph',[]);params('deph',[]);end
end%isempty ii

%convert {start,end,delta} to list,
if iscell(t);t=[t{1}:t{3}:t{2}];end;ts=[];ist=~1;
%see if data is out of bounds
tl=d(ii).breaks(1);if min(t)<tl;ist=ist|t<tl;warning('Min time for %d ephem spline is %s\n',bib,datestr(tl+730486.5));end
tl=d(ii).breaks(end);if max(t)>tl;ist=ist|t>tl;warning('Max time for %d ephem spline is %s\n',bib,datestr(tl+730486.5));end
if any(ist);ts=t(ist);t(ist)=[];end%ist is indeces of data outside of spline range

X=zeros(6,numel(t));
if bvi~=2||cb~=10;X(1:3,:)=ppval(d(ii,1),t);end
if bvi>1;X(4:6,:)=ppval(d(ii,2),t);end
if bvi==4;X=[X;ppval(d(ii,3),t);ppval(d(ii,4),t)];end

if cb~=10%gives better accuracy for Moons, but takes ~twice time
dcm=d(ii,1).dcm;
c=cos(X(2,:));s=sin(X(2,:));r=X(1,:);if bvi~=2;X(1:3,:)=dcm*[r.*c;r.*s;X(3,:)];end
if bvi>1;dr=X(4,:);dq=X(5,:);rdq=r.*dq;X(4:6,:)=dcm*[dr.*c-rdq.*s;dr.*s+rdq.*c;X(6,:)];end
if bvi==4;d2r=X(7,:);d2q=X(8,:);r_=d2r-rdq.*dq;q_=r.*d2q+2*dr.*dq;X(7:9,:)=dcm*[r_.*c-q_.*s;r_.*s+q_.*c;X(9,:)];
    r_=X(10,:)-3*dr.*dq.^2-3*rdq.*d2q;q_=r.*X(11,:)+3*d2r.*dq+3*dr.*d2q-rdq.*dq.^2;
    X(10:12,:)=dcm*[r_.*c-q_.*s;r_.*s+q_.*c;X(12,:)];end
end

%get data outside of spline range from Horizons
if any(ist);sav=params('save');params('save',0);if bvi==4;warning('no accel or jerk from Horizons');end
    Xs=mkeph(bib,ts);X(:,~ist)=X;X(1:6,ist)=Xs;params('save',sav);
end
return

function pp=RVspline6(X)
%5th order spline fiting R,V at segment endpoints, continuous A,J
%nonuniform step ok
t=X(1,:);R=X(2:4,:);V=X(5:7,:)*86400;
dt=diff(t);n=numel(t)-1;o=[1;1;1];
c0=R(:,1:end-1);c1=V(:,1:end-1);%first two coeff
%4 coeff remain
R_=(c0+c1.*dt(o,:)-R(:,2:end))./dt(o,:).^2;V_=(c1-V(:,2:end))./dt(o,:)/2;
%R_+c2+c3.*dt+c4.*dt2,V_+c2+3/2*c3.*dt+2*c4.*dt2,%satisfy R and V at end of segment
%c5(i)*dt^3+6*R_-6*V_+c2(i)=c2(i+1) , (3*c5*dt^3+8*R_-6*V_+2*c2)/dt=(c5*dt^3-4*R_+2*V_-2*c2)/dt%continuous 2nd & 3rd derivative
%c5(i)=(6*(V_-R_)-c2(i)+c2(i+1))/dt^3 , (12*V_-10*R_-c2(i-1)+3*c2(i))/dt=(8*V_-10*R_-3*c2(i)+c2(i+1))/dt%combine into single eqn in c2
%(-15*R_+16*V_-2*c2(i-1)+3*c2(i))/dt^2=(15*R_-14*V_+3*c2(i)-2*c2(i+1))/dt^2 %continuous 4th derivative at end points
ii=2:n;b=(12*V_(:,ii-1)-10*R_(:,ii-1))./dt(o,ii-1)-(8*V_(:,ii)-10*R_(:,ii))./dt(o,ii);
s=[1./dt(ii-1) -3./dt(ii-1)-3./dt(ii) 1./dt(ii)];%s*[c2(i-1);c2(i);c2(i+1)]=b, tridiagonal
jj=[ii-1 ii ii+1];ii=[ii ii ii];%sparsity pattern x,y,z
%continuous 4th der at 2
b_=(-15*R_(:,1)+16*V_(:,1))./dt(o,1).^2-(15*R_(:,2)-14*V_(:,2))./dt(o,2).^2;
s_=[2./dt(1).^2 -3./dt(1).^2+3./dt(2).^2 -2./dt(2).^2];
b=[b_ b];ii=[1 1 1 ii];jj=[1 2 3 jj];s=[s_ s];
%continuous 4th der at n
b_=(-15*R_(:,n-1)+16*V_(:,n-1))./dt(o,n-1).^2-(15*R_(:,n)-14*V_(:,n))./dt(o,n).^2;
s_=[2./dt(n-1).^2 -3./dt(n-1).^2+3./dt(n).^2 -2./dt(n).^2];
b=[b b_];ii=[ii n+[1 1 1]];jj=[jj n+[-1 0 1]];s=[s s_];
s=sparse(jj,ii,s);c2=b/s;%solve for c2
dt=dt(o,:);dt2=dt.*dt;dt3=dt2.*dt;
c22=c2(:,2:n+1);c2=c2(:,1:n);c5=(6*(V_-R_)-c2+c22)./dt3;%solve for c5
c5dt3=c5.*dt3;c3=(2*V_-4*R_-2*c2+c5dt3)./dt;c4=(3*R_-2*V_+c2-2*c5dt3)./dt2;%solve for c3 and c4
pp.form='pp';pp.breaks=t;pp.coefs=[c5(:) c4(:) c3(:) c2(:) c1(:) c0(:)];pp.pieces=numel(t)-1;pp.order=6;pp.dim=3;pp.dcm=[];pp1=pp;
pp(2)=pp1;pp(2).coefs=[5*c5(:) 4*c4(:) 3*c3(:) 2*c2(:) c1(:)]/86400;pp(2).order=5;%derivatives
pp(3)=pp1;pp(3).coefs=[20*c5(:) 12*c4(:) 6*c3(:) 2*c2(:)]/86400^2;pp(3).order=4;
pp(4)=pp1;pp(4).coefs=[60*c5(:) 24*c4(:) 6*c3(:)]/86400^3;pp(4).order=3;
return

function pp=RVspline5(X)
%4th order spline fiting R,V at segment endpoints, continuous A
%nonuniform step ok
%only one boundary cond'n, not as accurate
t=X(1,:);R=X(2:4,:);V=X(5:7,:)*86400;
dt=repmat(diff(t),3,1);dt2=dt.*dt;
c0=R(:,1:end-1);c1=V(:,1:end-1);
R_=(c0+c1.*dt-R(:,2:end))./dt2;V_=(c1-V(:,2:end))./dt/2;
%R_+c2+c3.*dt+c4.*dt2,V_+c2+3/2*c3.*dt+2*c4.*dt2,6*(R_-V_)+c2(i)=c2(i+1),
c2=((V_(:,2)-2*R_(:,2))./dt(:,2)-(10*R_(:,1)-9*V_(:,1))./dt(:,1)).*dt(:,1)/2;
RV6=6*(R_-V_);c2=cumsum([c2 RV6],2);c2=c2(:,1:end-1);
%c2=((V_(:,end)-2*R_(:,end))./dt(:,end)-(3*V_(:,end-1)-2*R_(:,end-1))./dt(:,end-1)).*dt(:,end-1)/2;
%c2=cumsum([c2 -RV6(:,end-1:-1:1)],2);c2=c2(:,end:-1:1);%/2+c2_/2;
c3=(V_-2*R_-c2)*2./dt;c4=(3*R_-2*V_+c2)./dt2;
pp.form='pp';pp.breaks=t;pp.coefs=[c4(:) c3(:) c2(:) c1(:) c0(:)];pp.pieces=numel(t)-1;pp.order=5;pp.dim=3;pp.dcm=[];
pp(2)=pp;pp(2).coefs=[4*c4(:) 3*c3(:) 2*c2(:) c1(:)]/86400;pp(2).order=4;
return

function pp=RVAspline(X)
%5th order spline using R,V,A at beginning and end of segments
%approximates A with central diff of R and V
%assumes uniform time step
t=X(1,:);R=X(2:4,:);V=X(5:7,:)*86400;
dt=diff(t);if any(abs(diff(dt))>1e-9);error('nonuniform timestep for spline');end
dt=dt(1);V=V*dt;n=numel(t);%uniform timestep, scale V

c2=R(:,1:n-2)-2*R(:,2:n-1)+R(:,3:n)+(V(:,1:n-2)-V(:,3:n))/4;%central diff for A (c2 = coeff on t^2 term)
c21=-23/4*R(:,1)+4*R(:,2)+7/4*R(:,3)-3*V(:,1)-4*V(:,2)-V(:,3)/2;%initial A (uses first 3 R,V)
c2n=7/4*R(:,n-2)+4*R(:,n-1)-23/4*R(:,n)+V(:,n-2)/2+4*V(:,n-1)+3*V(:,n);%final A (uses last 3 R,V)
c0=R(:,1:n-1);R2=c0-R(:,2:n);c1=V(:,1:n-1);V2=V(:,2:n);A2=[c2 c2n];c2=[c21 c2];%coeff from initial R,V,A and temp final R,V,A on segs

c3=-10*R2-6*c1-4*V2-3*c2+A2;c4= 15*R2+8*c1+7*V2+3*c2-2*A2;c5=-6*R2-3*c1-3*V2-c2+A2;%solve for R2,V2,A2 at seg endpts
c=zeros(3*n-3,6);c(:,6)=c0(:);c(:,5)=c1(:)/dt;c(:,4)=c2(:)/dt^2;c(:,3)=c3(:)/dt^3;c(:,2)=c4(:)/dt^4;c(:,1)=c5(:)/dt^5;%unscale time

pp.form='pp';pp.breaks=t;pp.coefs=c;pp.pieces=numel(t)-1;pp.order=6;pp.dim=3;pp.dcm=[];
d=diag(5:-1:1,1);d=d(:,2:end);pp(2)=pp;pp(2).coefs=c*d/86400;pp(2).order=5;
return

function [xx t]=orient(bb,tt,ll)
%ll=1: Pole
%ll=2: [Node;Q;Pole] where Node=cross([0;0;1],Pole), Q=cross(Pole,Node)
%ll=3: Prime Meridian
%ll=4: [Node;Q;Pole] where Node=Prime Meridian
%ll=5: Equinox, cross(Pole,H)
%ll=6: [Node;Q;Pole] where Node=Equinox

sb=size(bb);if all(sb>2);bb=bb(1,:);elseif sb(1)>2;bb=bb';end

%get saved data, idori is list of body numbers, dori is spline data
id=params('idori');d=params('dori');if isempty(d);clear d;end%clear double so can overwrite as struct

tb=0;if ~iscell(tt);tt=tt(:)';tb=numel(tt)==numel(bb);end%match t and b 1-to-1 if numel(t)==size(b,2)
[bb ~, it]=unique(bb);%do for each unique body
for ibu=1:numel(bb)
i0=ibu==it;bib=bb(ibu);if tb;t=tt(i0);else;t=tt;end%match times for each unique body

%small body orientation not in Horizons, may be in pck file. write xx 1-to-1 if tb
if bib>1e3;X=qorient(bib,t,ll);if tb;xx(:,i0)=X;else;xx(:,:,i0)=repmat(X,[1 1 sum(i0)]);end;continue;end

%cl2 is Pole|PM|eqx output, ii is cl2 data exists, ij is any data exists
cl2=ceil(ll/2);[ij jj]=find(id==bib);ii=ij(jj==cl2);ssd=params('ssd');
if isempty(ii)||ssd%read from file or horizons
if ssd;te=t;else;te={-.5,60*365.25,.1};end%te is saved time span

fls={'pol' 'pm' 'eqx'};err=0;%file names
file=sprintf('%sephem/%s%d',params('bdir'),fls{cl2},bib);rf=~ssd&&exist(file,'file');
if rf;f=fopen(file);n=fread(f,2,'double')';di=fread(f,n(1:2),'double');fclose(f);%read from file if rf, 1st two entries is data size
else;switch cl2;%read from horizons
    case 1;[X td err]=mkeph(bib,te,[0 90]);%pole
    sp=X(3,:)<0;if any(sp);X(:,sp)=-X(:,sp);if any(~sp);warning('pole crosses ecliptic');end;end%z component is positive for sqrt
    case 2;[X td err]=mkeph(bib,te,[0 0]);P=orient(bib,te,1);%prime meridian
    case 3;sav=params('save');params('save',0);[X td]=ephem1(bib,te,3);params('save',sav);%get H
    P=orient(bib,te,1);X=cross(P,cross(X(1:3,:),X(4:6,:)));%equinox
    end;X=X(1:3,:);X=X./([1;1;1]*sqrt(sum(X.*X)));%unit vector
    if cl2==1;di=[td;X(1:2,:)];if ll==2;P=X;X=[P(2,:);-P(1,:);0*td];X=X./([1;1;1]*sqrt(sum(X.*X)));end%X is pole node, Pole x [0;0;1]
    %if cl2==1;di=[td;X(1:2,:)];if ll==2;P=X;X=cross(qorient(399,0*td,1),P);X=X./([1;1;1]*sqrt(sum(X.*X)));end%X is pole node, Pole x [0;0;1]
    %set up basis orthogonal to pole, N is pole node (Pole x [0;0;1]) Q is P x N
    else;p12=sqrt(P(1,:).^2+P(2,:).^2);N=[P(2,:)./p12;-P(1,:)./p12;0*p12];Q=[-P(3,:).*N(2,:);P(3,:).*N(1,:);-p12];%Pole x [0;0;1]
    %else;N=cross(qorient(399,0*td,1),P);N=N./repmat(sqrt(sum(N.*N)),3,1);Q=cross(P,N);%Earth pole x pole
    w=unwrap(atan2(sum(X.*Q),sum(X.*N)));di=[td;w];%clock angle (period has to be longer than 2X sample for unwrap)
    end%cl2
end%rf

di(:,isnan(di(2,:)))=[];%data to save, don't save nans
if isempty(di);if iscell(t);t=[t{1}:t{3}:t{2}];end;X=nan*ones(3+6*(cl2==ll/2),numel(t));%no data, output nans
if tb;xx(:,i0)=X;else;xx(:,:,i0)=repmat(X,[1 1 sum(i0)]);end;continue;end

%see if times are out of horizons range if ssd
if ssd&err(1);warning('Min Horizons orient time for %d is %s\n',bib,datestr(err(1)+730486.5));
elseif ssd&err(2);warning('Max Horizons orient time for %d is %s\n',bib,datestr(err(2)+730486.5));
%save data if no error and didn't read from existing file
elseif ~rf&params('save')&~err;f=fopen(file,'w');fwrite(f,size(di),'double');fwrite(f,di,'double');fclose(f);end

if ~ssd%make spline if not ssd flag
    if ~isempty(ij);ii=ij(1);else;ii=size(id,1)+1;end;id(ii,cl2)=bib;%see if body already exists in idori
    d(ii,cl2)=spline(di(1,:),di(2:end,:));
    if params('keeps');params('idori',id);params('dori',d);end%save data
end
end%ii|ssd

%convert {start,end,delta} to list
if iscell(t);t=[t{1}:t{3}:t{2}];end;ts=[];ist=~1;
if ~ssd
%see if data is out of bounds
tl=d(ii,cl2).breaks(1);if min(t)<tl;ist=ist|t<tl;warning('Min time for %d orient spline is %s\n',bib,datestr(tl+730486.5));end
tl=d(ii,cl2).breaks(end);if max(t)>tl;ist=ist|t>tl;warning('Max time for %d orient spline is %s\n',bib,datestr(tl+730486.5));end
if any(ist);ts=t(ist);t(ist)=[];end%ist is indeces of data outside of spline range

if ~id(ii,1);P=orient(bib,t,1);else;P=ppval(d(ii,1),t);end%get new pole data if not already in spline
p12=P(1,:).^2+P(2,:).^2;P=[P(1,:);P(2,:);sqrt(1-p12)];X=P;%pole vector
if ll>1;p12=sqrt(p12);X=[P(2,:)./p12;-P(1,:)./p12;0*p12];%X is pole node: Pole x [0;0;1]
if cl2>1;w=ppval(d(ii,cl2),t);X=[1;1;1]*cos(w).*X+[1;1;1]*sin(w).*[-P(3,:).*X(2,:);P(3,:).*X(1,:);-p12];%get angle and rotate
%if ll>1;X=cross(qorient(399,0*t,1),P);X=X./repmat(sqrt(sum(X.*X)),3,1);%X is pole node: Earth pole x Pole
%if cl2>1;w=ppval(d(ii,cl2),t);X=[1;1;1]*cos(w).*X+[1;1;1]*sin(w).*cross(P,X);%get angle and rotate
end;end
end%ssd
if ll==2&bb~=399;X=cross(qorient(399,0*t,1),P);X=X./([1;1;1]*sqrt(sum(X.*X)));end
if cl2==ll/2;X=[X;cross(P,X);P];end%dcm

if any(ist);sav=params('save');params('save',0);params('ssd',1);%reset flags
Xs=orient(bib,ts,ll);X(:,~ist)=X;X(:,ist)=Xs;params('save',sav);params('ssd',ssd);end%get horizons data for out of range times

if tb;xx(:,i0)=X;else;xx(:,:,i0)=repmat(X,[1 1 sum(i0)]);end%match 1-to-one, or all times per body
end%for ibu

return

function [X tt err]=mkeph(n,t,ll)
if n<1e3||n>1e6%read from horizons
err=0;tephf=getef;%get ephemeris time span
cb=cbfun(n);%get data wrt central body
if n>1e6&n<2e6;CAP='%3BCAP';else;CAP='';end%comets are dumb
if n>1e6;tlim=tephf.sb;DES='DES=';ieph=1;%time limits for small bodies
else;ieph=tephf.numbers==n;tlim=tephf.tl(ieph,:);DES='';end%time limits for bog bodies

%tt is expected time output from horizons, can only run 400 times in list, or 90,000 times in {start,end,delta}
if isnumeric(t);[tt ist]=sort(t);mnt=4e2;else;tt=t{1}:t{3}:t{2};mnt=9e4;
if numel(tt)<2;warning('invalid {start time, end time, delta time} input to Horizons');end;end
if all(~ieph);warning('no horizons ephemeris file for %d',n);X=nan*ones(6,numel(tt));err=1;return;end
nt1=sum(tt<tlim(1));nt2=sum(tt>tlim(2));err=~~[nt1 nt2].*tlim;%nt1 and nt2 are out of time range
n2t=numel(tt)-nt2;t12=round(linspace(nt1,n2t,ceil((n2t-nt1)/mnt)+1));%t12 breaks time span into different runs
if n>1e6;[~,tnow]=min(abs(now-730486.5-tt));t12=unique([tnow-1 t12]);end

X=zeros(7,0);n12=numel(t12);
fprintf('Retrieving data for Body %d from Horizons\n',n);
for ti=1:n12-1;%only do mnt at a time
%time input for url
if iscell(t);tstr=sprintf('&START_TIME=''JD%%20%.9f''&STOP_TIME=''JD%%20%.9f''&STEP_SIZE=''%d''',...
tt(t12(ti)+1)+2451545,tt(t12(ti+1))+2451545,t12(ti+1)-t12(ti)-1);
else;tstr=['&TLIST=''' sprintf('%.9f%%0A',tt(t12(ti)+1:t12(ti+1))+51544.5) ''''];end
%set target as center of body and observer as ll (lat,lon) for orientation
if nargin>2;cb=n;coord='coord';llstr=sprintf('&COORD_TYPE=''GEODETIC''&SITE_COORD=''%.8f,%.8f,0''',ll);else;coord='';llstr='';end

url=sprintf(['https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1'...horizons url mumbojumbo
'&COMMAND=''%s%d%s'''...
'&CENTER=''%s@%d'''...
'%s',...
'&MAKE_EPHEM=''YES'''...
'&TABLE_TYPE=''VECTORS'''...
'%s',...
'&OUT_UNITS=''KM-S'''...
'&VECT_TABLE=''2'''...
'&REF_PLANE=''ECLIPTIC'''...
'&REF_SYSTEM=''J2000'''...
'&VECT_CORR=''NONE'''...
'&VEC_LABELS=''NO'''...
'&CSV_FORMAT=''NO'''...
'&OBJ_DATA=''NOPE'''],DES,n,CAP,coord,cb,llstr,tstr);
ss=curl(url);i1=strfind(ss,'$$SOE')+6;i2=strfind(ss,'$$EOE')-1;%run horizons url command and extract goodies
if ~isempty(i1);X=[X sscanf(ss(i1:i2),'%f%*40c%e%e%e%e%e%e\n',[7 inf])];%successful run
%tri=regexp(ss,'Center radii\s*: ([\d.]+) x ([\d.]+) x ([\d.]+) ','tokens');putsbmb(cb,{'tri'},{str2num(cat(1,tri{1}{:}))});%triaxial also in horizons
else;err_str=ss,warning('Horizons error');X=nan*ones(6,numel(tt));return;end
if ti~=n12-1;fprintf('Horizons progress: %.0f%%\n',(t12(ti+1)-nt1)/(n2t-nt1)*100);end
end%ti
tt=[tt(1:nt1) X(1,:)-2451545 tt(end-nt2+1:end)];
X=[nan*ones(6,nt1) (1-2*(nargin>2))*X(2:7,:) nan*ones(6,nt2)];%pad out of range data with nans, if orient X=-X
if ~iscell(t);X(:,ist)=X;tt(ist)=tt;end%output times in same order as input

else%Lagrange point
  err=0;if isnumeric(t);tt=t;ist=1:numel(t);else;tt=t{1}:t{3}:t{2};end
  X=Lagrange(n,tt);
end

if ~exist([params('bdir') 'ephem'],'dir');mkdir([params('bdir') 'ephem']);end%make "ephem" directory
if params('save')&&nargin<3%save ephemeris data
f=fopen(sprintf('%sephem/%d',params('bdir'),n),'w');
d=[tt;X];if ~iscell(t);d=d(:,ist);end;d(:,isnan(X(1,:)))=[];%save times in order, don't save nans
fwrite(f,size(d),'double');fwrite(f,d,'double');fclose(f);%write data, 1st two entries is data size
if n>1e6;[n e]=getsb(n);putsbmb(n{1},{'ephref' 'ephdate'},{e now-730486.5});%always update small bodies
elseif n<1e3;putsbmb(n,{'ephref' 'ephdate'},{deblank(tephf.f(ieph,:)) now-730486.5});%update reference file and time
else;putsbmb(n,{'ephref' 'ephdate'},[getref(floor(n/10)) now-730486.5]);end
end%save
return

function X=Lagrange(b,t)
i=mod(b,10);%Lagrange pt number
b=floor(b/10);%central body
m1=getgm(cbfun(b));m2=getgm(b);%masses

if i<4%L1,L2,L3
%calculate distance for circular (a), then scale by distance (r)
m=m1+m2;m1=m1/m;m2=m2/m;%non-dim
m1i=m1;m2i=m2;if i==1||i==3;m2i=-m2;end;if i==3;m1i=-m1i;end%account for direction of pull
if i<3;L=exp(log(m2/m1/3)/3);if i==1;L=-L;end;elseif i==3;L=-2;end%initial guess
L=L+1e-42i;%complex step
for cc=1:99
e=m1+L-m1i./(L+1).^2-m2i./L.^2;%balance dynamics, rotation - gravity
de=imag(e)/1e-42;e=real(e);if all(abs(e)<9*eps);break;end
%dL=e./de;ii=abs(dL)>xL;if any(ii);dL(ii)=sign(dL(ii))*xL;end
L=L-e./de;
end
L=real(L);if any(abs(e)>9*eps);warning('Lagrange not converged');end
%e=m1+L-m1i./(L+1).^2-m2i./L.^2;dL=-imag(e)/1e-42./de;L=L+dL*1e-42i;X=L*ephem1(b,t+1e-42i,1);X=[real(X);imag(X)/1e-42/86400];
X=L*ephem1(b,t,3);%point is along position of secondary wrt primary

else%L4,L5
%triangular points
X=ephem1(b,t+1e-42i,3);R=X(1:3,:);r=sqrt(sum(R.^2));V=X(4:6,:);H=cross(R,V);
Y=unit(cross(H,R));y=r.*sqrt(.75);%distance from axis connecting bodies
if i==5;y=-y;end%L5 trails
X=-R/2+Y.*[y;y;y];X=[real(X);imag(X)/1e-42/86400];%complex step
end%i

return

function ed=getephdat(b)
%ephemeris header data
ef=getef;ii=ef.numbers==b;ef=lower(regexp(ef.f(ii,:),'^[^-._ ]*','match'));%only need first file if merged
if ~isempty(ef);ef=ef{1};else;warning('no ephemeris header file for %d',b);ed=[];return;end
ed=params(ef);%see if file data already saved
if isempty(ed)&b<400;%inner planets and barycenters
    %read in DE### header data
    ef=char(regexp(ef,'de\d+','match'));%deXXX format
    efdir=['ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/' ef '/'];%directory
    [d err]=curl([efdir 'header.' ef(3:end)]);%standard filename is "header.XXX"
    if err;fils=regexp(curl(efdir),['header.' ef(3:end) '_?\d*'],'match');%sometimes in "header.XXX_YYY" format
    if isempty(fils);warning('ephemeris header file "%s" for %d not found',upper(ef),b);delete(hout);return
    else;[d err]=curl([efdir fils{1}]);end
    end

    ii=strfind(d,'GROUP   1040')+12;[nv,~,~,c]=sscanf(d(ii+[0:9]),'%d',1);   
    vars=reshape(sscanf(d(ii+c:end),'  %6c',nv),[6 nv]);%read in nv variable names (each 6 chars)
    ii=strfind(d,'GROUP   1041')+12;[nv,~,~,c]=sscanf(d(ii+[0:9]),'%d',1); 
    vals=sscanf(d(ii+c:end),'%eD%e',[2 nv]);vals=vals(1,:).*10.^vals(2,:);%read in values, assumes "D" exponent
    ed.vars=vars';ed.vals=vals';
    %write gm values
    au=ed.vals(uniqstr(ed.vars,'AU'));params('AU',au);
    ed.gm10=ed.vals(uniqstr(ed.vars,'GMS'))*au^3/86400^2;%sun
    ed.gm3=ed.vals(uniqstr(ed.vars,'GMB'))*au^3/86400^2;%E-M bary
    emr=ed.vals(uniqstr(ed.vars,'EMRAT'));%E/M ratio
    ed.gm301=ed.gm3/(1+emr);ed.gm399=ed.gm3*emr/(1+emr);%Moon and Earth
    %barycenters
    ed.gm0=ed.gm10+ed.gm3;
    for ii=[1:2 4:9];GM=['gm' num2str(ii)];gm=ed.vals(uniqstr(ed.vars,GM))*au^3/86400^2;eval(['ed.' GM '=gm;']);ed.gm0=ed.gm0+gm;end
    ed.gm199=ed.gm1;ed.gm299=ed.gm2;%Mercury and Venus are Trouble
    %radii and J2
    ed.rad10=ed.vals(uniqstr(ed.vars,'ASUN'));ed.j2_10=ed.vals(uniqstr(ed.vars,'J2SUN'));
    ed.rad399=ed.vals(uniqstr(ed.vars,'RE'));ed.j2_399=ed.vals(uniqstr(ed.vars,'J2E'));
    ed.rad301=ed.vals(uniqstr(ed.vars,'AM'));ed.j2_301=ed.vals(uniqstr(ed.vars,'J2M'));
    ed.rad199=ed.vals(uniqstr(ed.vars,'RAD1'));ed.rad299=ed.vals(uniqstr(ed.vars,'RAD2'));
    params(ef,ed);%save data
elseif isempty(ed);%outer planets and satellites
    %check default directory and eph filename
    [ss,err]=curl(['ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/' ef '.txt']);
    %check default directory for similar filename
    if err;fils=char(regexp(curl('ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/'),'\S+(?=.txt)','match'));
        ii=uniqstr(fils,ef);if isempty(ii);fils=fils(end:-1:1,:);ii=uniqstr(fils,ef(1:3));end%look for first 3 letters of highest # file
        if isempty(ii);err=1;else;[ss err]=curl(['ftp://ssd.jpl.nasa.gov/pub/eph/satellites/nio/LINUX_PC/' deblank(fils(ii,:)) '.txt']);end
    end;%could also check ftp://ssd.jpl.nasa.gov/pub/eph/satellites/rckin/rckin.*.log, but GM likely 0
    %ss is ephemeris header file 
    if ~err;
    %get GM from "Bodies on the File" table: skip space,keep 3 #s,skip some space, keep some #s, decimal, & non-space,skip space
    gm=regexp(ss,'\s(\d{3})\s+((\d+\.\S*))\s','tokens');gm=cat(1,gm{:});
    ed=cell2struct(cellfun(@str2num,gm(:,2),'un',0),strcat('gm',gm(:,1)));%make struct of ed.gm(body#)=GM
    b=floor(b/100)*100+99;%if mod(b,100)==99%get J2 and radius of planet
    ed.(['j2_' num2str(b)])=sscanf(ss(strfind(ss,['J' num2str(floor(b/100)) '02']):end),'%*s%e',1);
    ed.(['rad' num2str(b)])=sscanf(ss(strfind(ss,'RADIUS'):end),'%*s%e',1);
    %end%if
    else;warning('ephemeris header file "%s" for %d not found',upper(ef),b);end%err
    params(ef,ed);
end
return

function tephf=getef
%Ephemeris file and time spans used by Horizons
tephf=params('tephf');if isempty(tephf)
ss=curl('https://ssd.jpl.nasa.gov/eph_spans.cgi?id=A');%Planets
%Read in bodynumber, begin time, " not " or " to " flag, end time, and file
ss=regexp(ss,'<td.*?>(\d+)</td>.*?<\w\w?>(.*?) (not|to) (.*?)</\w\w?>&nbsp;.*?&nbsp;(.*?)&nbsp;','tokens');sss=cat(1,ss{:});
%Mercury and Venus are Trouble (no ephemeris file, but still in table)
for ii=find(strcmp(sss(:,3),'not'))';jj=find(strcmp(sss(:,1),[sss{ii,1} '99']));sss(ii,:)=sss(jj,:);sss{ii,1}=sss{jj,1}(1);end;sss(:,3)=[];
ss=curl('https://ssd.jpl.nasa.gov/eph_spans.cgi?id=B');%Satellites
%Read in bodynumber, begin time, end time, and file (no need to flag " not " or " to "
ss=regexp(ss,'<td.*?>(\d+)</td>.*?<tt>(.*?) to (.*?)</tt>&nbsp;.*?&nbsp;(.*?)&nbsp;','tokens');ss=cat(1,sss,ss{:});
ss=regexprep(ss,{'B.C. \d{4}' 'A.D. '},{'0000' ''});%matlab doesn't do B.C., convert to days from J2000
%ss=strrep(ss,'--','-Aug-');%fix error on website
tephf.numbers=cellfun(@str2num,ss(:,1));tephf.f=char(ss(:,4));tephf.tl=[datenum(ss(:,2),'yyyy-mmm-dd') datenum(ss(:,3),'yyyy-mmm-dd')]-730486.5;
ss=curl('https://ssd.jpl.nasa.gov/eph_spans.cgi?id=D');%Small body time spans (eph file read from getsb)
ss=regexp(ss,'&nbsp;(\S+) to (\S+)&nbsp;','tokens');tephf.sb(1:2)=datenum(ss{:},'yyyy-mmm-dd')-730486.5;
params('tephf',tephf);end
return

function sd=getsatdat
%Satellite data page, has GM and mean radius (& density, magnitude, albedo)
sd=params('satdat');
if isempty(sd);sd=curl('https://ssd.jpl.nasa.gov/?sat_phys_par');params('satdat',sd);end
return

function d=getpck(b,v)
%NAIF pck file, analytic orientation & contains some small body radii and orientation data not on Horizons
pck=params('pck');if isempty(pck)
pck=regexp(curl('ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/'),'pck\d+.tpc','match');%get latest file in directory
pck=curl(['ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/' pck{end}]);%read contents to string
pck=strrep(strrep(strrep(pck,sprintf('\n'),';'),'(','['),')',']');%turn newlines into ; and parens into brackets for matlabese
pck=strrep(strrep(pck,'2431010','2000243'),'9511010','2000951');%change # of Ida and Gaspra to SPK-ID
params('pck',pck);end%save
v=sprintf('BODY%d_%s',b,upper(v));%variable name
d=regexp(pck,['(?<=\\begindata.*?)' v '.*?];'],'match');%just single param
%d=regexp(pck,['\\begindata[;\s]+(' v{1} '.*?;)[;\s]+\\begintext'],'tokens');%everything in data block
if ~isempty(d);eval([d{end} 'd=' v ';']);else;d=[];end
return

function o=getx(b,s)
%get data from mb or sb structure
x=getsbmb(b);ii=find(b==x.numbers);%get datastructure and see if body is there
if ~isempty(ii)&&~params('ssd')&&isfield(x,s)%see if should skip saved data and if datafield exists
    %see if entry exists in field for body
    if isnumeric(x.(s))&&size(x.(s),2)>=ii&&~isnan(x.(s)(1,ii));o=x.(s)(:,ii);
    elseif ischar(x.(s))&&size(x.(s),1)>=ii&&~isspace(x.(s)(ii,1));o=deblank(x.(s)(ii,:));
else;o=[];end
else;o=[];end
return

function x=getmb
%major body list
nn=curl('https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=MB');
nn=regexp(nn,'\n\s*(\d){1,3}\s+(.*?)\s\s','tokens');nn=cat(1,nn{:});%keep 1--3 numbers, skip spaces, keep stuff until two spaces
nums=cellfun(@str2num,nn(:,1))';nams=nn(:,2);%numbers and names
ii=nums>9;nums=[nums(ii) nums(~ii)];nams=[nams(ii);nams(~ii)];%match barycenters last
x=getsbmb('mb');b=x.numbers;[ii jj]=ismember(nums,b);%see if body numbers already exist in data structure
nums=[b nums(~ii)];n(jj(ii),1)=nams(ii);nams=[n;nams(~ii)];%put new bodies last
x.numbers=nums;x.names=char(nams);params('mb',x);%params('mbmod',true);%save data
return

function [nnf eph s]=getsb(n)
%small body web pages
if ischar(n);srch=upper(urlencode(n));else;srch=num2str(n);end%name or number
s=curl(['https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=' srch]);

nam=regexp(s,'+1"><b>[^<]+','match');%Name is bigger font "+1"
if ~isempty(nam);nam=nam{1}(8:end);
num=regexp(s,'>\d{7}<','match');num=str2double(num{1}(2:end-1));
if ~ischar(n)&&n~=num;fprintf('SPK ID of %d has been changed to %d\n',n,num);n=[n num];else;n=num;end
nnf={num nam};%{number name gm rad rot h occ type}
%read in data off table
a=regexp(s,'>GM<.*?>(\d.*?)<','tokens');if isempty(a);a=0;else;a=str2num(a{1}{1});end;nnf=[nnf a];
a=regexp(s,'>diameter<.*?>([\d.]+)<','tokens');if isempty(a);a=0;else;a=str2num(a{1}{1})/2;end;nnf=[nnf a];
a=regexp(s,'>rot_per<.*?>([\d.]+)<','tokens');if isempty(a);a=0;else;a=str2num(a{1}{1})/24;end;nnf=[nnf a];
a=regexp(s,'>H<.*?>([\d.]+)<','tokens');if isempty(a);a=nan;else;a=str2num(a{1}{1});end;nnf=[nnf a];
a=regexp(s,'>condition code<.*?sp;(\d)&nb','tokens');if isempty(a);a=nan;else;a=str2num(a{1}{1});end;nnf=[nnf a];
a=regexp(s,'>spec_T<.*?>([a-z_A-Z:\)\( ]+)<','tokens');if ~isempty(a);a=a{1}{1};end
b=regexp(s,'>spec_B<.*?>([a-z_A-Z:\)\( ]+)<','tokens');if isempty(b);b='_';elseif isempty(a);b=b{1}{1};else;b=['/' b{1}{1}];end;nnf=[nnf [a b]];
eph=regexp(s,'>Reference: <.*?>.*?>.*?>(.*?)<','tokens');eph=eph{1}{1};%Ephemeris ID
putsbmb(n,{'numbers' 'names' 'gm' 'rad' 'rot' 'h' 'occ' 'type'},nnf);%save data
else;nam=regexp(s,'">([^<]+)</a></td>','tokens');%see if returned multiple matches
nnf={nan char([nam{:}])};eph='';return;end
return

function nn=uniqsb(b)
b=regexprep(b,'[()]','');
bs={b,[b '*'],['*' b '*']};%search exact, beginning, fragment
for b=bs
nn=getsb(b{1});%check ssd
if isnan(nn{1})&~isempty(nn{2});matches=nn{2},warning([b{1} ' returned multiple matches'])
%n is {number name provisional} of firt match, keep first of those
n=regexp(nn{2}(1,:),'(\d*\s?)([^\(]*)(.*)','tokens');n=n{1}{find(~cellfun('isempty',n{1}),1)};
nn=getsb(deblank(n));if isnan(nn{1});warning('%s didn''t work',n);end%get data from webpage
end%if
if ~isnan(nn{1});break;end%found one!
end%for
return

function [bo bi]=getsbmb(n)
%get small or major body structure
if isnumeric(n);if n<1e6;bi='mb';else;bi='sb';end%convert number input to sb mb
else;bi=n;end%should be 'mb' or 'sb'
bo=params(bi);if isempty(bo)%read from saved data
bdf=[params('bdir') 'bodydata.mat'];
if exist(bdf,'file');bo=load(bdf,bi);bo=bo.(bi);%read from file
else;bo.names='';bo.numbers=[];end%noob
params(bi,bo);end%save
return

function x=putsbmb(n,f,v)
%write data to structure
if ischar(f);f={f};v={v};end%convert to cell
[x bi]=getsbmb(n(1));%get structure
bb=find(n(1)==x.numbers);if isempty(bb)&&numel(n)>1;n=n(2);bb=find(n==x.numbers);else;n=n(1);end
if isempty(bb);bb=numel(x.numbers)+1;f=['numbers' f];v=[n v];end%see if body exists in field
for ii=1:numel(f)%different fields
if ischar(v{ii});if ~isfield(x,f{ii});x.(f{ii})='';end%character, initialize if necessary
%pad with spaces, then write
if numel(v{ii})<size(x.(f{ii}),2);v{ii}(size(x.(f{ii}),2))=' ';else;x.(f{ii})(bb,numel(v{ii}))=' ';end;x.(f{ii})(bb,:)=v{ii};
else;if ~isfield(x,f{ii});x.(f{ii})=[];end%number, initialize if necessary
%write, fill with nans
nn=size(x.(f{ii}),2)+1;x.(f{ii})(:,bb)=v{ii};if bb>nn;x.(f{ii})(:,nn:bb-1)=nan;end;end
end%numel(f)
params([bi 'mod'],true);params(bi,x);%save
return

function in=uniqstr(xni,bi,nb)
%find a unique match of bi in xni
xn=xni';%bi=strtrim(lower(bi));xn=char(' ',lower(xni)');
nxn=size(xn,1);xn=xn(:)';%nxn is width of each xn entry
if nargin<3;nb=0;end%1 for only whole word matches, 2 for only fragment matches, 0 for either
bi=regexprep(bi,'([()])','\\$1');%replace ( with \(
if nb<2;in=regexpi(xn,['\<' bi '\>']);else;in=[];end%in=strfind(xn,[' ' bi ' ']);%whole word
if numel(in)>1%multiple entries with same word
    %in is strfind index, xm tracks current list of multiple matches
    xm=xni(ceil(in/nxn),:);in=in(strcmpi(bi,cellstr(xni(ceil(in/nxn),:))));%exact match
    if numel(in)>1;matches=xni(ceil(in/nxn),:),warning('identical entries');in=in(1);
    elseif isempty(in);matches=xm,warning('ambiguous string match');end
elseif isempty(in)&&(nb~=1)
    in=regexpi(xn,bi);%fragment anywhere
    if numel(in)>1
        in1=in;xm=xni(ceil(in/nxn),:);in=in(in==floor(in/nxn)*nxn+1);%beginning of entry
        if isempty(in);in=in(regexp(xn(in-1),'\W'));%beginning of word
            if numel(in)>1;xm=xni(ceil(in/nxn),:);elseif isempty(in);in=in1;end
        elseif numel(in)>1;xm=xni(ceil(in/nxn),:);
        end
%         in1=in;xm=xni(ceil(in/nxn),:);in=in(isspace(xn(in-1)));%beginning of word
%         if numel(in)>1
%             in1=in;xm=xni(ceil(in/nxn),:);in=in(in==floor(in/nxn)*nxn+2);%beginning of entry
%             if numel(in)>1;xm=xni(ceil(in/nxn),:);elseif isempty(in);in=in1;end
%         elseif isempty(in);in=in1;
%         end
    end
    if numel(in)>1;matches=xm,warning('ambiguous string match');in=in(1);end
end
in=ceil(in/nxn);
return

function cb=cbfun(b)
if numel(b)==1;
    cb=10;
    if imag(b);cb=real(b);
    elseif b>10&&b<1000&&mod(b,100)~=99;cb=floor(b/100)*100+99;
    elseif b>1e3&&b<1e4;cb=floor(b/10);
    end

else
    cb=repmat(10,size(b));
    ii=b>10&b<1000&mod(real(b),100)~=99;cb(ii)=floor(b(ii)/100)*100+99;
    ii=find(imag(b));cb(ii)=real(b(ii));
    ii=b>1e3&b<1e4;cb(ii)=floor(b(ii)/10);
end
return


%function [d e]=curl(url);e=false;try;d=webread(url);catch;e=true;end;return
function [d e]=curl(url,fout)
%use curl, "urlread" is slow over vpn, maybe adjust proxy?
%could also use wget -O- 
%http://curl.haxx.se/dlwiz/?type=bin&os=Win32&flav=-&ver=-
persistent cdir
if isempty(cdir);%see if curl exists in boddat dir or on path
cdir=mfilename('fullpath');ii=regexp(cdir,'/|\');cdir=cdir(1:ii(end));
if ~exist([cdir 'curl.exe'])&&~exist([cdir 'curl']);cdir=' ';[e d]=system('which curl');
if e;if ispc;dn='nul';else;dn='/dev/null';end
    [e d]=system([' curl -s "https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=2142" > ' dn]);%run it & see what happens
if e;error('curl is not a recognized system command from boddat.m\nTry copying a curl executable to %s',cdir);
end;end;end
end%isempty cdir 
%read the url, write to file and optionally read contents
if nargin<2;fout=tempname;end
if ispc;[e d]=system(['"', cdir, 'curl" "', url, '" > ', fout]);
else;[e d]=system([cdir 'curl "' url '" > ' fout]);end
if nargin<2;if ~e;f=fopen(fout);d=fscanf(f,'%c');fclose(f);end;delete(fout);end
return


function o=params(p,v)
%variables global to function
%make p persistent and write v to p, then retrieve p from any function via o 
persistent vs
if nargin==2;vs.(p)=v;elseif isfield(vs,p);o=vs.(p);else;o=[];end
return
