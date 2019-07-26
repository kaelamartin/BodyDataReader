"""
	pp = RVspline6!(X, pp)

5th order spline fitting the position and velocity of segment endpoints when there
is continuous acceleration and jerk for a specific body. Creates 5th, 4th, 3rd,
and 2nd order splines.
"""
function RVspline6!(X::AbstractArray{T},pp::Dict{String,Array{Float64}}) where T
    #5th order spline fiting R,V at segment endpoints, continuous A,J
    #nonuniform step ok
    t = view(X, 1, :)
    R = view(X, 2:4, :)
    X[5:7,:] *= 86400
    V = view(X, 5:7, :)
    n = size(X, 2) - 1
    R_ = Array{Float64}(undef, 3,n); V_ = Array{Float64}(undef, 3,n); dt = Array{Float64}(undef, n)
    b_ = Array{Float64}(undef, 3,n)
    for ii in 1:n
        dt[ii] = t[ii+1]-t[ii]
        for jj in 1:3
            R_[jj,ii] = (R[jj,ii]+V[jj,ii]*dt[ii]-R[jj,ii+1])/dt[ii]^2
            V_[jj,ii]=(V[jj,ii]-V[jj,ii+1])/dt[ii]/2.
        end
    end
    b = Array{Float64}(undef, 3,n-1)
    s = Array{Float64}(undef, 3*(n-1))
    for ii in 2:n
        for jj in 1:3
            b[jj,ii-1]=(12 .*V_[jj,ii-1]-10 .*R_[jj,ii-1])/dt[ii-1]-
            (8 .*V_[jj,ii]-10 .*R_[jj,ii])/dt[ii]
        end
        s[ii-1]= 1 ./dt[ii-1]
        s[ii+n-2] = -3 ./dt[ii-1]-3 ./dt[ii]
        s[ii+2*n-3] = 1 ./dt[ii]
    end
    #s*[c2(i-1);c2(i);c2(i+1)]=b, tridiagonal
    ii = 2:n
    jj=[ii .- 1; ii; ii .+ 1] #sparsity pattern x,y,z
    ii=[ii; ii; ii]
    b_ = zeros(3) #continuous 4th der at 2
    for kk in 1:3
        b_[kk]=(-15 .*R_[kk,1]+16 .*V_[kk,1])/dt[1]^2-
        (15 .*R_[kk,2]-14 .*V_[kk,2])/dt[2]^2
    end
    s_=[2 ./dt[1]^2; -3 ./dt[1]^2+3 ./dt[2]^2; -2 ./dt[2]^2]
    b=[b_ b]
    s=[s_; s]

    ii=[1; 1; 1; ii]
    jj=[1; 2; 3; jj]
    b_ = zeros(3) #continuous 4th der at n
    for kk in 1:3
        b_[kk]=(-15 .*R_[kk,n-1]+16*V_[kk,n-1])/dt[n-1]^2-
        (15 .*R_[kk,n]-14 .*V_[kk,n])/dt[n]^2
    end
    s_=[2 ./dt[n-1]^2; -3 ./dt[n-1]^2+3 ./dt[n]^2; -2 ./dt[n]^2]
    b=[b b_]
    s=[s; s_]
    ii=[ii; n+1; n+1; n+1]
    jj=[jj; n-1; n; n+1]

    s=sparse(jj,ii,s) #solve for c2
    c2=b/s
    c22=c2[:,2:n+1]
    c2=c2[:,1:n]
    #solve for coefficients
    c5 = Array{Float64}(undef, 3,n); c3 = Array{Float64}(undef, 3,n)
    c4 = Array{Float64}(undef, 3,n)
    c0 = Array{Float64}(undef, 3,n); c1 = Array{Float64}(undef, 3,n)
    for ii in 1:n
        for jj in 1:3
            c5[jj,ii] = (6 .*(V_[jj,ii]-R_[jj,ii])-c2[jj,ii]+c22[jj,ii])/dt[ii]^3
            c5dt3 = c5[jj,ii]*dt[ii]^3
            c3[jj,ii]=(2 .*V_[jj,ii]-4 .*R_[jj,ii]-2 .*c2[jj,ii]+c5dt3)/dt[ii]
            c4[jj,ii]=(3 .*R_[jj,ii]-2 .*V_[jj,ii]+c2[jj,ii]-2 .*c5dt3)/dt[ii]^2
            c0[jj,ii] = R[jj,ii]
            c1[jj,ii] = V[jj,ii]
        end
    end
    pp["breaks"] = t

    pp["dcm6"] = Float64[]
    pp["dcm5"] = Float64[]
    pp["dcm4"] = Float64[]
    pp["dcm3"] = Float64[]
    #derivatives
    pp["coefs6"] = zeros(length(c5),6)
    pp["coefs5"] = zeros(length(c5),5)
    pp["coefs4"] = zeros(length(c5),4)
    pp["coefs3"] = zeros(length(c5),3)
    for ii in 1:length(c5)
        pp["coefs6"][ii,1] = c5[ii]
        pp["coefs6"][ii,2] = c4[ii]
        pp["coefs6"][ii,3] = c3[ii]
        pp["coefs6"][ii,4] = c2[ii]
        pp["coefs6"][ii,5] = c1[ii]
        pp["coefs6"][ii,6] = c0[ii]

        pp["coefs5"][ii,1] =  5 .*c5[ii]/86400.
        pp["coefs5"][ii,2] =  4 .*c4[ii]/86400.
        pp["coefs5"][ii,3] =  3 .*c3[ii]/86400.
        pp["coefs5"][ii,4] =  2 .*c2[ii]/86400.
        pp["coefs5"][ii,5] =  c1[ii]/86400.

        pp["coefs4"][ii,1] = 20 .*c5[ii]/86400 .^2
        pp["coefs4"][ii,2] = 12 .*c4[ii]/86400 .^2
        pp["coefs4"][ii,3] = 6 .*c3[ii]/86400 .^2
        pp["coefs4"][ii,4] = 2 .*c2[ii]/86400 .^2

        pp["coefs3"][ii,1] = 60 .*c5[ii]/86400 .^3
        pp["coefs3"][ii,2] = 24 .*c4[ii]/86400 .^3
        pp["coefs3"][ii,3] = 6 .*c3[ii]/86400 .^3
    end
    return pp
end

"""
	RVAspline!(X, pp)

5th order spline that uses position, velocity, and acceleration to approximate
acceleration with central difference of position and velocity. Creates a 6th order
spline.
"""
function RVAspline!(X::AbstractArray{T},pp::Dict{String,Array{Float64}}) where T
    #5th order spline using R,V,A at beginning and end of segments
    #approximates A with central diff of R and V
    #assumes uniform time step
    t = view(X, 1, :)
    R = view(X, 2:4, :)
    X[5:7,:] *= 86400
    V = view(X, 5:7, :)
    n = size(X, 2)
    dt = Array{Float64}(undef, n-1) #uniform timestep, scale V
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

    c2 = Array{Float64}(undef, 3, n-2)#central diff for A (c2 = coeff on t^2 term)
    for ii in 1:n-2
        for jj in 1:3
            c2[jj,ii]=R[jj,ii]-2 .*R[jj,ii+1]+R[jj,ii+2]+(V[jj,ii]-V[jj,ii+2])/4
        end
    end
    c21 = Array{Float64}(undef, 3)#initial A (uses first 3 R,V)
    c2n = Array{Float64}(undef, 3)#final A (uses last 3 R,V)
    for jj in 1:3
        c21[jj]=-23 ./4 .*R[jj,1]+4 .*R[jj,2]+7 ./4 .*R[jj,3]-3 .*V[jj,1]-
        4 .*V[jj,2]-V[jj,3]/2.
        c2n[jj]=7 ./4 .*R[jj,n-2]+4 .*R[jj,n-1]-23 ./4 .*R[jj,n]+V[jj,n-2]/2 .+
        4 .*V[jj,n-1]+3 .*V[jj,n]
    end
    # coeff from initial R,V,A and temp final R,V,A on segs
    c0 = Array{Float64}(undef, 3, n-1); c1 = Array{Float64}(undef, 3, n-1)
    R2 = Array{Float64}(undef, 3, n-1); V2 = Array{Float64}(undef, 3, n-1)

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
    c3=-10 .*R2-6 .*c1-4 .*V2-3 .*c2+A2
    c4= 15 .*R2+8 .*c1+7 .*V2+3 .*c2-2 .*A2
    c5=-6 .*R2-3 .*c1-3 .*V2-c2+A2
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
    pp["dcm6"] = Float64[]
    pp["coefs6"] = c
    pp["dc"] = Float64[]
    d=[5 0 0 0 0; 0 4 0 0 0; 0 0 3 0 0; 0 0 0 2 0; 0 0 0 0 1.; 0 0 0 0 0]
    pp["coefs5"] = c*d/86400
    pp["dcm5"] = Float64[]
end

"""
	spline!(x, y, pp)

Cubic spline that creates third order polynomials between timed segments
segments (interpolation for data) for pole and pm calculations. Used under the
orient function.
"""
function spline!(x::AbstractArray{P},y::AbstractArray{P},
                pp::Dict{String,Array{Float64}}) where P # Cubic spline with x1 <= x2 <= .. <= xn
    # x is a 1 x n or n x 1 array
    # y is a 2 X n, n X 2, 1 x n, n x 1 array

    (size(x,2)>size(x,1)) && (x = x')
    dx = diff(x)

    # Find y difference
    (size(y,2)>size(y,1)) && (y = y')
    dy = diff(y, dims=1)
    n = length(x)

    # values in tridiagonal
    av = zeros(n-1); bv = zeros(n); cv = zeros(n-1)
    cv[1] = x[3]-x[1]
    bv[1] = dx[2]; bv[n] = dx[n-2]
    av[n-1] = x[n] - x[n-2]

    for ii in 1:n-2
        cv[ii+1] = dx[ii]
        bv[ii+1] = 2*(dx[ii]+dx[ii+1])
        av[ii] = dx[ii+1]
    end
    A = spdiagm(-1 => av, 0 => bv, 1 => cv)

    # Create vector so that A D = C
    if size(y,2) == 2
        C = zeros(n,2)
        C[1,1] = ((dx[1]+2*(x[3]-x[1]))*dx[2]*dy[1,1]/dx[1]+
            dx[1]^2*dy[2,1]/dx[2])/(x[3]-x[1])
        C[1,2] = ((dx[1]+2*(x[3]-x[1]))*dx[2]*dy[1,2]/dx[1]+
            dx[1]^2*dy[2,2]/dx[2])/(x[3]-x[1])
        C[end,1] = (dx[n-1]^2*dy[n-2,1]/dx[n-2]+
            (2*(x[n]-x[n-2])+dx[n-1])*dx[n-2]*dy[n-1,1]/dx[n-1])/(x[n]-x[n-2])
        C[end,2] = (dx[n-1]^2*dy[n-2,2]/dx[n-2]+
            (2*(x[n]-x[n-2])+dx[n-1])*dx[n-2]*dy[n-1,2]/dx[n-1])/(x[n]-x[n-2])
        for ii in 2:n-1
            C[ii,1] = 3*(dx[ii]*dy[ii-1,1]/dx[ii-1]+dx[ii-1]*dy[ii,1]/dx[ii])
            C[ii,2] = 3*(dx[ii]*dy[ii-1,2]/dx[ii-1]+dx[ii-1]*dy[ii,2]/dx[ii])
        end

        D = \(A,C)
        # Determine the coefficients of the third order polynomial
        coefs = zeros(2*(n-1),4)
        for jj in 1:n-1
            coefs[jj,4] = y[jj,1]
            coefs[n-1+jj,4] = y[jj,2]
            coefs[jj,3] = D[jj,1]
            coefs[n-1+jj,3] = D[jj,2]
            coefs[jj,2] = 3*dy[jj,1]/dx[jj]^2 - 2*D[jj,1]/dx[jj] - D[jj+1,1]/dx[jj]
            coefs[n-1+jj,2] = 3*dy[jj,2]/dx[jj]^2 - 2*D[jj,2]/dx[jj] - D[jj+1,2]/dx[jj]
            coefs[jj,1] = -2*dy[jj,1]/dx[jj]^3 + D[jj,1]/dx[jj]^2 + D[jj+1,1]/dx[jj]^2
            coefs[n-1+jj,1] = -2*dy[jj,2]/dx[jj]^3 + D[jj,2]/dx[jj]^2 + D[jj+1,2]/dx[jj]^2
        end
    elseif size(y,2) == 1
        C = zeros(n)
        C[1] = ((dx[1]+2*(x[3]-x[1]))*dx[2]*dy[1,1]/dx[1]+
            dx[1]^2*dy[2,1]/dx[2])/(x[3]-x[1])
        C[end] = (dx[n-1]^2*dy[n-2,1]/dx[n-2]+
            (2*(x[n]-x[n-2])+dx[n-1])*dx[n-2]*dy[n-1,1]/dx[n-1])/(x[n]-x[n-2])
        for ii in 2:n-1
            C[ii] = 3*(dx[ii]*dy[ii-1,1]/dx[ii-1]+dx[ii-1]*dy[ii,1]/dx[ii])
        end
        D = \(A,C)

        # Determine the coefficients of the third order polynomial
        coefs = zeros(n-1,4)
        for jj in 1:n-1
            coefs[jj,4] = y[jj]
            coefs[jj,3] = D[jj]
            coefs[jj,2] = 3*dy[jj,1]/dx[jj]^2 - 2*D[jj]/dx[jj] - D[jj+1]/dx[jj]
            coefs[jj,1] = -2*dy[jj,1]/dx[jj]^3 + D[jj]/dx[jj]^2 + D[jj+1]/dx[jj]^2
        end
    end

    # send out the results
    pp["breaks"] = x
    pp["coefs4"] = coefs
end

"""
	pval!(pp, xx, order, v)


Evaluates a polynomial by getting the coefficients from the spline and evaluating
the polynomial over many steps. v has to be a zero input.
"""
function pval!(pp::Dict{String,Array{Float64}},xx::AbstractArray{T},
                order::Int64,v::AbstractArray{T}) where T
    # Assumes that xx is sorted from smallest to largest
    # pp is dictionary
    # xx is a 1 x n or n x 1
    # Uses Horner's method to calculate the polynomial
    xs = vec(xx)
    lx = length(xs)
    b = pp["breaks"]
    lb = length(b)
    c = pp[string("coefs",order)]
    k = order
    (size(v,2)==lx) ? (sizev = size(v,1)) : (sizev = size(v,2))

    index = zeros(Int64, lx)

    for ii in 1:lx
        index[ii] = searchsortedlast(b,real(xs[ii]))
        (index[ii] >= lb) && (index[ii]= lb-1)
        (isnan(xs[ii])) && (index[ii] = 1)  # Remove NaNs
        (iszero(index[ii])) && (index[ii] = 1)
    end

    for ii in 1:lx
        xT = xs[ii]-b[index[ii]]
        for jj in 1:sizev
            iT = sizev*(index[ii]-1)+jj
            for kk in 1:k
                v[jj,ii] = xT*v[jj,ii] + c[iT,kk]
            end
        end
    end
end

function pvalorient!(pp::Dict{String,Array{Float64}},xx::AbstractArray{T},
                order::Int64,v::AbstractArray{T}) where T
    # For evaluationg pole and pm splines
    # pm and pole splines are ordered by first element wrt to time and then second
    # pp is dictionary
    # xx is a 1 x n or n x 1
    # Uses Horner's method to calculate the polynomial
    xs = vec(xx)
    lx = length(xs)
    b = pp["breaks"]
    lb = length(b)
    c = pp[string("coefs",order)]
    k = order
    (size(v,2)==lx) ? (sizev = size(v,1)) : (sizev = size(v,2))

    index = zeros(Int64,lx)

    for ii in 1:lx
        index[ii] = searchsortedlast(b,real(xs[ii]))
        (index[ii] >= lb) && (index[ii] = lb-1)
        (isnan(xs[ii])) && (index[ii] = 1)  # Remove NaNs
        (iszero(index[ii])) && (index[ii] = 1)
    end

    for ii in 1:lx
        xT = xs[ii]-b[index[ii]]
        for jj in 1:sizev
            iT = (lb-1)*(jj-1)+index[ii]
            for kk in 1:k
                v[jj,ii] = xT*v[jj,ii] + c[iT,kk]
            end
        end
    end
end
