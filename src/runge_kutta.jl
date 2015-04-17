# Butcher Tableaus
# see Hairer & Wanner 1992, p. 134, 166

abstract Tableau{Name, S, T<:FloatingPoint}
# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# assumes fields
# order::(Int...) # order of the method(s)
# a::Matrix{T}  # SxS-1 matrix
# b::Matrix{T}  # 1 or 2 x S matrix (fixed step/ adaptive)
# c::Vector{T}  # S
#
# For a tableau:
#  c1  | a_11   ....   a_1s
#  .   | a_21 .          .           
#  .   | a_31     .      .           
#  .   | ....         .  .           
#  c_s | a_s1  ....... a_ss
# -----+--------------------
#      | b_1     ...   b_s 
#      | b'_1    ...   b'_s 

Base.eltype{N,S,T}(b::Tableau{N,S,T}) = T

# #@doc "Test whether it's an explicit method"->
isexplicit(b::Tableau) = istril(b.a)
isadaptive(b::Tableau) = size(b.b, 1)==2
order(b::Tableau) = b.order

# The advantage of each having its own type, makes it possible to have
# specialized methods for a particular tablau
immutable TableauExplicit{Name, S, T} <: Tableau{Name, S, T}
    order::(Int...) # the order of the methods
    a::Matrix{T}    # Todo: make this into a LowerTriangular matrix in 0.4
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function TableauExplicit(order,a,b,c)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert c[1]==0
        @assert istril(a)
        @assert S==length(c)==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        new(order,a,b,c)
    end
end
function TableauExplicit{T}(name::Symbol, order::(Int...),
                   a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    TableauExplicit{name,length(c),T}(order, a, b, c)
end
function TableauExplicit(name::Symbol, order::(Int...), T::Type,
                   a::Matrix, b::Matrix, c::Vector)
    TableauExplicit{name,length(c),T}(order, convert(Matrix{T},a),
                                        convert(Matrix{T},b), convert(Vector{T},c) )
end


## Tableaus for explicit methods
# Fixed step:
bt_feuler = TableauExplicit(:feuler,(1,), Float64,
                             zeros(Int,1,1),
                            [1]',
                            [0]
                             )
bt_midpoint = TableauExplicit(:midpoint,(2,), Float64,
                               [0  0
                                .5  0],
                              [0, 1]',
                              [0, .5]
                              )
bt_heun = TableauExplicit(:heun,(2,), Float64,
                           [0  0
                            1  0],
                          [1//2, 1//2]',
                          [0, 1])

bt_rk4 = TableauExplicit(:rk4,(4,),Float64,
                          [0 0 0 0
                           1//2 0 0 0
                           0 1//2 0 0
                           0 0 1 0],
                         [1//6, 1//3, 1//3, 1//6]',
                         [0, 1//2, 1//2, 1])
                         
# Adaptive step:

# Fehlberg https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
bt_rk45 = TableauExplicit(:fehlberg,(4,5),Float64,
                          [  0           0           0            0         0     0
                             1//4        0           0            0         0     0
                             3//32       9//32       0            0         0     0
                          1932//2197 -7200//2197  7296//2197      0         0     0
                           439//216     -8        3680//513    -845//4104   0     0
                            -8//27       2       -3544//2565   1859//4104 -11//40 0 ],
                           [25//216      0        1408//2565   2197//4104  -1//5  0
                            16//135      0        6656//12825 28561//56430 -9//50 2//55],
                            [0,          1//4,       3//8,       12//13,    1,    1//2])

                    

# Dormand-Prince https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
bt_dopri5 = TableauExplicit(:dopri, (5,4), Float64,
                     [0   0 0 0 0 0 0
                      1//5 0 0 0 0 0 0
                      3//40     9//40 0 0 0 0 0 
                      44//45     -56//15     32//9 0 0 0 0
                      19372//6561     -25360//2187     64448//6561     -212//729 0 0 0
                      9017//3168     -355//33     46732//5247     49//176     -5103//18656 0 0
                      35//384     0     500//1113     125//192     -2187//6784     11//84  0],
                     [35//384         0      500//1113      125//192      -2187//6784         11//84      0
                      5179//57600     0     7571//16695     393//640     -92097//339200     187//2100     1//40],
                     [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
                     )

# Fehlberg 7(8) coefficients
# Values from pag. 65, Fehlberg, Erwin. "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control".
# National Aeronautics and Space Administration.
bt_feh78 = TableauExplicit(:feh78, (7,8), Float64,
                            [     0      0      0       0        0         0       0       0     0      0    0 0 0
                                  2//27   0      0       0        0         0       0       0     0      0    0 0 0
                                  1//36   1//12   0       0        0         0       0       0     0      0    0 0 0
                                  1//24   0      1//8     0        0         0       0       0     0      0    0 0 0
                                  5//12   0    -25//16   25//16     0         0       0       0     0      0    0 0 0
                                  1//20   0      0       1//4      1//5       0       0       0     0      0    0 0 0
                                -25//108  0      0     125//108  -65//27    125//54    0       0     0      0    0 0 0
                                 31//300  0      0       0       61//225    -2//9    13//900   0     0      0    0 0 0
                                  2      0      0     -53//6    704//45   -107//9    67//90    3     0      0    0 0 0
                                -91//108  0      0      23//108 -976//135   311//54  -19//60   17//6  -1//12   0    0 0 0
                               2383//4100 0      0    -341//164 4496//1025 -301//82 2133//4100 45//82 45//164 18//41 0 0 0
                                  3//205  0      0       0        0        -6//41   -3//205  -3//41  3//41   6//41 0 0 0
                              -1777//4100 0      0    -341//164 4496//1025 -289//82 2193//4100 51//82 33//164 12//41 0 1 0],
                              [41//840 0 0 0 0 34//105 9//35 9//35 9//280 9//280 41//840 0 0
                               0 0 0 0 0 34//105 9//35 9//35 9//280 9//280 0     41//840 41//840],
                               [0,    2//27, 1//9, 1//6 , 5//12, 1//2 , 5//6 , 1//6 , 2//3 , 1//3 , 1 , 0, 1]
                            )



# make a Runge-Kutta method for a given Butcher tableau.  Follows
# Hairer & Wanner 1992 p.134, p.165-169
export oderk_fixed, oderk_adapt

# to put ys into the vector of vector format:
function transformys{T}(ys::Array{T})
    if size(ys,1)==1
        squeeze(ys,1)
    elseif length(size(ys))!=1
        Vector{T}[ys[:,i] for i=1:size(ys,2)]
    else
        ys
    end
end

# returns if t1 is before or at t2 (in direction of integration tdir)
before(t1, t2, tdir) = tdir*t1 < tdir*t2

# Fixed step Runge-Kutta method
# TODO: iterator method
function oderk_fixed{N,S,T}(fn, y0, tspan::AbstractVector,
                            btab::TableauExplicit{N,S,T})
    dof = length(y0)
    tsteps = length(tspan)
    ys = Array(T, dof, tsteps)
    ys[:,1] = y0'
    tspan = convert(Vector{T}, tspan)
    # work arrays:
    ks = zeros(T, dof, S)
    ytmp = zeros(T, dof)
    # time stepping:
    for i=1:length(tspan)-1
        dt = tspan[i+1]-tspan[i]
        ys[:,i+1] = ys[:,i]
        for s=1:S
            ytmp[:] = ys[:,i]
            calc_next_k!(ks, ytmp, s, fn, tspan[i], dt, dof, btab)
            for d=1:dof
                ys[d,i+1] += dt * btab.b[s]*ks[d,s]
            end
        end
    end
    return tspan, transformys(ys)
end
# calculates k[s]
function calc_next_k!{N,S,T}(ks::Matrix, ytmp::Vector, s, fn, t, dt, dof, btab::TableauExplicit{N,S,T})
    # Calculates the next ks and puts it into ks[
    # - ytmp needs to be set to the current y on entry.
    #
    # - ks and ytmp are modified inside this function.
    #
    # TODO:
    # - allow reuse as in Dormand-Price
    # - calculate all k in one go
    # - k as vector of vector
    for ss=1:s-1
        for d=1:dof
            #@show dt*ks[d,ss] * btab.a[s,ss]
            ytmp[d] += dt * ks[d,ss] * btab.a[s,ss]
        end
    end
    ks[:,s] = fn(t + btab.c[s]*dt, ytmp) # TODO update k to be vector of vector k[s][:] = ...
    nothing
end

ode4_v2(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_rk4)
ode1_euler(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_feuler)
ode2_midpoint(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_midpoint)
ode2_heun(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_heun)


# Adaptive Runge-Kutta method
function oderk_adapt{N,S,T}(fn, y0, tspan, btab::TableauExplicit{N,S,T};
                     reltol = 1.0e-5, abstol = 1.0e-8,
                     norm=Base.norm,
                     minstep=abs(tspan[end] - tspan[1])/1e9,
                     maxstep=abs(tspan[end] - tspan[1])/2.5,
                     initstep=0.)

    !isadaptive(btab) && error("Can only use this solver with an adpative Butcher table")

    # TODO: fix magic numbers
    const large = 1.0e5
    facmax = large
    order = minimum(btab.order)

    ## Initialization
    dof = length(y0)
    tspan = convert(Vector{T}, tspan)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # If tspan is a more than a length two vector: return solution at
    # those points only
    tstepsgiven = length(tspan)>2

    if tstepsgiven
        nsteps = length(tspan)
        ys = zeros(T, dof, nsteps)
        ys[:,1] = y0
        y0 = ys[:,1] # now of right type
        iter = 1 # the index into tspan
    else
        ys = zeros(T, dof)
        ys[:] = y0
        y0 = ys[:] # now of right type
        tspan = [tstart]
        # TODO check that this makes a difference:
        sizehint_size = 100
        sizehint(tspan, sizehint_size)
        sizehint(ys, sizehint_size)
    end
    # Initialize work arrays:
    y      = copy(y0)  # the y at current time t
    ytrial = zeros(T, dof)
    yerr   = zeros(T, dof)
    ks     = zeros(T, dof, S) 
    ytmp   = zeros(T, dof)
    # Initialize time:
    tdir = sign(tend-tstart)
    if initstep == 0.
        # initial guess at a step size
        dt, f0 = hinit(fn, y0, tstart, tdir, order, reltol, abstol)
    else
        dt = sign(initstep)==tdir ? initstep : error("initstep has wrong sign.")
    end
    # Diagnostics
    dts = Float64[]
    errs = Float64[]
    steps = [0,0]  # [accepted, rejected]
    laststep = false
    while true
        # do one step
        ytrial[:] = 0
        yerr[:] = 0
        for s=1:S
            ytmp[:] = y # needed for calc_next_k!
            calc_next_k!(ks, ytmp, s, fn, t, dt, dof, btab)
            for d=1:dof
                ytrial[d] += btab.b[1,s]*ks[d,s]
                yerr[d]   += btab.b[2,s]*ks[d,s]
            end
        end
        for d=1:dof
            yerr[d]   = dt * (ytrial[d]-yerr[d])
            ytrial[d] = y[d] + dt * ytrial[d]
        end
        
        # Check error and find a new step size:
        err, newdt = stepsize_hw92(dt, tdir, y, ytrial, yerr,
                                   abstol, reltol, order, facmax, dof,
                                   maxstep)

        if err<=1.0
            # diagnostics
            steps[1] +=1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            if tstepsgiven
                # interpolate onto given output points
                f0 = ks[:,1]
                f1 = fn(t+dt, ytrial) # TODO: this could be used in next step
                while iter<nsteps && tdir*tspan[iter+1]<=tdir*(t+dt)  # output at all new times which are ≤ t+dt
                    iter += 1
                    # ys[:,iter] = hermite_interp(tspan[iter], t, dt,
                    #                             y, ytrial, f0, f1)
                    hermite_interp!(ys, iter, tspan[iter], t, dt,
                                                y, ytrial, f0, f1)
                end
            else
                # output every step
                append!(ys, ytrial)
                push!(tspan, t+dt)
            end
            # Break if this was the last step:
            if laststep
                break
            end

            # Swap bindings of y and ytrial, avoids one copy
            tmp = y
            y = ytrial
            ytrial = tmp

            # Update t to the time at the end of current step:
            t += dt
            dt = newdt

            # Hit end point exactly if within 1% of end:
            if !before(t+dt*1.01, tend, tdir)
                dt = tend-t
                laststep = true # next step is the last, if it succeeds
            end
            facmax = large
        elseif abs(dt)<minstep  # minimum step size reached
            @show length(ys), t, dt
            println("Warning: dt < minstep.  Stopping.")
            break
        else # redo step with smaller dt
            laststep = false
            steps[2] +=1 
            dt = newdt
            facmax = 1.0 # forbids dt increases in the next step
                         # TODO: make that several steps
        end
    end
    if tstepsgiven
        ys = reshape(ys, dof, nsteps)
    end
    #    xerrs = reshape(xerrs, dof, length(dts))
    @show steps
    return tspan, transformys(ys)
end
ode45_v2(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk45; kwargs...)
ode54_v2(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_dopri5; kwargs...)
ode78_v2(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_feh78; kwargs...)

# Helper functions:
function stepsize_hw92(dt, tdir, x0, xtrial, xerr, abstol, reltol, order,
                       facmax, dof, maxstep)
    # Estimates new best step size following
    # Hairer & Wanner 1992, p167 (with some modifications)

    # TODO: parameters to lift out of this function 
    fac = [0.8, 0.9, 0.25^(1/(order+1)), 0.38^(1/(order+1))][1]
    facmax_def = 2.0 # maximal step size increase. 1.5-5
    facmax = min(facmax_def, facmax)
    facmin = 1./facmax_def  # maximal step size decrease. ?

    any(isnan(xtrial)) && return 10., dt*facmin

    # figure out a new step size
    tol = abstol + max(abs(x0), abs(xtrial)).*reltol # 4.10
    err = sqrt(1/dof * sum((xerr./tol).^2) )     # 4.11
    #err = norm(xerr./tol,Inf )     # 4.11

    newdt = dt * min(facmax, max(facmin, fac*(1/err)^(1/(order+1))))
    return err, tdir*min(tdir*newdt, maxstep)
end

function stepsize_wh96(dt, tdir, x0, xtrial, xerr, abstol, reltol, order,
                       facmax, dof, maxstep)
    # Estimates new best step size following
    # Wanner & Hairer 1996, p.112 (with some modifications)

    # TODO: lift parameters out of here:
    c1 = min(facmax, 6.)  # max step size increase
    c2 = 0.2 # max step size decrease
    c3 = 0.8 # safety factor, lower==saver
    
    if any(isnan(xtrial))
        # reduce step size by max if solution contains NaN
        return dt*c2
    end

    
    tol = abstol + max(abs(x0), abs(xtrial)).*reltol
    err = maximum(abs(xerr./tol))
    # e = 0.9*(1/err)^(1/(order+1))
    # @show e

    # TODO maxstep
    return err, min( dt * min(c1, max(c2, c3*(1/err)^(1/(order+1))) ), maxstep)
end

function stepsize_ODE(dt, tdir, x0, xtrial, xerr, abstol, reltol, order,
                      facmax, dof, maxstep)
    # estimate the local truncation error
    gamma1 = xerr
    
    # Estimate the error and the acceptable error
    delta = norm(gamma1, Inf)              # actual error
    tau   = max(reltol*norm(x0,Inf),abstol) # allowable error
    err = delta/tau
    newdt = tdir*min(maxstep, 0.8*abs(dt)*(1/err)^(1/(order+1)))
    return err, newdt
end

# For dense output see Hairer & Wanner p.190 using Hermite interpolation
function hermite_interp(tquery,t,dt,y0,y1,f0,f1)
    # f_0 = f(x_0 , y_0) , f_1 = f(x_0 + h, y_1 )
    # this is O(3). TODO for higher order.
    theta = (tquery-t)/dt
    return (1-theta)*y0 + theta*y1 + theta*(theta-1) *
           ((1-2*theta)*(y1-y0) + (theta-1)*dt*f0 + theta*dt*f1)
end

function hermite_interp!(ys, iter, tquery,t,dt,y0,y1,f0,f1)
    # f_0 = f(x_0 , y_0) , f_1 = f(x_0 + h, y_1 )
    # this is O(3). TODO for higher order.
    theta = (tquery-t)/dt
    for i=1:length(y0)
        ys[i,iter] = ((1-theta)*y0[i] + theta*y1[i] + theta*(theta-1) *
                      ((1-2*theta)*(y1[i]-y0[i]) + (theta-1)*dt*f0[i] + theta*dt*f1[i]) )
    end
    nothing
end
