# Butcher Tableaus
# see Hairer & Wanner 1992, p. 134, 166

using Compat

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
#      | b_1     ...   b_s   this is the one used for stepping
#      | b'_1    ...   b'_s  this is the one used for error-checking

Base.eltype{N,S,T}(b::Tableau{N,S,T}) = T

# Test whether it's an explicit method
isexplicit(b::Tableau) = istril(b.a)
isadaptive(b::Tableau) = size(b.b, 1)==2
order(b::Tableau) = b.order

# The advantage of each having its own type, makes it possible to have
# specialized methods for a particular tablau
immutable TableauRKExplicit{Name, S, T} <: Tableau{Name, S, T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function TableauRKExplicit(order,a,b,c)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert c[1]==0
        @assert istril(a)
        @assert S==length(c)==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        new(order,a,b,c)
    end
end
function TableauRKExplicit{T}(name::Symbol, order::(Int...),
                   a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    TableauRKExplicit{name,length(c),T}(order, a, b, c)
end
function TableauRKExplicit(name::Symbol, order::(Int...), T::Type,
                   a::Matrix, b::Matrix, c::Vector)
    TableauRKExplicit{name,length(c),T}(order, convert(Matrix{T},a),
                                        convert(Matrix{T},b), convert(Vector{T},c) )
end

# First same as last, c.f. H&W p.167
isFSAL(btab::TableauRKExplicit) = btab.a[end,:]==btab.b[1,:] && btab.c[end]==1 # the latter is not needed really

## Tableaus for explicit methods
# Fixed step:
const bt_feuler = TableauRKExplicit(:feuler,(1,), Float64,
                             zeros(Int,1,1),
                            [1]',
                            [0]
                             )
const bt_midpoint = TableauRKExplicit(:midpoint,(2,), Float64,
                               [0  0
                                .5  0],
                              [0, 1]',
                              [0, .5]
                              )
const bt_heun = TableauRKExplicit(:heun,(2,), Float64,
                           [0  0
                            1  0],
                          [1//2, 1//2]',
                          [0, 1])

const bt_rk4 = TableauRKExplicit(:rk4,(4,),Float64,
                          [0 0 0 0
                           1//2 0 0 0
                           0 1//2 0 0
                           0 0 1 0],
                         [1//6, 1//3, 1//3, 1//6]',
                         [0, 1//2, 1//2, 1])
                         
# Adaptive step:

# Fehlberg https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
const bt_rk45 = TableauRKExplicit(:fehlberg,(4,5),Float64,
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
const bt_dopri5 = TableauRKExplicit(:dopri, (5,4), Float64,
                     [0   0 0 0 0 0 0
                      1//5 0 0 0 0 0 0
                      3//40     9//40 0 0 0 0 0 
                      44//45     -56//15     32//9 0 0 0 0
                      19372//6561     -25360//2187     64448//6561     -212//729 0 0 0
                      9017//3168     -355//33     46732//5247     49//176     -5103//18656 0 0
                      35//384         0     500//1113       125//192      -2187//6784         11//84      0],
                     [35//384         0     500//1113       125//192      -2187//6784         11//84      0
                      5179//57600     0     7571//16695     393//640     -92097//339200     187//2100     1//40],
                     [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
                     )

# Fehlberg 7(8) coefficients
# Values from pag. 65, Fehlberg, Erwin. "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control".
# National Aeronautics and Space Administration.
const bt_feh78 = TableauRKExplicit(:feh78, (7,8), Float64,
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
export oderk_fixed, ode_adapt

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

# Fixed step Runge-Kutta method
# TODO: iterator method
function oderk_fixed{N,S,T}(fn, y0, tspan::AbstractVector,
                            btab::TableauRKExplicit{N,S,T})
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
            calc_next_k!(ks, ytmp, ytmp, s, fn, tspan[i], dt, dof, btab)
            for d=1:dof
                ys[d,i+1] += dt * btab.b[s]*ks[d,s]
            end
        end
    end
    return tspan, transformys(ys)
end
# calculates k[s]
function calc_next_k!{N,S,T}(ks::Matrix{T}, ytmp::Vector, y, s, fn, t, dt, dof, btab::TableauRKExplicit{N,S,T})
    # Calculates the next ks and puts it into ks[:,s]
    # - ks and ytmp are modified inside this function.
    ytmp[:] = y
    for ss=1:s-1, d=1:dof
        ytmp[d] += dt * ks[d,ss] * btab.a[s,ss]
    end
    ks[:,s] = fn(t + btab.c[s]*dt, ytmp) # ::Vector{T}
    nothing
end

ode4_v2(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_rk4)
ode1_euler(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_feuler)
ode2_midpoint(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_midpoint)
ode2_heun(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_heun)

# Does one embedded R-K step updating ytrial, yerr and ks.
function rk_embedded_step!{N,S}(ytrial, yerr, ks, ytmp, y, fn, t, dt, dof, btab::TableauRKExplicit{N,S})
    # Assumes that ks[:,1] is already calculated!
    #
    # Modifies ytrial, yerr, ks, and ytmp
    ytrial[:] = 0
    yerr[:] = 0
    for d=1:dof
        ytrial[d] += btab.b[1,1]*ks[d,1]
        yerr[d]   += btab.b[2,1]*ks[d,1]
    end
    for s=2:S
        calc_next_k!(ks, ytmp, y, s, fn, t, dt, dof, btab)
        for d=1:dof
            ytrial[d] += btab.b[1,s]*ks[d,s]
            yerr[d]   += btab.b[2,s]*ks[d,s]
        end
        
    end
    for d=1:dof
        yerr[d]   = dt * (ytrial[d]-yerr[d])
        ytrial[d] = y[d] + dt * ytrial[d]
    end
end


# Adaptive ODE time stepper
function ode_adapt{N,S,T}(fn, y0, tspan, btab::Tableau{N,S,T};
                     reltol = 1.0e-5, abstol = 1.0e-8,
                     norm=Base.norm,
                     minstep=abs(tspan[end] - tspan[1])/1e9,
                     maxstep=abs(tspan[end] - tspan[1])/2.5,
                     initstep=0.)

    !isadaptive(btab) && error("Can only use this solver with an adpative RK Butcher table")

    # parameters
    timeout_const = 5 # after step reduction do not increase step for
                      # timeout_const steps
    
    order = minimum(btab.order)

    ## Initialization
    dof = length(y0)
    tspan = convert(Vector{T}, tspan)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]
    timeout = 0 # for step-control

    # work arrays:
    y = zeros(T, dof)      # y at time t
    y[:] = y0
    ytrial = zeros(T, dof) # trial solution at time t+dt
    yerr   = zeros(T, dof) # error of trial solution
    ks     = zeros(T, dof, S) 
    ytmp   = zeros(T, dof)

    # If tspan is a more than a length two vector: return solution at
    # those points only
    tstepsgiven = length(tspan)>2
    if tstepsgiven
        nsteps = length(tspan)
        ys = zeros(T, dof, nsteps)
        ys[:,1] = y
        iter = 1 # the index into tspan
    else
        ys = copy(y)
        tspan = [tstart]
    end
    # Time:
    dt, tdir, ks[:,1] = hinit(fn, y, tstart, tend, order, reltol, abstol)
    if initstep!=0
        dt = sign(initstep)==tdir ? initstep : error("initstep has wrong sign.")
    end
    # Diagnostics
    dts = Float64[]
    errs = Float64[]
    steps = [0,0]  # [accepted, rejected]
    laststep = false
    # Integration loop
    ii = 1
    while true
        # do one step (assumes ks[:,1]==f0)
        rk_embedded_step!(ytrial, yerr, ks, ytmp, y, fn, t, dt, dof, btab)

        # Check error and find a new step size:
        err, newdt, timeout = stepsize_hw92(dt, tdir, y, ytrial, yerr, abstol,
                                            reltol, order, timeout, dof, maxstep)

        if err<=1.0 # accept step
            # diagnostics
            steps[1] +=1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            if tstepsgiven
                # interpolate onto given output points
                f0 = slice(ks, 1:dof, 1)
                f1 = isFSAL(btab) ? slice(ks, 1:dof, S) : fn(t+dt, ytrial)
                while iter<nsteps && (tdir*tspan[iter+1]<=tdir*(t+dt) || laststep) # output at all new times which are ≤ t+dt
                    iter += 1
                    hermite_interp!(ys, iter, tspan[iter], t, dt, y, ytrial, f0, f1) # TODO: 3rd order only!
                end
                ks[:,1] = f1 # load ks[:,1] for next step
            else
                # output every step taken
                append!(ys, ytrial)
                push!(tspan, t+dt)
                # load ks[:,1] for next step
                ks[:,1] = isFSAL(btab) ? ks[:,end] : fn(t+dt, ytrial)
            end
            # Break if this was the last step:
            laststep && break

            # Swap bindings of y and ytrial, avoids one copy
            y, ytrial = ytrial, y

            # Update t to the time at the end of current step:
            t += dt
            dt = newdt

            # Hit end point exactly if next step within 1% of end:
            if tdir*(t+dt*1.01) >= tdir*tend
                dt = tend-t
                laststep = true # next step is the last, if it succeeds
            end
        elseif abs(newdt)<minstep  # minimum step size reached, break
            println("Warning: dt < minstep.  Stopping.")
            break
        else # redo step with smaller dt
            laststep = false
            steps[2] +=1 
            dt = newdt
            timeout = timeout_const
        end
    end
    if !tstepsgiven
        ys = reshape(ys, dof, length(tspan))
    end
#    @show steps
    return tspan, transformys(ys)
end
ode45_v2(fn, y0, tspan; kwargs...) = ode_adapt(fn, y0, tspan, bt_rk45; kwargs...)
ode54_v2(fn, y0, tspan; kwargs...) = ode_adapt(fn, y0, tspan, bt_dopri5; kwargs...)
ode78_v2(fn, y0, tspan; kwargs...) = ode_adapt(fn, y0, tspan, bt_feh78; kwargs...)

# Helper functions:
function stepsize_hw92(dt, tdir, x0, xtrial, xerr, abstol, reltol, order,
                       timeout, dof, maxstep)
    # Estimates new best step size following
    # Hairer & Wanner 1992, p167 (with some modifications)
    #
    # If timeout>0 no step size increase is allowed.

    # TODO:
    #  - parameters to lift out of this function
    #  - make type-agnostic
    timout_after_nan = 5
    fac = [0.8, 0.9, 0.25^(1/(order+1)), 0.38^(1/(order+1))][1] 
    facmax = 5.0 # maximal step size increase. 1.5-5
    facmin = 1./facmax  # maximal step size decrease. ?

    # if NaN present make step size smaller by maximum
    any(isnan(xtrial)) && return 10., dt*facmin, timout_after_nan

    # Eq 4.10:
    tol = abstol + max(abs(x0), abs(xtrial)).*reltol # 4.10

    # Eq. 4.11:
    # err = norm(1/dof^2*xerr./tol, 2)     
    err = norm(xerr./tol, 2)     # sans 1/dof
    #err = norm(xerr./tol, Inf)

    # Eq 4.13:
#    newdt = tdir*dt * min(facmax, max(facmin, fac*(1/err)^(1/(order+1))))
    newdt = min(maxstep, tdir*dt*max(facmin, fac*(1/err)^(1/(order+1)))) # modified
    if timeout>0
        newdt = min(newdt, dt)
        timeout -= 1
    end
    return err, tdir*min(newdt, maxstep), timeout
end

# For dense output see Hairer & Wanner p.190 using Hermite interpolation
function hermite_interp!(ys, iter, tquery,t,dt,y0,y1,f0,f1)
    # f_0 = f(x_0 , y_0) , f_1 = f(x_0 + h, y_1 )
    # this is O(3). TODO for higher order.
    #
    # Updates ys[:,iter] in-place
    theta = (tquery-t)/dt
    for i=1:length(y0)
        ys[i,iter] = ((1-theta)*y0[i] + theta*y1[i] + theta*(theta-1) *
                      ((1-2*theta)*(y1[i]-y0[i]) + (theta-1)*dt*f0[i] + theta*dt*f1[i]) )
    end
    nothing
end
