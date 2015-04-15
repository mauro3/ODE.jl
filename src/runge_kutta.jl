# Butcher Tableaus
# see Hairer & Wanner 1992, p. 134, 166

abstract Tableau{T<:FloatingPoint}
# assumes fields
# order::(Int...) # order of the method(s)
# a::Matrix{T}  # Stage x Stage-1 matrix
# b::Matrix{T}  # 1 or 2 x Stage matrix (fixed step/ adaptive)
# c::Vector{T}  # Stage
#
# For a tableau:
#  c1  | a_11   ....   a_1(s-1)             
#  .   | a_21 .          .           
#  .   | a_31     .      .           
#  .   | ....         .  .           
#  c_s | a_s1  ....... a_s(s-1)     
# -----+-----------------------------
#      | b_1           b_(s-1)   b_s 
#      | b'_1          b'_(s-1)  b'_s 

Base.eltype{T}(b::Tableau{T}) = T

# #@doc "Test whether it's an explicit method"->
isexplicit(b::Tableau) = istril(b.a)
isadaptive(b::Tableau) = size(b.b, 1)==2
order(b::Tableau) = b.order

# The advantage of each having its own type, makes it possible to have
# specialized methods for a particular tablau
immutable TableauExplicit{Name, S, T} <: Tableau{T} # S is the number of stages (TODO needed as type parameter?)
    order::(Int...) # the order of the methods
    a::Matrix{T} # make this into a LowerTriangular matrix in 0.4
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function RKTableauExplicit(order,a,b,c)
        @assert c[1]==0
        @assert istril(a)
        @assert S==length(c)==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        new(order,a,b,c)
    end
end
function RKTableauExplicit{T}(name::Symbol, order::(Int...),
                   a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    RKTableauExplicit{name,length(c),T}(order, a, b, c)
end

update order cs,as,bs to a,b,c

# make a Runge-Kutta method for a given Butcher tableau.  Follows
# Hairer & Wanner 1992 p.165-169
export rk_runner, rk_runner_adaptive


function calc_ks{N,S,T}(fn, t0, y0::Vector, dt, btab::RKTableauExplicit{N,S,T})
    # calculate all ks[stage, dof]
    dof = length(y0)
    ks = zeros(T, S, dof)
    for i=1:S
        a = zeros(T,dof)
        for j=1:i-1
            for d =1:dof
                a[d] += btab.as[i,j]*ks[j,d]
            end
        end
        ks[i,:] = fn(t0 + dt*btab.cs[i], y0 + dt*a)
    end
    return ks
end

function rkstep_naive{N,S,T}(fn, t0, y0::Number, dt, btab::RKTableauExplicit{N,S,T})
    # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
    #
    # Only for scalar problems

    ks = calc_ks(fn, t0, [y0], dt, btab)
    y = zero(T)
    for i=1:S
        y += btab.bs[1,i]*ks[i]
    end
    return y0 + dt*y
end

function rkstep_naive{N,S,T}(fn, t0, y0::Vector, dt, btab::RKTableauExplicit{N,S,T})
    # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
    #
    # For vector problems

    ks = calc_ks(fn, t0, y0, dt, btab)

    dof = length(y0)
    y = zeros(T,dof)
    for d=1:dof
        for i=1:S
            y[d] += btab.bs[1,i].*ks[i,d]
        end
    end
    return y0 + dt*y
end

function rkstep_embedded_naive{N,S,T}(fn, t0, y0::Vector, dt, btab::RKTableauExplicit{N,S,T})
    # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
    #
    # For vector problems. Returns y and an error estimate

    ks = calc_ks(fn, t0, y0, dt, btab)

    dof = length(y0)
    y = zeros(T,dof)
    yerr = zeros(T,dof)
    for d=1:dof
        for i=1:S
            y[d]    += btab.bs[1,i].*ks[i,d]
            yerr[d] += btab.bs[2,i].*ks[i,d]
        end
    end
    return y0 .+ dt*y, dt*(y-yerr), squeeze(ks[1,:],1) # this is f0
end


function rk_runner(fn, ts, y0::Number, btab::RKTableauExplicit)
    # scalar y0
    ys = Array(typeof(y0), length(ts))
    ys[1] = y0
    for i=1:length(ts)-1
        dt = ts[i+1]-ts[i]
        ys[i+1] = rkstep_naive(fn, ts[i], ys[i], dt, btab)
    end
    return ys
end

function rk_runner(fn, ts, y0::Vector, btab::RKTableauExplicit)
    # vector y0
    ys = Array(eltype(y0), length(y0), length(ts))
    ys[:,1] = y0'
    for i=1:length(ts)-1
        dt = ts[i+1]-ts[i]
        ys[:,i+1] = rkstep_naive(fn, ts[i], ys[:,i], dt, btab)
    end
    return ys
end


function rk_runner_adaptive(fn, ts, y0::Vector, btab::RKTableauExplicit;
                            reltol=1e-6, abstol=1e-6, dt0=0.0, mindt=1e-5)
    # uses Hairer et al 1992 p.167
    if !isadaptive(btab)
        error("can only use this solver with an adpative Butcher table")
    end

    #########
    # helper functions:

    ###### end helper functions

    const large = 1.0e5

    ts, tstart, tend, tsgiven, t, ys, yold, dof, facmax = init(y0, ts, dt0)
    dt = get_intial_step(fn, y0, dt0)
    dts = Float64[]
    xerrs = Float64[]
    iter = 1
    steps = [0,0]  # [accepted, rejected]
    while t<tend
        ytrial, yerr, f0 =  rkstep_embedded_naive(fn, t, yold, dt, btab)
        newdt = stepsize_hw(dt, yold, ytrial, yerr,
                            abstol, reltol, order(btab), facmax, dof)

        if newdt>=dt # accept step new dt is larger than old one
            steps[1] +=1 
            if tsgiven
                # interpolate onto given output points
                f1 = fn(t+dt, ytrial)
                while iter<length(ts) && ts[iter+1]<=t+dt  # output at all new times which are ≤ t+dt
                    iter += 1
                    ys[:,iter] = hermite_interp(ts[iter], t, dt, yold, ytrial, f0, f1)
                end
            end

            yold = ytrial
            t += dt
            dt = newdt
            if (t+dt) > (tend + dt*0.01)
                # hit end point exactly
                dt = tend-t
            end
            push!(dts, dt)
            append!(xerrs, yerr)

            if !tsgiven
                append!(ys, ytrial)
                push!(ts, t)
            end
            facmax = large
        elseif dt<mindt  # minimum step size reached
            @show length(ys), t, dt
            error("dt < mindt")
            
        else # redo step with smaller dt
            steps[2] +=1 
            dt = newdt
            facmax = 1.0 # forbids dt increases in the next step
        end
    end
    if !tsgiven
        ys = reshape(ys, dof, length(ts))
    end
    xerrs = reshape(xerrs, dof, length(dts))
    return ys, ts, steps, dts, xerrs
end

# For dense output see p.188 using Hermite interpolation
function hermite_interp(tquery, t,dt,y0,y1,f0,f1)
    # this is O(3)
    theta = (tquery-t)/dt
    return (1-theta)*y0 + theta*y1 + theta*(theta-1) *
           ((1-2*theta)*(y1-y0) + (theta-1)*dt*f0 + theta*dt*f1)
end
