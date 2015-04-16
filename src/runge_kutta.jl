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
function TableauExplicit(name::Symbol, order::(Int...), T,
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
                     [35//384     0     500//1113     125//192     -2187//6784     11//84     0
                      5179//57600     0     7571//16695     393//640     -92097//339200     187//2100     1//40],
                     [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
                     )


# make a Runge-Kutta method for a given Butcher tableau.  Follows
# Hairer & Wanner 1992 p.134, p.165-169
export oderk_fixed, oderk_adapt

# to put ys into the vector of vector format:
transformys(T,ys) = size(ys,1)==1 ? squeeze(ys,1) : Vector{T}[ys[:,i] for i=1:size(ys,2)]

# Fixed step Runge-Kutta method
# function oderk_fixed(fn, y0::Number, tspan, btab::TableauExplicit)
#     # scalar y0
#     ys = Array(typeof(y0), length(tspan))
#     ys[1] = y0
#     for i=1:length(tspan)-1
#         dt = tspan[i+1]-tspan[i]
#         ys[i+1] = rkstep_naive(fn, tspan[i], ys[i], dt, btab)
#     end
#     return vcat(tspan), ys
# end
function oderk_fixed{N,S,T}(fn, y0, tspan::AbstractVector,
                            btab::TableauExplicit{N,S,T})
    dof = length(y0)
    tsteps = length(tspan)
    ys = Array(T, dof, tsteps)
    ys[:,1] = y0'
    # work arrays
    ks = zeros(T, dof, S)
    ytmp = zeros(T, dof)
    for i=1:length(tspan)-1
        dt = tspan[i+1]-tspan[i]
        ys[:,i+1] = ys[:,i]
        for s=1:S
            calc_next_k!(ks, ytmp, s, fn, tspan[i], dt, ys, i, dof, btab)
            for d=1:dof
                ys[d,i+1] += dt * btab.b[s]*ks[d,s]
            end
        end
    end
    return tspan, transformys(T,ys)
end
# calculates k[s]
function calc_next_k!{N,S,T}(ks::Matrix, ytmp::Vector, s, fn, t, dt, ys, i, dof, btab::TableauExplicit{N,S,T})
    ytmp[:] = ys[:,i]
    for ss=2:s
        for d=1:dof
            ytmp[d] += dt * ks[d,ss] * btab.a[s,ss]
        end
    end
    @show 1
    ks[:,s] = fn(t + btab.c[s]*dt, ytmp)
end

ode4_v2(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_rk4)

# # Helper functions:S
# function calc_ks!{N,S,T}(ks, fn, t0, y0::Vector, dt, btab::TableauExplicit{N,S,T})
#     # calculate all k[dof, stage]
#     dof = length(y0)
#     for s=1:S
#         a = zeros(T,dof)
#         for j=1:s-1
#             for d =1:dof
#                 a[d] += btab.a[s,j]*ks[j,d]
#             end
#         end
#         ks[s,:] = fn(t0 + dt*btab.c[s], y0 + dt*a)
#     end
#     return ks
# end

# function rkstep_naive{N,S,T}(fn, t0, y0::Number, dt, btab::TableauExplicit{N,S,T})
#     # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
#     #
#     # Only for scalar problems

#     ks = calc_ks(fn, t0, [y0], dt, btab)
#     y = zero(T)
#     for i=1:S
#         y += btab.b[1,i]*ks[i]
#     end
#     return y0 + dt*y
# end

# function rkstep_naive{N,S,T}(fn, t0, y0::Vector, dt, btab::TableauExplicit{N,S,T})
#     # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
#     #
#     # For vector problems

#     ks = calc_ks(fn, t0, y0, dt, btab)

#     dof = length(y0)
#     y = zeros(T,dof)
#     for d=1:dof
#         for i=1:S
#             y[d] += btab.b[1,i].*ks[i,d]
#         end
#     end
#     return y0 + dt*y
# end


# # Adaptive Runge-Kutta method
# function oderk_adapt(fn, y0::Vector, tspan, btab::TableauExplicit;
#                      reltol = 1.0e-5, abstol = 1.0e-8,
#                      norm=Base.norm,
#                      minstep=abs(tspan[end] - tspan[1])/1e9,
#                      maxstep=abs(tspan[end] - tspan[1])/2.5, # TODO
#                      initstep=0.)

#     !isadaptive(btab) && error("can only use this solver with an adpative Butcher table")

#     const large = 1.0e5

#     tspan, tstart, tend, tspangiven, t, ys, yold, dof, facmax = init(y0, tspan, initstep)
#     dt = get_intial_step(fn, y0, initstep)
#     dts = Float64[]
#     xerrs = Float64[]
#     iter = 1
#     steps = [0,0]  # [accepted, rejected]
#     while t<tend
#         ytrial, yerr, f0 =  rkstep_embedded_naive(fn, t, yold, dt, btab)
#         newdt = stepsize_hw(dt, yold, ytrial, yerr,
#                             abstol, reltol, order(btab), facmax, dof)

#         if newdt>=dt # accept step new dt is larger than old one
#             steps[1] +=1 
#             if tspangiven
#                 # interpolate onto given output points
#                 f1 = fn(t+dt, ytrial)
#                 while iter<length(tspan) && tspan[iter+1]<=t+dt  # output at all new times which are ≤ t+dt
#                     iter += 1
#                     ys[:,iter] = hermite_interp(tspan[iter], t, dt, yold, ytrial, f0, f1)
#                 end
#             end

#             yold = ytrial
#             t += dt
#             dt = newdt
#             if (t+dt) > (tend + dt*0.01)
#                 # hit end point exactly
#                 dt = tend-t
#             end
#             push!(dts, dt)
#             append!(xerrs, yerr)

#             if !tspangiven
#                 append!(ys, ytrial)
#                 push!(tspan, t)
#             end
#             facmax = large
#         elseif dt<minstep  # minimum step size reached
#             @show length(ys), t, dt
#             error("dt < minstep")
            
#         else # redo step with smaller dt
#             steps[2] +=1 
#             dt = newdt
#             facmax = 1.0 # forbids dt increases in the next step
#         end
#     end
#     if !tspangiven
#         ys = reshape(ys, dof, length(tspan))
#     end
#     xerrs = reshape(xerrs, dof, length(dts))
#     return ys, tspan, steps, dts, xerrs
# end

# # Helper functions:
# function calc_ks{N,S,T}(fn, t0, y0::Vector, dt, btab::TableauExplicit{N,S,T})
#     # calculate all k[stage, dof]
#     dof = length(y0)
#     ks = zeros(T, S, dof)
#     for i=1:S
#         a = zeros(T,dof)
#         for j=1:i-1
#             for d =1:dof
#                 a[d] += btab.a[i,j]*ks[j,d]
#             end
#         end
#         ks[i,:] = fn(t0 + dt*btab.c[i], y0 + dt*a)
#     end
#     return ks
# end

# function rkstep_naive{N,S,T}(fn, t0, y0::Number, dt, btab::TableauExplicit{N,S,T})
#     # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
#     #
#     # Only for scalar problems

#     ks = calc_ks(fn, t0, [y0], dt, btab)
#     y = zero(T)
#     for i=1:S
#         y += btab.b[1,i]*ks[i]
#     end
#     return y0 + dt*y
# end

# function rkstep_naive{N,S,T}(fn, t0, y0::Vector, dt, btab::TableauExplicit{N,S,T})
#     # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
#     #
#     # For vector problems

#     ks = calc_ks(fn, t0, y0, dt, btab)

#     dof = length(y0)
#     y = zeros(T,dof)
#     for d=1:dof
#         for i=1:S
#             y[d] += btab.b[1,i].*ks[i,d]
#         end
#     end
#     return y0 + dt*y
# end

# function rkstep_embedded_naive{N,S,T}(fn, t0, y0::Vector, dt, btab::TableauExplicit{N,S,T})
#     # Does an S-stage explicit Runge-Kutta step for RHS f(t,y)
#     #
#     # For vector problems. Returns y and an error estimate

#     ks = calc_ks(fn, t0, y0, dt, btab)

#     dof = length(y0)
#     y = zeros(T,dof)
#     yerr = zeros(T,dof)
#     for d=1:dof
#         for i=1:S
#             y[d]    += btab.b[1,i].*ks[i,d]
#             yerr[d] += btab.b[2,i].*ks[i,d]
#         end
#     end
#     return y0 .+ dt*y, dt*(y-yerr), squeeze(ks[1,:],1) # this is f0
# end

# # For dense output see Hairer & Wanner p.188 using Hermite interpolation
# function hermite_interp(tquery, t,dt,y0,y1,f0,f1)
#     # this is O(3)
#     theta = (tquery-t)/dt
#     return (1-theta)*y0 + theta*y1 + theta*(theta-1) *
#            ((1-2*theta)*(y1-y0) + (theta-1)*dt*f0 + theta*dt*f1)
# end
