# Rosenbrock-Wanner methods
###########################
#
# Main references: Wanner & Hairrere 1996, PETSc http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSROSW.html
@show "NEW"
export ode_rosw, ode_rosw_fixed

# Rosenbrock-W methods are typically specified for autonomous DAE:
# Mẋ = f(x)
#
# by the stage equations
# M kᵢ = hf(x₀ + + Σⱼ aᵢⱼkⱼ) + h J Σⱼ γᵢⱼkⱼ
#
# and step completion formula
# x₁ = x₀ + \Sigmaⱼ bⱼ kⱼ
#
# The method used here uses transformed equations as done in PETSc.


# Tableaus
##########

immutable TableauRosW{Name, S, T} <: Tableau{Name, S, T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}
    γ::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    function TableauRosW(order,a,γ,b)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert S==size(γ,1)==size(a,2)==size(γ,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        new(order,a,γ,b)
    end
end
function TableauRosW{T}(name::Symbol, order::(@compat(Tuple{Vararg{Int}})),
                   a::Matrix{T}, γ::Matrix{T}, b::Matrix{T})
    TableauRosW{name,size(b,2),T}(order, a, γ, b)
end
function TableauRosW(name::Symbol, order::(@compat(Tuple{Vararg{Int}})), T::Type,
                   a::Matrix, γ::Matrix, b::Matrix)
    TableauRosW{name,size(b,2),T}(order, convert(Matrix{T},a),
                                  convert(Matrix{T},γ), convert(Matrix{T},b) )
end
conv_field{T,N}(D,a::Array{T,N}) = convert(Array{D,N}, a)
function Base.convert{Tnew<:Real,Name,S,T}(::Type{Tnew}, tab::TableauRosW{Name,S,T})
    # Converts the tableau coefficients to the new type Tnew
    newflds = ()
    @compat for n in fieldnames(tab)
        fld = getfield(tab,n)
        if eltype(fld)==T
            newflds = tuple(newflds..., conv_field(Tnew, fld))
        else
            newflds = tuple(newflds..., fld)
        end
    end
    TableauRosW{Name,S,Tnew}(newflds...) # TODO: could this be done more generically in a type-stable way?
end



# Transformed Tableau, used only internally
immutable TableauRosW_T{Name, S, T} <: Tableau{Name, S, T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}  # this is TableauRosW.a transformed
    γinv::Matrix{T}
    b::Matrix{T} # this is TableauRosW.b transformed
    # derived quantities:
    γii::T
    c::Matrix{T} # = (tril(btab.γinv)-diagm(diag(btab.γinv)))
end
function tabletransform{Name,S,T}(rt::TableauRosW{Name,S,T})
    # the code only works if
    if !all(x->x==rt.γ[1],diag(rt.γ))
        error("This Rosenbrock implementation only works for tableaus with γ_ii==γ_jj for all i,j.")
    end
    γii = rt.γ[1,1]
    γinv = inv(rt.γ)
    ahat = rt.a * γinv
    bhat = similar(rt.b)
    bhat[1,:] = squeeze(rt.b[1,:]*γinv,1)
    bhat[2,:] = squeeze(rt.b[2,:]*γinv,1)
    c = (tril(γinv)-diagm(diag(γinv))) # negative of W&H definition
    return TableauRosW_T{Name,S,T}(rt.order, ahat, γinv, bhat, γii, c)
end


## tableau for ros34pw2
const bt_ros34pw2 = TableauRosW(:ros34pw2, (3,4), Float64,
                                [0  0  0  0
                                 8.7173304301691801e-01  0  0  0
                                 8.4457060015369423e-01  -1.1299064236484185e-01  0  0
                                 0  0  1.  0],
                                [4.3586652150845900e-01  0  0  0
                                 -8.7173304301691801e-01  4.3586652150845900e-01  0  0
                                 -9.0338057013044082e-01  5.4180672388095326e-02  4.3586652150845900e-01  0
                                 2.4212380706095346e-01  -1.2232505839045147e+00  5.4526025533510214e-01  4.3586652150845900e-01],
                                 [2.4212380706095346e-01  -1.2232505839045147e+00  1.5452602553351020e+00  4.3586652150845900e-01
                                  3.7810903145819369e-01  -9.6042292212423178e-02  5.0000000000000000e-01  2.1793326075422950e-01]
                                  )

###################
# Fixed step solver
###################
ode_rosw_fixed(fn, Jfn, x0, tspan) = oderosw_fixed(fn, Jfn, x0, tspan, bt_ros34pw2)
function oderosw_fixed{N,S}(fn, Jfn, x0::AbstractVector, tspan,
                            btab::TableauRosW{N,S})
    # TODO: refactor with oderk_fixed
    Et, Exf, Tx, btab = make_consistent_types(fn, x0, tspan, btab)
    btab = tabletransform(btab)
    dof = length(x0)

    xs = Array(Tx, length(tspan))
    allocate!(xs, x0, dof)
    xs[1] = deepcopy(x0)

    tspan = convert(Vector{Et}, tspan)
    # work arrays:
    ks = Array(Tx, S)
    # allocate!(ks, x0, dof) # no need to allocate as fn is not in-place
    xtmp = similar(x0, Exf, dof)
    for i=1:length(tspan)-1
        dt = tspan[i+1]-tspan[i]
        xs[i+1] = rosw_step!(fn, Jfn, xs[i], dt, btab, 2)
    end
    return tspan, xs
end


ode_rosw(fn, Jfn, x0, tspan;kwargs...) = oderosw_adapt(fn, Jfn, x0, tspan, bt_ros34pw2; kwargs...)
function oderosw_adapt{N,S}(fn, Jfn, x0::AbstractVector, tspan, btab::TableauRosW{N,S};
                            reltol = 1.0e-5, abstol = 1.0e-8,
                            norm=Base.norm,
                            minstep=abs(tspan[end] - tspan[1])/1e9,
                            maxstep=abs(tspan[end] - tspan[1])/2.5,
                            initstep=0.,
                            points=:all
                            )
    # TODO: refactor with oderk_adapt

    ## Figure types
    fn_expl = (t,x)->fn(x, x*0)
    Et, Exf, Tx, btab = make_consistent_types(fn, x0, tspan, btab)
    btab = tabletransform(btab)
    # parameters
    order = minimum(btab.order)
    timeout_const = 5 # after step reduction do not increase step for
                      # timeout_const steps

    ## Setup
    dof = length(x0)
    tspan = convert(Vector{Et}, tspan)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # work arrays:
    x      = similar(x0, Exf, dof) # x at time t (time at beginning of step)
    x[:]   = x0                    # fill with IC
    xtrial = similar(x0, Exf, dof) # trial solution at time t+dt
    xerr   = similar(x0, Exf, dof) # error of trial solution
    k = Array(Tx, S) # stage variables
    allocate!(k, x0, dof)
    ks = zeros(Exf,dof) # work vector for one k
    h_store = zeros(Exf,dof) # work vector 
    jac_store = zeros(Exf,dof,dof) # Jacobian storage
    u = zeros(Exf,dof)    # work vector 
    udot = zeros(Exf,dof) # work vector 

    # output xs
    nsteps_fixed = length(tspan) # these are always output
    xs = Array(Tx, nsteps_fixed)
    allocate!(xs, x0, dof)
    xs[1] = x0

    # Option points determines where solution is returned:
    if points==:all
        tspan_fixed = tspan
        tspan = Et[tstart]
        iter_fixed = 2 # index into tspan_fixed
        sizehint!(tspan, nsteps_fixed)
    elseif points!=:specified
        error("Unrecognized option points==$points")
    end
    # Time
    dt, tdir, k[1] = hinit(fn_expl, x, tstart, tend, order, reltol, abstol) # sets k[1]=f0
    if initstep!=0
        dt = sign(initstep)==tdir ? initstep : error("initstep has wrong sign.")
    end
    # Diagnostics
    dts = Et[]
    errs = Float64[]
    steps = [0,0]  # [accepted, rejected]

    ## Integration loop
    laststep = false
    timeout = 0 # for step-control
    iter = 2 # the index into tspan and xs
    while true
        rosw_step!(xtrial, fn, Jfn, x, dt, dof, btab,
                   k, h_store, jac_store, ks, u, udot, 1)
        # Completion again for embedded method, see line 927 of
        # http://www.mcs.anl.gov/petsc/petsc-current/src/ts/impls/rosw/rosw.c.html#TSROSW
        xerr[:] = 0.0
        for s=1:S
            for d=1:dof
                xerr[d] += (btab.b[2,s]-btab.b[1,s])*k[s][d]
            end
        end
        err, newdt, timeout = stepsize_hw92!(dt, tdir, x, xtrial, xerr, order, timeout,
                                            dof, abstol, reltol, maxstep, norm)
        if err<=1.0 # accept step
            # diagnostics
            steps[1] +=1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            f0 = k[1]
            f1 = fn_expl(t+dt, xtrial)
            if points==:specified
                # interpolate onto given output points
                while iter-1<nsteps_fixed && (tdir*tspan[iter]<tdir*(t+dt) || laststep) # output at all new times which are < t+dt
                    hermite_interp!(xs[iter], tspan[iter], t, dt, x, xtrial, f0, f1) # TODO: 3rd order only!
                    iter += 1
                end
            else
                # first interpolate onto given output points
                while iter_fixed-1<nsteps_fixed && tdir*t<tdir*tspan_fixed[iter_fixed]<tdir*(t+dt) # output at all new times which are < t+dt
                    xout = hermite_interp(tspan_fixed[iter_fixed], t, dt, x, xtrial, f0, f1)
                    index_or_push!(xs, iter, xout) # TODO: 3rd order only!
                    push!(tspan, tspan_fixed[iter_fixed])
                    iter_fixed += 1
                    iter += 1
                end
                # but also output every step taken
                index_or_push!(xs, iter, copy(xtrial))
                push!(tspan, t+dt)
                iter += 1
            end
            k[1] = f1 # load k[1]==f0 for next step

            # Break if this was the last step:
            laststep && break

            # Swap bindings of x and xtrial, avoids one copy
            x, xtrial = xtrial, x

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
    return tspan, xs

end
function rosw_step!{N,S}(xtrial, g, gprime, x, dt, dof, btab::TableauRosW_T{N,S},
                         k, h_store, jac_store, ks, u, udot, bt_ind=1)
    # This takes one step for a ode/dae system defined by
    # g(x,xdot)=0
    # gprime(x, xdot, α) = dg/dx + α dg/dxdot

    # first step
    # ks[:] = 0 assumes ks==0
    s = 1
    h!(h_store, ks, u, udot, g, x, btab, dt, k, s, dof)
    # It's sufficient to only update the Jacobian once per time step:
    # (Note that this uses u & udot calculated in h! above.)
    jac, tmp = hprime!(jac_store, x, gprime, u, udot, btab, dt)
    
    # calculate k
    for s=1:S-1
        # first step of Newton iteration with guess ks==0
        k[s][:] = ks - jac\h_store # TODO use A_ldiv_B!

        h!(h_store, ks, u, udot, g, x, btab, dt, k, s+1, dof)
    end
    # last step
    k[S][:] = ks - jac\h_store

    # completion:
    xtrial[:] = x
    for s=1:S
        for d=1:dof
            xtrial[d] += btab.b[bt_ind,s]*k[s][d]
        end
    end
    return nothing
end
function h!(res, ks, u, udot, g, x, btab, dt, k, s, dof)
    # h(ks)=0 to be solved for ks.
    #
    # Calculates h(ks) and h'(ks) (if s==1).
    # Modifies its first 3 arguments.
    #
    # ks -- guess for k[s] (usually==0) AND result g(u,udot)
    # u, udot -- work arrays for first and second argument to g
    # g -- obj function
    # x -- current state
    # btab
    # k -- stage vars (jed's yi)
    # s -- current stage to calculate
    # dof
    #
    # TODO: could re-use ks for res

    # stage independent and first stage:
    ss = 1
    for d=1:dof
        u[d] = (ks[d] + x[d]
                + btab.a[s,ss]*k[ss][d])
        udot[d] = (1/(dt*btab.γii).*ks[d] # usually ==0 as ks==0
                   + btab.c[s,ss]/dt*k[ss][d])
    end
    # later stages:
    for ss=2:s-1
        for d=1:dof
            u[d] += btab.a[s,ss]*k[ss][d]
            udot[d] += btab.c[s,ss]/dt*k[ss][d]
        end
    end
    res[:] = g(u, udot)
    return nothing
end
function hprime!(res, xi, gprime, u, udot, btab, dt)
    # The Jacobian of h!  Note that u and udot should be calculated by
    # running h! first.
    res[:,:] = gprime(u, udot, 1./(dt*btab.γii))
    tmp = copy(res)
    return lufact!(res), tmp
end
