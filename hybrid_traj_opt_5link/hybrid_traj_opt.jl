import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
import MathOptInterface as MOI
import Ipopt 
import ForwardDiff as FD
import Convex as cvx 
import ECOS
import Plots

using LinearAlgebra
using Plots
using Random
using JLD2
using Test
using StaticArrays
using Printf

plotlyjs()

## 

include(joinpath(@__DIR__, "utils","fmincon.jl"))
include(joinpath(@__DIR__, "utils","walker.jl"))
include("DynamicSystems/biped5link.jl")
include("DynamicSystems/constraints.jl")
include("DynamicSystems/animate_walker.jl")

##

# feel free to solve this problem however you like, below is a template for a 
# good way to start. 

function create_idx(nx,nu,N)
    # create idx for indexing convenience
    # x_i = Z[idx.x[i]]
    # u_i = Z[idx.u[i]]
    # and stacked dynamics constraints of size nx are 
    # c[idx.c[i]] = <dynamics constraint at time step i>
    #
    # feel free to use/not use this 
    
    # our Z vector is [x0, u0, x1, u1, …, xN]
    nz = (N-1) * nu + N * nx # length of Z 
    x = [(i - 1) * (nx + nu) .+ (1 : nx) for i = 1:N]
    u = [(i - 1) * (nx + nu) .+ ((nx + 1):(nx + nu)) for i = 1:(N - 1)]
    
    # constraint indexing for the (N-1) dynamics constraints when stacked up
    c = [(i - 1) * (nx) .+ (1 : nx) for i = 1:(N - 1)]
    nc = (N - 1) * nx # (N-1)*nx 
    
    return (nx=nx, nu=nu, N=N, nz=nz, nc=nc, x=x, u=u, c=c)
end

function walker_cost(params::NamedTuple, Z::Vector)::Real
    # cost function 
    idx, N, xg = params.idx, params.N, params.xg
    Q, R, Qf = params.Q, params.R, params.Qf
    Xref,Uref = params.Xref, params.Uref
    
    # TODO: input walker LQR cost 
    
    J = 0 
    for i = 1:(N-1)
        xi = Z[idx.x[i]]
        ui = Z[idx.u[i]]
       
        J += 0.5*(xi - Xref[i])'*Q*(xi - Xref[i])
        J += 0.5*(ui - Uref[i])'*R*(ui - Uref[i]) 
    end
    
    xn = Z[idx.x[N]]
    J += 0.5*(xn - xg)'*Qf*(xn - xg)
        
    return J 
end

function walker_dynamics_constraints(params::NamedTuple, Z::Vector)::Vector
    idx, N, dt = params.idx, params.N, params.dt
    M1, M2 = params.M1, params.M2 
    J1, J2 = params.J1, params.J2 
    model = params.model 
        
    # create c in a ForwardDiff friendly way (check HW0)
    c = zeros(eltype(Z), idx.nc)
    
    # TODO: input walker dynamics constraints (constraints 3-6 in the opti problem)
    for k = 1:(N-1)
        xk   = Z[idx.x[k]]
        uk   = Z[idx.u[k]]
        xkp1 = Z[idx.x[k+1]]

        if (k in M1) # (not in J1) is implied 
            c[idx.c[k]] = discrete_unconstrained_dynamics(model, xk, uk, dt) - xkp1 
        elseif (k in M2) # (not in J1) is implied 
            c[idx.c[k]] = discrete_unconstrained_dynamics(model, xk, uk, dt) - xkp1 
        end

    end

    return c 
end

function walker_stance_constraint(params::NamedTuple, Z::Vector)::Vector
    idx, N, dt = params.idx, params.N, params.dt
    M1, M2 = params.M1, params.M2 
    J1, J2 = params.J1, params.J2 
    
    model = params.model 

    # create c in a ForwardDiff friendly way (check HW0)
    c = zeros(eltype(Z), 2*N)
    
    # TODO: add walker stance constraints (constraints 7-8 in the opti problem)
    for i = 1:2:2*(N-1) 
        k = Int((i+1)/2)
        x = Z[idx.x[k]]
        xkp1 = Z[idx.x[k+1]]
        q = x[1:7]
        qkp1 = xkp1[1:7]
        r = biped5link_kinematics(q, model)
        rkp1 = biped5link_kinematics(qkp1, model)
        if (k in J1) || (k in J2)
            c[i] = r[1,2] - height_stairs(r[1,1])
            c[i+1] = r[5,2] - height_stairs(r[5,1])
        elseif (k in M1)
            foot_diff = r[1,:] - rkp1[1,:]
            c[i] = foot_diff[1]
            c[i+1] = r[1,2] - height_stairs(r[1,1])
        elseif (k in M2)
            foot_diff = r[5,:] - rkp1[5,:]
            c[i] = foot_diff[1]
            c[i+1] = r[5,2] - height_stairs(r[5,1])
        end
    end

    return c
end
    
function walker_equality_constraint(params::NamedTuple, Z::Vector)::Vector
    N, idx, xic, xg = params.N, params.idx, params.xic, params.xg 
    
    # TODO: stack up all of our equality constraints 
    
    # should be length 2*nx + (N-1)*nx + N 
    # inital condition constraint (nx)       (constraint 1)
    # terminal constraint         (nx)       (constraint 2)
    # dynamics constraints        (N-1)*nx   (constraint 3-6)
    # stance constraint           N          (constraint 7-8)
    [   
      Z[idx.x[1]] - xic;
      Z[idx.x[N]] - xg;
      walker_dynamics_constraints(params, Z);
      walker_stance_constraint(params, Z)
    ]
end

function walker_inequality_constraint(params::NamedTuple, Z::Vector)::Vector
    idx, N, dt = params.idx, params.N, params.dt
    M1, M2 = params.M1, params.M2 
    
    cons = 6

    # create c in a ForwardDiff friendly way (check HW0)
    c = zeros(eltype(Z), cons*N)
    
    # TODO: add the length constraints shown in constraints (9-10)
    # there are 2*N constraints here 
    
    iter = 1 
    for k = 1:cons:cons*N
        xk = Z[idx.x[Int((k+cons-1)/cons)]]
        r = biped5link_kinematics(xk[1:7], params.model)
        r1y = r[1,2] - height_stairs(r[1,1])
        r2y = r[2,2] - height_stairs(r[2,1])
        r4y = r[4,2] - height_stairs(r[4,1])
        r5y = r[5,2] - height_stairs(r[5,1])
        c[k:k+3] = vcat([r1y; r2y; r4y; r5y]...)

        px, py, θ1, θ2, θ3, θ4, θ5 = xk[1:7]

        c[k+4] = θ2 - θ1
        c[k+5] = θ4 - θ5

        # c[k+8] = θ1 - (θ2 - π)
        # c[k+9] = θ5 - (θ4 - π)

        # c[k+6] = θ3 - (π/2 - π/6) 
        # c[k+7] = (π/2 + π/6) - θ3
    end
    return c
end

##

# dynamics parameters 
model = (m1 = .5,  m2 = .5,  m3 = 1,  m4 = .5,  m5 = .5,  m6 = .01,
            l12 = 1, l23 = 1, l34 = 1, l45 = 1, l36 = 1, g = 9.81)

# problem size 
nx = 14 
nu = 4 
# tf = 4.4
tf = 6 
dt = 0.1
t_vec = 0:dt:tf 
N = length(t_vec)

# initial and goal states
# determine the fixed joint angles
theta_offset = 20 # determines, how spread out we want the legs

q1 = 280 * (π/180)
q2 = 300 * (π/180)
q3 = 90 * (π/180)
q4 = 250 * (π/180)
q5 = 220 * (π/180)

x0 = - model.l23 * cos(q2) - model.l12 * cos(q1)
y0 = - model.l23 * sin(q2) - model.l12 * sin(q1)

dx0 = dy0 = dq1 = dq2 = dq3 = dq4 = dq5 = 0

# #determine the fixed x and y positions

xic = [ x0;  y0;  q1;  q2;  q3;  q4;  q5; 
        dx0; dy0; dq1; dq2; dq3; dq4; dq5]
# dx = 5 # suppose our goal is to move like 5 meters forward
dx = 5 # suppose our goal is to move like 5 meters forward

# D = [x0 + dx, y0]
# D_norm = norm(D)
# ϕ = atan(D[2], D[1])

# A = acos( (D_norm^2 + model.l12^2 - model.l23^2) / (2 * D_norm * model.l12) )
# q1 = ϕ + π - A
# C = acos( (D_norm^2 + model.l23^2 - model.l12^2) / (2 * D_norm * model.l23) )
# q2 = C + π + ϕ

xg = [x0 + dx;  y0 + height_stairs(x0 + dx);  q1;  q2;  q3;  q4;  q5; 
            dx0; dy0; dq1; dq2; dq3; dq4; dq5]

# xg = [x0 + dx;  y0 + height_stairs(x0 + dx);  q5;  q4;  q3;  q2;  q1; 
#             dx0; dy0; dq1; dq2; dq3; dq4; dq5]

# xg = [x0 + dx;  y0;  q5;  q4;  q3;  q2;  q1; 
#           dx0; dy0; dq1; dq2; dq3; dq4; dq5]

# index sets 
# M1 = vcat([1:20, 41:60, 81:100]...)
# M2 = vcat([21:40, 61:80, 101:121]...)
# J1 = [20, 60, 100]
# J2 = [40, 80] 

M1 = vcat([1:10,  26:38, 53:61]...)
M2 = vcat([11:25, 39:52]...)
J1 = [10, 38]
J2 = [25, 52] 

# M1 = vcat([1:30]...)
# M2 = vcat([31:61]...)
# J1 = [30]
# J2 = [62] 

# reference trajectory 
Xref, Uref = reference_trajectory(model, xic, xg, dt, N, M1, tf)

# animate_walker(Xref, model)

# LQR cost function (tracking Xref, Uref)
# Q = diagm([1; 10; fill(1.0, 5); 1; 10; fill(1.0, 5)]);
# TODO: change this ↓ to maximize cg position along trajectory
Q = diagm(fill(1.0,14))
# Q[2,2] = 100
R = diagm(fill(1e-3,4))
Qf = 1*Q;

# create indexing utilities 
idx = create_idx(nx,nu,N)

# put everything useful in params 
params = (
    model = model, 
    nx = nx, nu = nu,
    tf = tf, dt = dt, 
    t_vec = t_vec,
    N = N, 
    M1 = M1, M2 = M2, J1 = J1, J2 = J2,
    xic = xic, xg = xg,
    idx = idx,
    Q = Q, R = R, Qf = Qf,
    Xref = Xref, 
    Uref = Uref
)

# TODO: primal bounds (constraint 11)
x_l = -Inf*ones(idx.nz)
x_u =  Inf*ones(idx.nz)
# for i = 1:N
#     x_l[idx.x[i][[2]]] .= 0 
# end

# TODO: inequality constraint bounds
cons = 6
c_l = 0*ones(cons*N)
# c_l = -Inf*ones(cons*N)
c_u = Inf*ones(cons*N)

# TODO: initial guess, initialize z0 with the reference Xref, Uref 
z0 = zeros(idx.nz)

for i = 1:N 
    z0[idx.x[i]] = 1*Xref[i]
end
for i = 1:N-1 
    z0[idx.u[i]] = 1*Uref[i]
end

# adding a little noise to the initial guess is a good idea 
z0 = z0 + (1e-3)*randn(idx.nz)

diff_type = :auto 

Z = fmincon(walker_cost,walker_equality_constraint,walker_inequality_constraint,
            x_l,x_u,c_l,c_u,z0,params, diff_type;
            tol = 1e-6, c_tol = 1e-6, max_iters = 10_000, verbose = true)

# pull the X and U solutions out of Z 
X = [Z[idx.x[i]] for i = 1:N]
U = [Z[idx.u[i]] for i = 1:(N-1)]

# ------------plotting--------------
Xm = hcat(X...)
Um = hcat(U...)

animate_walker(X, model)

display(plot(Xm[1,:],Xm[2,:], label = "body"))
display(plot(t_vec[1:end-1], Um',xlabel = "time (s)", ylabel = "U",
                label = ["u1" "u2" "u3" "u4"], title = "Controls"))