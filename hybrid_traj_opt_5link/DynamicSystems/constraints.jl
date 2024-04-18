include("biped5link.jl") # Defines hybrid system functions

function discrete_stance1_dynamics(model::NamedTuple, x::Vector, u::Vector, dt::Real)::Vector
    # dynamics when foot 1 is in contact with the ground 

    # x, y, t1, t2, t3, t4, t5 = q
    # xd, yd, t1d, t2d, t3d, t4d, t5d = q̇
    
    q, q̇ = x[1:7], x[8:14]
    r = biped5link_kinematics(q, model)
    fpos = r[1,:]
    q̇ₖ₊₁ = kkt_newton_step(q, q̇, u, model, fpos, left_foot_constraint, dt)
    qₖ₊₁ = q + q̇ₖ₊₁*dt
    xₖ₊₁ = [qₖ₊₁; q̇ₖ₊₁]

    return xₖ₊₁
end

function discrete_stance2_dynamics(model::NamedTuple, x::Vector, u::Vector, dt::Real)::Vector
    # dynamics when foot 1 is in contact with the ground 

    # x, y, t1, t2, t3, t4, t5 = q
    # xd, yd, t1d, t2d, t3d, t4d, t5d = q̇
    
    q, q̇ = x[1:7], x[8:14]
    r = biped5link_kinematics(q, model)
    fpos = r[5,:]
    q̇ₖ₊₁ = kkt_newton_step(q, q̇, u, model, fpos, _foot_constraint, dt)
    qₖ₊₁ = q + q̇ₖ₊₁*dt
    xₖ₊₁ = [qₖ₊₁; q̇ₖ₊₁]
    
    return xₖ₊₁
end