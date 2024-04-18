include("biped5link.jl") # Defines hybrid system functions

function discrete_stance1_dynamics(params::NamedTuple, x::Vector, u::Vector, dt::Real)::Vector
    # dynamics when foot 1 is in contact with the ground 

    # x, y, t1, t2, t3, t4, t5 = q
    # xd, yd, t1d, t2d, t3d, t4d, t5d = q̇
    
    q, q̇ = x[1:7], x[8:14]

    kkt_newton_step(q, q̇, u, params, fpos, left_foot_constraint, dt)
    
    # rb  = x[1:2]   # position of the body
    # rf1 = x[3:4]   # position of foot 1
    # rf2 = x[5:6]   # position of foot 2
    # v   = x[7:12]  # velocities
    
    
    # ℓ1x = (rb[1]-rf1[1])/norm(rb-rf1)
    # ℓ1y = (rb[2]-rf1[2])/norm(rb-rf1)
    # ℓ2x = (rb[1]-rf2[1])/norm(rb-rf2)
    # ℓ2y = (rb[2]-rf2[2])/norm(rb-rf2)
    
    # B = [ℓ1x  ℓ2x  ℓ1y-ℓ2y;
    #      ℓ1y  ℓ2y  ℓ2x-ℓ1x;
    #       0    0     0;
    #       0    0     0;
    #       0  -ℓ2x  ℓ2y;
    #       0  -ℓ2y -ℓ2x]
    
    # v̇ = [0; -g; 0; 0; 0; -g] + M\(B*u)
    
    # ẋ = [v; v̇]
    
    return ẋ
end