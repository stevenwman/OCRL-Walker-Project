include("biped5link.jl") # Defines hybrid system functions

function discrete_stance1_dynamics(model::NamedTuple, x::Vector, u::Vector, dt::Real)::Vector
    # dynamics when foot 1 is in contact with the ground 

    # x, y, t1, t2, t3, t4, t5 = q
    # xd, yd, t1d, t2d, t3d, t4d, t5d = q̇
    
    q, q̇ = x[1:7], x[8:14]
    r = biped5link_kinematics(q, model)
    fpos = [r[1,:]; r[5,:]]
    newton_step = kkt_newton_step(q, q̇, u, model, fpos, left_foot_constraint, dt)
    q̇ₖ₊₁ = newton_step[1:7]
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
    fpos = [r[1,:]; r[5,:]]
    newton_step = kkt_newton_step(q, q̇, u, model, fpos, right_foot_constraint, dt)
    q̇ₖ₊₁ = newton_step[1:7]
    qₖ₊₁ = q + q̇ₖ₊₁*dt
    xₖ₊₁ = [qₖ₊₁; q̇ₖ₊₁]
    
    return xₖ₊₁
end


function reference_trajectory(model, xic, xg, dt, N, M1, dx, tf)
    # creates a reference Xref and Uref for the walker 
    # for the reference, we're going to assume the walker is a rigid
    # body, and then we linearly interpolate a trajectory that has it move
    # horizonally at a constant velocity
    
    #assume Uref is just zeros for now, TODO: we can adjust this later to make it a sinusoial input or something
    Uref = [[0; 0; 0; 0] for i = 1:(N-1)]
    
    Xref = [zeros(14) for i = 1:N]
    
    # set first and last timesteps
    Xref[1] = 1*xic 
    Xref[N] = 1*xg
    
    # interpolate the x position values in between
    x_start = xic[1]
    x_end = xg[1]
    xs = range(x_start, x_end, length = N)

    #determine the fixed joint angles
    # theta_offset = 20 # determines, how spread out we want the legs
    # q1 = (270+theta_offset)*pi/180
    # q2 = (270+theta_offset)*pi/180
    # q3 = pi/2 
    # q4 = (270-theta_offset)*pi/180
    # q5 = (270-theta_offset)*pi/180

    q1 = 280 * (π/180)
    q2 = 300 * (π/180)
    q3 = 90 * (π/180)
    q4 = 250 * (π/180)
    q5 = 220 * (π/180)

    x0 = - model.l23 * cos(q2) - model.l12 * cos(q1)
    y0 = - model.l23 * sin(q2) - model.l12 * sin(q1)


    # #determine the fixed y-height
    # total_link_length = 2 # assume each link is 1 for now, need to adjust later
    # ys = total_link_length*(cos(20*pi/180)) # y position of body

    #determine the fixed x-horizontal velocity
    horiz_v = x_end / tf
    
    #now construct the Xref vector of vectors
    for i = 2:(N-1) 
        if i in M1
            Xref[i] = [x0 + horiz_v*i*dt , y0, q1, q2, q3, q4, q5, horiz_v, 0, 0, 0 ,0, 0, 0]
        else
            Xref[i] = [x0 + horiz_v*i*dt , y0, q5, q4, q3, q2, q1, horiz_v, 0, 0, 0 ,0, 0, 0]
        end
    end
        
    return Xref, Uref
end