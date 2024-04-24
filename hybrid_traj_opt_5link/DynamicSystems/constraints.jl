include("biped5link.jl") # Defines hybrid system functions

function discrete_unconstrained_dynamics(model::NamedTuple, x::Vector, u::Vector, dt::Real)::Vector
    q, q̇ = x[1:7], x[8:14]
    return unconstrained_dynamics(q, q̇, u, model, dt)
end

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


function reference_trajectory(model, xic, xg, dt, N, M1, tf)
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

    #determine the fixed joint angles

    q1 = 280 * (π/180)
    q2 = 300 * (π/180)
    q3 = 90 * (π/180)
    q4 = 250 * (π/180)
    q5 = 220 * (π/180)

    x0 = - model.l23 * cos(q2) - model.l12 * cos(q1)
    y0 = - model.l23 * sin(q2) - model.l12 * sin(q1)

    #determine the fixed x-horizontal velocity
    horiz_v = (x_end - x_start)/ tf
    #now construct the Xref vector of vectors
    for i = 2:(N-1) 
        # if i in M1
        #     Xref[i] = [x0 + horiz_v*i*dt , y0, q1, q2, q3, q4, q5, horiz_v, 0, 0, 0 ,0, 0, 0]
        # else
        #     Xref[i] = [x0 + horiz_v*i*dt , y0, q5, q4, q3, q2, q1, horiz_v, 0, 0, 0 ,0, 0, 0]
        # end
        if i < N/2
            Xref[i] = xic
        else
            Xref[i] = xg
        end
    end
        
    return Xref, Uref
end

function height_stairs(x_distance)
# function to increase stairs step size by 20cm, for every 20cm horizontal distance
    
    # step height increase occurs every 20 cm of horizontal distance traveled
    step_height_increase = 0.2  # 20 cm in m
    # calculate the number of step height increases based on the x_distance traveled
    num_increases = floor(x_distance / step_height_increase)
    # calculate the current height of the staircase
    height = num_increases * step_height_increase

    return height
end