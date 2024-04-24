include("biped5link.jl") # Defines hybrid system functions

function discrete_unconstrained_dynamics(model::NamedTuple, x::Vector, u::Vector, dt::Real)::Vector
    q, q̇ = x[1:7], x[8:14]
    return unconstrained_dynamics(q, q̇, u, model, dt)
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
        new_x = x0 + horiz_v*i*dt
        q = [new_x, y0, q1, q2, q3, q4, q5]
        r = biped5link_kinematics(q, model)
        new_y = y0 + height_stairs(r[1,1])
        if i in M1
            Xref[i] = [new_x , new_y, q1, q2, q3, q4, q5, 0, 0, 0, 0 ,0, 0, 0]
        else
            Xref[i] = [new_x , new_y, q5, q4, q3, q2, q1, 0, 0, 0, 0 ,0, 0, 0]
        end
    end
        
    return Xref, Uref
end

function height_stairs(x_distance)
# function to increase stairs step size by 20cm, for every 20cm horizontal distance
    
    # step height increase occurs every 20 cm of horizontal distance traveled
    step_height_increase = 1  # 20 cm in m
    # calculate the number of step height increases based on the x_distance traveled
    num_increases = floor(x_distance / step_height_increase)
    # calculate the current height of the staircase

    if x_distance < 0.2
        height = 0
    elseif x_distance < 1.5
        height = num_increases * step_height_increase / 2
    else
        height = 0.5
    end

    # height = 0

    return height
end