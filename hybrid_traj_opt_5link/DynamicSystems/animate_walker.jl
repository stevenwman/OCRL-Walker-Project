include("biped5link.jl")  # assume this file is in the same folder as 'biped5link.jl'

function animate_walker(state_matrix, params, guard=nothing, saveAnimation=true)

    """
    animate_walker(states, guard=nothing, saveAnimation=False)

    # Arguments
    - `state_matrix`: a 14 x N (# timesteps) matrix of states 
        - assume every row looks like [q, q_dot] where
            x, y, t1, t2, t3, t4, t5 = q
            xd, yd, t1d, t2d, t3d, t4d, t5d = q̇
    - 'params': tuple that includes 6 mass weights, 5 link lengths, and gravity constant
    - `guard`: currently a float that represents a slop angle, should probably be a time varying function in the future
    - 'saveAnimation' : boolean that determines whether or not to save a copy of the animation

    # Returns
    Nothing, but produces an animation that can optionally be saved locally

    """

    # extract timesteps


    N = size(state_matrix)[1]

    # define guard
    if guard === nothing
        #guard is not provided, use default of 0 slope ground
        slope = 0.0
    else
        #guard is provided
        slope = guard
    end

    # start animation
    anim = Animation()

    lim = 3.5 # likely needs to be tuned

    timesteps = N

    # Loop to create frames
    for i in 1:timesteps

        #extract joint positions at i-th timestep
        q_i = state_matrix[i]

        joint_pos = biped5link_kinematics(q_i, params)
        r1, r2, r3, r4, r5, r6 = eachrow(joint_pos)

        # Plot the data

        # set limits for the plot
        l_bound = -lim
        r_bound = 2*lim

        # calc guard and plot
        xs = l_bound:0.01:r_bound
        ys = -slope * xs # y = -m*x, passing through (0, 0)
        plot(xs, ys, linestyle=:solid, label="Ground", color="black")

        # now plot leg1
        leg1 = vcat(r1', r2', r3')
        plot!(leg1[:, 1], leg1[:, 2], xlim=(l_bound,r_bound),
                                ylim=(l_bound,r_bound), 
                                label="leg1", 
                                color="blue", 
                                linestyle=:solid, 
                                marker=:circle, 
                                aspect_ratio=:equal)

        # now plot leg2
        leg2 = vcat(r5', r4', r3')
        plot!(leg2[:, 1], leg2[:, 2], xlim=(l_bound,r_bound),
                                ylim=(l_bound,r_bound), 
                                label="leg2", 
                                color="red", 
                                linestyle=:solid, 
                                marker=:circle, 
                                aspect_ratio=:equal)

        # and plot torso
        torso = vcat(r3', r6')
        plot!(torso[:, 1], torso[:, 2], xlim=(l_bound,r_bound),
                                ylim=(l_bound,r_bound), 
                                label="torso", 
                                color="green", 
                                linestyle=:solid, 
                                marker=:circle, 
                                aspect_ratio=:equal)

        title!("Timestep $i")
        
        # Add frame to animation
        frame(anim)
    end

    if saveAnimation
        # Save the animation as a gif
        gif(anim, "5_link_walker.gif", fps = 30)
    end
end