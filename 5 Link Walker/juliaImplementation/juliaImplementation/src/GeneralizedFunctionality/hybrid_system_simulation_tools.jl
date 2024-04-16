#These functions are weird things to take intuitive dynamics
#functions into the forms needed for the integrator
function guard_vector_function(guards)
    
    return function(out,u,t,integrator)
        for i = 1:size(guards,1)
            out[i] = guards[i](u,t)
        end
    end
end

function dynamics_applier_function(dynamics,dim)
    return function(dx,x,p,t)
        res = dynamics(x,p,t)
        for i = 1:dim
            dx[i] = res[i]
        end
    end
end

function guard_condition_affect!(integrator,idx)
    terminate!(integrator)
end

#Simulate a single step in a trajectory
function dynamic_step(x,mode,dynamics,resets,guards, t, dt, p)
    curr_mode = deepcopy(mode)
    curr_t = deepcopy(t)
    curr_state = deepcopy(x)

    impact_state = x*NaN
    impact_time = NaN
    while abs(curr_t - (t+dt)) > 1e-16
        tspan = (curr_t,t+dt)
        dyn = dynamics_applier_function(dynamics[curr_mode],size(curr_state,1))
        cb_condition = guard_vector_function(guards[curr_mode])
        guard_cb = VectorContinuousCallback(cb_condition,nothing,guard_condition_affect!,size(dynamics,1))
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(dyn,curr_state,tspan,p)
        sol = solve(prob,Tsit5(),callback=guard_cb,dt=1e-3,adaptive=false)
        curr_state = deepcopy(sol[end])
        curr_t = deepcopy(sol.t[end])
        if abs(curr_t - (t+dt)) > 1e-16
            guardresults = [f(curr_state,t) for f in guards[curr_mode]]
            new_mode = findfirst(abs.(guardresults) .< 1e-8)
            impact_state = deepcopy(curr_state)
            curr_state = resets[curr_mode][new_mode](curr_state)
            curr_mode = new_mode
            impact_time = curr_t
        end
    end
    return curr_state, curr_mode, impact_state,impact_time
    
end 

#Simulate an entire trajectory
function sim_system(x,mode,dynamics,resets,guards,W,V,C,t,dt)
    times = 0:dt:t
    states = [x for i = 1:size(times,1)]
    measurements = 0
    # measurements = [C*x for i = 1:size(times,1)-1]
    modes = [mode for i = 1:size(times,1)]
    impact_states = [x*NaN for i = 1:size(times,1)-1]
    num_impacts = [0 for i = 1:size(times,1)]
    noise = [x*0 for i = 1:size(times,1)-1]
    for i = 2:size(times,1)
        # noise[i-1] = sqrt(W)*randn(size(states[1])) #gaussian
        # noise[i-1] = 0
        states[i], modes[i],impact_states[i-1],_ = dynamic_step(states[i-1],modes[i-1],dynamics,resets,guards,times[i-1],dt,noise[i-1])
        # measurements[i-1] = C*states[i] + sqrt(V)*randn(size(measurements[1])) #gaussian
        if !isnan(impact_states[i-1][1])
            num_impacts[i] = num_impacts[i-1] + 1
        else
            num_impacts[i] = num_impacts[i-1]
        end
    end
    return states, measurements, times, modes,impact_states,num_impacts,noise
end