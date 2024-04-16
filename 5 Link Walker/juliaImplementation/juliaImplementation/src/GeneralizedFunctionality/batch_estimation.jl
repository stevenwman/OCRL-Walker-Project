function runMultipleKF(dynamics,resets,guards,W,V,C,init_state,Σ,init_mode,t,dt; num_iters=100, randseed=1)
    param = ESTParameters(C,V,W,dt,dynamics,resets,guards)
    errors = [0.0 for i=1:num_iters]
    simtimes = [0.0 for i=1:num_iters]
    avg_step_errors = [0.0 for i=1:round(1+t/dt)]
    avg_error = [0.0.*init_state for i=1:round(1+t/dt)]
    step_errors = [[0.0 for i=1:round(1+t/dt)] for j = 1:num_iters]
    Random.seed!(randseed)
    for i = 1:num_iters
        real_state = deepcopy(init_state)
        real_state = real_state+sqrt(Σ)*randn(size(init_state))
        states, measurements, times = sim_system(real_state,init_mode,dynamics,resets,guards,W,V,C,t,dt)

        res,time = @timed skf(param,init_state,Σ,init_mode,measurements;verbose=false)
        Xf = res[1]
        modes = res[2]
        simtimes[i] = time
        err, stepErrors, error = calcErrors(states[1:end],Xf[1:end])
        errors[i] = err
        step_errors[i] = stepErrors
        avg_step_errors = avg_step_errors .+ abs.(stepErrors)/num_iters
        avg_error = avg_error .+ [abs.(err)/num_iters for err in error]

    end
    avgErr = Statistics.mean(errors)
    avgTime = Statistics.mean(simtimes)
    return avgErr, avg_step_errors, avg_error, avgTime, step_errors
end

function calcErrors(ref,est)
    error = est .- ref
    errorMag = LinearAlgebra.norm.(error)
    avgError = Statistics.mean(errorMag)

    return avgError,errorMag,error
    
end
    