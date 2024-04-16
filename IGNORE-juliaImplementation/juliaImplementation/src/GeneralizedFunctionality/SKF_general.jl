function skf(param,x0,Σ0,mode0,y;verbose=false)
	N = size(y,1)+1
    X = [zeros(size(x0,1)) for i = 1:N]
    modes = [mode0 for i = 1:N]

    X[1] = x0
    Σ = Σ0

    for i = 1:N-1
    	t = param.dt*(i-1)
    	X[i+1],Σ,modes[i+1],impact = skfStep(param,X[i],Σ,modes[i],y[i],t)
		
    	if verbose
    		display(X[i+1])
    	end
    end
    return X,modes
end

function skfPredict(param,x,Σ,mode,t)

	#simulate the mean forward, storing the impact time and the mode it changes to
	xn, newmode, impactstate, impacttime = dynamic_step(x,mode,param.dynamics,param.resets,param.guards,t,param.dt,x*0)

	#update the cov based on the dynamic linearization, the process noise
	#and the saltation matrix if there is an impact
	impact = false
	#case 1: no impact
	if isnan(impacttime)
		A = calc_A_disc(param.dynamics[mode],x,x*0,t,param.dt)
		Σn = A*Σ*A' + param.W*param.dt^2
	#case 2: impact
	else
		impact = true
		A1 = calc_A_disc(param.dynamics[mode],x,x*0,t,impacttime-t)
		x_post = param.resets[mode][newmode](impactstate)
		A2 = calc_A_disc(param.dynamics[newmode],x_post,x_post*0,impacttime,t+param.dt-impacttime)
		Ξ = calc_salt(impactstate,impactstate*0,t,mode,newmode,param.dynamics,param.resets,param.guards)
		Σn = A2*Ξ*A1*Σ*A1'*Ξ'*A2'
	end

	return xn,Σn,newmode,impact
end

function skfUpdate(param,x,Σ,mode,meas,t,predimpact)
	impact = false
	# Kalman gain calculation
	S = param.C*Σ*param.C' + param.V
	K = Σ*param.C'*inv(S)
	res = meas - param.C*x

	#apply kalman update to pos & cov
	xn = x + K*res
	Σn = Σ - K*param.C*Σ

	newmode = mode

	#if a guard is crossed in the update, apply the reset to the resulting pos & cov
	guardresults = [f(xn,t) for f in param.guards[mode]]
	if any(guardresults.<0)
		impact = true
		newmode = findfirst(guardresults.<0)
		Ξ = calc_salt(xn,xn*0,t,mode,newmode,param.dynamics,param.resets,param.guards)
		xn = param.resets[mode][newmode](xn)
		Σn = Ξ*Σn*Ξ'
	end

	return xn,Σn,newmode,impact

end

function skfStep(param,x,Σ,mode,meas,t)
	impact = false
	xint,Σint,modeint,predimpact = skfPredict(param,x,Σ,mode,t)
	xn,Σn,newmode,updateimpact = skfUpdate(param,xint,Σint,modeint,meas,t,predimpact)
	if predimpact || updateimpact
		impact = true
	end
	return xn,Σn,newmode,impact
end