function calc_DxR(reset,x)
    return FD.jacobian(reset,x)
end

function calc_Dxg(g,x,t)
    return FD.gradient(dx -> g(dx, t), x)
end

function calc_A(dyn,x,u,t)
    return FD.jacobian(dx -> dyn(dx,u,t),x)
end

function calc_B(dyn,x,u,t)
    return FD.jacobian(du -> dyn(x,du,t),u)
end

function calc_A_disc(dyn,x,u,t,dt)
	out = dt*calc_A(dyn,x,u,t)
	for i = 1:size(x,1)
		out[i,i] = out[i,i] + 1.0
	end
    return out
end

function calc_B_disc(dyn,x,u,t,dt)
    return dt*calc_B(dyn,x,u,t)
end

function calc_salt(x,u,t,mode1,mode2,dynamics,resets,guards)
    f1 = dynamics[mode1](x,u,t)
    f2 = dynamics[mode2](resets[mode1][mode2](x),u,t)
    DxR = calc_DxR(resets[mode1][mode2],x)
    Dxg = calc_Dxg(guards[mode1][mode2],x,t)
    salt = DxR + ((f2-DxR*f1)*Dxg')/(Dxg'*f1)
    return salt
end