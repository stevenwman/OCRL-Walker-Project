###### 2-D Bouncing Ball Dynamics
## Dynamics function
function balldynamics(x,p,t)
    dx = [x[3] + p[1],x[4] + p[2],0 + p[3],-9.8 + p[4]]
end

# Reset Functions
function groundreset(x)
    xplus = deepcopy(x)
    xplus[4] = -0.8*x[4]
    return xplus
end

function airreset(x)
    xplus = deepcopy(x)
    return xplus
end
function idreset(x)
    xplus = deepcopy(x)
    return xplus
end

# Guard Functions - identity guard exists only
# to never trigger & act as a placeholder for
# unreachable transitions

function groundguard(x,t)
    return x[2]
end
function airguard(x,t)
    return x[4]
end
function idguard(x,t)
    return 1
end