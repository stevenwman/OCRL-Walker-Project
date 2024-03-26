###### 2-D Four Link Biped Walker
## Dynamics functions
# refer to thesis (https://groups.csail.mit.edu/robotics-center/public_papers/Hsu07.pdf)


# maps coords to joint positions

function q2joints(q)
    params = set_params()
    L, a1, b1, a2, b2, mh, mt, ms, g = params[:L], params[:a1], params[:b1], params[:a2], params[:b2], params[:mh], params[:mt], params[:ms], params[:g]
    q1, q2, q3 = q
    j1 = [0,0]
    j2 = L*[cos(q1 + π/2), sin(q1 + π/2)]
    j3 = j2 + (L - b1 - a1)*[cos(-π/2 + q2), sin(-π/2 + q2)]
    j4 = j3 + (b1 + a1)*[cos(-π/2 + q3), sin(-π/2 + q3)]
    joints = [j1, j2, j3, j4]
    return joints
end

function unlocked_dynamics(state, p, t)
    params = set_params()
    L, a1, b1, a2, b2, mh, mt, ms, g = params[:L], params[:a1], params[:b1], params[:a2], params[:b2], params[:mh], params[:mt], params[:ms], params[:g]

    lt = a2 + b2
    ls = a1 + b1
    
    # State variables - angles relative to vertical axis
    q1, q2, q3 = state[1:3]
    dq1, dq2, dq3 = state[4:6]

    H11 = ms * a1^2 + mt * (ls + a2)^2 + (mh + ms + mt) * L^2
    H12 = -(mt * b2 + ms * lt) * L * cos(q2 - q1)
    H13 = -ms * b1 * L * cos(q3 - q1)
    H22 = mt * b2^2 + ms * lt^2
    H23 = ms * lt * b1 * cos(q3 - q2)
    H33 = ms * b1^2

    h122 = -(mt * b2 + ms * lt) * L * sin(q1 - q2)
    h133 = -ms * b1 * L * sin(q1 - q3)
    h211 = -h122
    h233 = ms * lt * b1 * sin(q3 - q2)
    h311 = -h133
    h322 = -h233

    H = [H11 H12 H13; H12 H22 H23; H13 H23 H33]
    B = [0 (h122 * dq2) (h133 * dq3); (h211 * dq1) 0 (h233 * dq3); (h311 * dq1) (h322 * dq2) 0]
    G = [-(ms * a1 + mt * (ls + a2) + (mh + ms + mt) * L) * g * sin(q1); 
          (mt * b2 + ms * lt) * g * sin(q2); 
          ms * b1 * g * sin(q3)]

    ddq = H \ -(B * [dq1; dq2; dq3] + G)
    dx = [dq1; dq2; dq3; ddq...]
    return dx
end

function locked_dynamics(state, p, t)
    params = set_params()
    L, a1, b1, a2, b2, mh, mt, ms, g = params[:L], params[:a1], params[:b1], params[:a2], params[:b2], params[:mh], params[:mt], params[:ms], params[:g]

    lt = a2 + b2
    ls = a1 + b1

    q1, q2, q3 = state[1:3]
    dq1, dq2, dq3 = state[4:6]

    H11 = ms * a1^2 + mt * (ls + a2)^2 + (mh + ms + mt) * L^2
    H12 = -(mt * b2 + ms * (lt + b1)) * L * cos(q2 - q1)
    H22 = mt * b2^2 + ms * (lt + b1)^2

    h = (mt * b2 + ms * (lt + b1)) * L * sin(q1 - q2)
    b = h * dq2

    H = [H11 H12; H12 H22]
    B = [0 b; -b 0]
    G = [-(ms * a1 + mt * (ls + a2) + (mh + ms + mt) * L) * g * sin(q1); 
         (mt * b2 + ms * (lt + b1)) * g * sin(q2)]

    ddq = H \ -(B * [dq1; dq2] + G)
    ddq = [ddq; ddq[2]]  # Assuming ddq3 = ddq2 for locked knee

    dx = [dq1; dq2; dq3; ddq...]
    return dx
end

# Reset Functions
function kneestrike_reset(x)
    params = set_params()
    L, a1, b1, a2, b2, mh, mt, ms, g = params[:L], params[:a1], params[:b1], params[:a2], params[:b2], params[:mh], params[:mt], params[:ms], params[:g]

    lt = a2 + b2
    ls = a1 + b1
    
    xplus = deepcopy(x)

    q1 = x[1]
    q2 = x[2]
    q3 = x[3]

    dq1 = x[4]
    dq2 = x[5]
    dq3 = x[6]

    q_minus = [q1, q2, q3]
    dq_minus = [dq1, dq2, dq3]

    alpha = cos(q1 - q2)
    beta = cos(q1 - q3)
    gamma_val = cos(q2 - q3)

    Q11_minus = -(ms * lt + mt * b2) * L * cos(alpha) - ms * b1 * L * cos(beta) + ms * a1^2 + mt * (ls + a2)^2
    Q12_minus = -(ms * lt + mt * b2) * L * cos(alpha) + ms * b1 * lt * cos(gamma_val) + mt * b2^2 + ms * lt^2
    Q13_minus = -ms * b1 * L * cos(beta) + ms * b1 * lt * cos(gamma_val) + mt * b2^2 + ms * lt^2
    Q21_minus = -(ms * lt + mt * b2) * L * cos(alpha) - ms * b1 * L * cos(beta)
    Q22_minus = ms * b1 * lt * cos(gamma_val) + ms * lt^2 + mt * b2^2
    Q23_minus = ms * b1 * lt * cos(gamma_val) + ms * b1^2

    Q21_plus = -(ms * (b1 + lt) + mt * b2) * L * cos(alpha)
    Q22_plus = ms * (lt + b1)^2 + mt * b2^2
    Q11_plus = Q21_plus + mt * (ls + a2)^2 + (ms + mt + mh) * L^2 + ms * a1^2
    Q12_plus = Q21_plus + ms * (lt + b1)^2 + mt * b2^2

    Q_plus = [Q11_plus Q12_plus; Q21_plus Q22_plus]
    Q_minus = [Q11_minus Q12_minus Q13_minus; Q21_minus Q22_minus Q23_minus]

    q_plus = q_minus
    dq_plus = Q_plus \ (Q_minus * q_minus)
    dq_plus = [dq_plus..., dq_plus[2]]

    xplus = [q_plus..., dq_plus...]
    
    return xplus
end

function heelstrike_reset(x)
    params = set_params()
    L, a1, b1, a2, b2, mh, mt, ms, g = params[:L], params[:a1], params[:b1], params[:a2], params[:b2], params[:mh], params[:mt], params[:ms], params[:g]

    lt = a2 + b2
    ls = a1 + b1
    
    xplus = deepcopy(x)

    q1 = x[1]
    q2 = x[2]
    q3 = x[3]  # equal to q2 for compass gait
    q_minus = [q1, q2]

    dq1 = x[4]
    dq2 = x[5]
    dq3 = x[6]  # equal to dq2 for compass gait

    # xplus[1] = q2
    # xplus[2] = q1
    # xplus[3] = q1

    # xplus[4] = dq2
    # xplus[5] = dq1
    # xplus[6] = dq1

    alpha = cos(q1 - q2)

    Q12_minus = -(ms * a1) * (lt + b1) + mt * b2 * (ls + a2)
    Q11_minus = Q12_minus + (mh * L + 2 * mt * (a2 + ls) + ms * a1) * L * cos(alpha)
    Q21_minus = Q12_minus
    Q22_minus = 0

    Q21_plus = -(ms * (b1 + lt) + mt * b2) * L * cos(alpha)
    Q22_plus = ms * (lt + b1)^2 + mt * b2^2
    Q11_plus = Q21_plus + (mh + ms + mt) * L^2 + ms * a1^2 + mt * (a2 + ls)^2
    Q12_plus = Q21_plus + ms * (b1 + lt)^2 + mt * b2^2

    Q_plus = [Q11_plus Q12_plus; Q21_plus Q22_plus]
    Q_minus = [Q11_minus Q12_minus; Q21_minus Q22_minus]

    q_plus = [0 1; 1 0; 1 0] * q_minus
    # q_plus = [1 0; 0 1; 0 1] * q_minus

    dq_plus = Q_plus \ (Q_minus * q_minus)
    dq_plus = [dq_plus..., dq_plus[2]]

    x_plus = [q_plus..., dq_plus...]

    return x_plus
end

function idreset(x)
    xplus = deepcopy(x)
    
    return xplus
end



# Guard functions - identity guard exists only
# to never trigger & act as a placeholder for
# unreachable transitions

function kneestrike_guard(x, t)
    q2 = x[2]
    q3 = x[3]
    
    return -(q3 - q2)
end

function heelstrike_guard(x, t)
    params = set_params()
    gamma = params[:gamma]

    q1 = x[1]
    q2 = x[2]

    # return (q2 + gamma)
    # return q1 + q2
    return gamma + (q1 + q2)/2
end

function idguard(x,t)
    # xplus = deepcopy(x)
    
    # return xplus
    return 1
end

# Function to set biped walker parameters - 
function set_params()
    return Dict(
        :L => 1,        
        :a1 => 0.375,
        :b1 => 0.125,  
        :a2 => 0.175, 
        :b2 => 0.325,   
        :mh => 0.5,    
        :mt => 0.5,    
        :ms => 0.05,   
        :g => 9.81,    
        :gamma => 0.0504
    )
end