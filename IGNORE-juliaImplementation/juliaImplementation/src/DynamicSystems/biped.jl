###### 2-D Biped Walker Dynamics
# Dynamics equation: H ̈q + B ̇q + G = 0

# p holds system parameters
# p = (L = 1, a1 = 0.375, b1 = 0.125, a2 = 0.175, b2 = 0.375, mH = 0.5, mt = 0.5, ms = 0.05, g = 9.81, slope = 0.06)

# maps coords to joint positions
function q2joints(q, j1=[0,0])
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    q1, q2, q3 = q
    j2 = j1 + L*[cos(q1 + π/2), sin(q1 + π/2)]
    j3 = j2 + (L - b1 - a1)*[cos(-π/2 + q2), sin(-π/2 + q2)]
    j4 = j3 + (b1 + a1)*[cos(-π/2 + q3), sin(-π/2 + q3)]
    joints = [j1, j2, j3, j4]
    return joints
end

function H_matrix_3link(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    lt = a2 + b2
    ls = a1 + b1
    q1, q2, q3 = q
    H11 = ms*a1^2 + mt*(ls+a2)^2 + (mH + ms + mt)*L^2
    H12 = -(mt*b2 + ms*lt)*L*cos(q2 - q1)
    H13 = -ms*b1*L*cos(q3 - q1)
    H22 = mt*b2^2 + ms*lt^2
    H23 = ms*b1*lt*cos(q3 - q2)
    H33 = ms*b1^2
    H = [H11 H12 H13; 
         H12 H22 H23; 
         H13 H23 H33] 
    return H
end

function B_matrix_3link(q, q̇)
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    q1, q2, q3 = q
    lt = a2 + b2
    ls = a1 + b1
    h122 = -(mt*b2 + ms*lt)*L*sin(q1-q2)
    h133 = -ms*b1*L*sin(q1-q3)
    h211 = -h122
    h233 = -ms*lt*b1*sin(q3-q2)
    h311 = -h133
    h322 = -h233
    B = [   0 h122 h133; 
         h211    0 h233; 
         h311 h322    0]
    B = B*q̇
    return B
end

function G_matrix_3link(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    q1, q2, q3 = q
    lt = a2 + b2
    ls = a1 + b1
    g1 = -(ms*a1 + mt*(ls+a2) + (mH+ms+mt)*L)*g*sin(q1)
    g2 = (mt*b2 + ms*lt)*g*sin(q2)
    g3 = ms*b1*g*sin(q3)
    G = [g1; g2; g3]
    return G
end

function H_matrix_2link(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    lt = a2 + b2
    ls = a1 + b1
    q1, q2 = q
    H11 = ms*a1^2 + mt*(ls+a2)^2 + (mH + ms + mt)*L^2
    H12 = -(mt*b2 + ms*(lt+b1))*L*cos(q2 - q1)
    H22 = mt*b2^2 + ms*(lt+b1)^2

    H = [H11 H12; 
         H12 H22] 
    return H
end

function B_matrix_2link(q, q̇)
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    q1, q2 = q
    lt = a2 + b2
    ls = a1 + b1
    h = -(mt*b2 + ms*(lt+b1))*L*sin(q1-q2)
    B = [ 0 h; 
         -h 0]
    B = B*q̇
    return B
end

function G_matrix_2link(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p
    q1, q2 = q
    lt = a2 + b2
    ls = a1 + b1
    g1 = -(ms*a1 + mt*(ls+a2) + (mH+ms+mt)*L)*g*sin(q1)
    g2 = (mt*b2 + ms*(lt+b1))*g*sin(q2)
    G = [g1; g2]
    return G
end

function biped_dynamics_3link(x, u, t)
    q, q̇ = x[1:3], x[4:6]
    H = H_matrix_3link(q)
    B = B_matrix_3link(q, q̇)
    G = G_matrix_3link(q)
    q̈ = -H\(B + G)
    ẋ = [q̇; q̈]
    return ẋ
end

function biped_dynamics_2link(x, u, t)
    q, q̇ = x[1:3], x[4:6]
    H = H_matrix_2link(q[1:2])
    B = B_matrix_2link(q[1:2], q̇[1:2])
    G = G_matrix_2link(q[1:2])
    q̈ = -H\(B + G)
    # q̈ = [q̈; q̈[2]]
    q̈ = [q̈; 0]
    q̇[3] = q̇[2]
    ẋ = [q̇; q̈]
    return ẋ
end

function Q⁺knee(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p   
    q1, q2, q3 = q
    lt = a2 + b2
    ls = a1 + b1
    α, β, γ  = cos(q1-q2), cos(q1-q3), cos(q2-q3)
    Q21 = -(ms*(b1+lt)+mt*b2)*L*cos(α)
    Q22 = mt*b2^2 + ms*(b1+lt)^2
    Q11 = Q21 + mt*(ls+a2)^2 + (mH+ms+mt)*L^2 + ms*a1^2
    Q12 = Q21 + ms*(lt+b1)^2 + mt*b2^2
    Q = [Q11 Q12; Q21 Q22]
    return Q
end

function Q⁻knee(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p   
    q1, q2, q3 = q
    lt = a2 + b2
    ls = a1 + b1
    α, β, γ  = cos(q1-q2), cos(q1-q3), cos(q2-q3)
    Q11 = -(ms*lt+mt*b2)*L*cos(α) - ms*b1*L*cos(β) + (mt+ms+mH)*L^2 + ms*a1^2 + mt*(ls+a2)^2
    Q12 = -(ms*lt+mt*b2)*L*cos(α) + ms*b1*lt*cos(γ) + mt*b2^2 + ms*lt^2
    Q13 = -ms*b1*L*cos(β) + ms*b1*lt*cos(γ) + ms*b1^2
    Q21 = -(ms*lt+mt*b2)*L*cos(α) - ms*b1*L*cos(β)
    Q22 = ms*b1*lt*cos(γ) + ms*lt^2 + mt*b2^2
    Q23 = ms*b1*lt*cos(γ) + ms*b1^2
    Q = [Q11 Q12 Q13; Q21 Q22 Q23]
    return Q
end

function Q⁺heel(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p   
    q1, q2, q3 = q
    lt = a2 + b2
    ls = a1 + b1
    α = cos(q1-q2)
    Q21 = -(ms*(b1+lt)+mt*b2)*L*cos(α)
    Q22 = ms*(b1+lt)^2 + mt*b2^2 
    Q11 = Q21 + mt*(ls+a2)^2 + (mH+ms+mt)*L^2 + ms*a1^2
    Q12 = Q21 + ms*(lt+b1)^2 + mt*b2^2
    Q = [Q11 Q12; Q21 Q22]
    return Q
end

function Q⁻heel(q)
    L, a1, b1, a2, b2, mH, mt, ms, g = p   
    q1, q2, q3 = q
    lt = a2 + b2
    ls = a1 + b1
    α = cos(q1-q2)
    Q12 = -ms*a1*(lt+b1) + mt*b2*(ls+a2)
    Q11 = Q12 + (mH*L+2*mt*(a2+ls)+ms*a1)*L*cos(α)
    Q = [Q11 Q12; Q12 0]
    return Q
end

function kneeReset(x)
    q⁻, q̇⁻ = x[1:3], x[4:6]
    Q⁺, Q⁻ = Q⁺knee(q⁻), Q⁻knee(q⁻)
    q⁺ = q⁻
    q̇⁺ = Q⁺\(Q⁻*q̇⁻)
    q̇⁺ = [q̇⁺; q̇⁺[2]]
    return [q⁺;q̇⁺]
end

function heelReset(x)
    q⁻, q̇⁻ = x[1:3], x[4:6]
    Q⁺, Q⁻ = Q⁺heel(q⁻), Q⁻heel(q⁻)
    q⁺ = [0 1; 1 0; 1 0] * q⁻[1:2]
    q̇⁺ = Q⁺\(Q⁻*q̇⁻[1:2])
    q̇⁺ = [q̇⁺; q̇⁺[2]]
    q̇⁺[2] = q̇⁺[3] = 0
    return [q⁺;q̇⁺]
end

function idReset(x)
    xplus = deepcopy(x)
    return xplus
end

function kneeGuard(x,t)
    q, q̇ = x[1:3], x[4:6]
    q1, q2, q3 = q
    kneeStrikeCheck = (q3 < q2)
    return kneeStrikeCheck
end

function heelGuard(x,t)
    L, a1, b1, a2, b2, mH, mt, ms, g, slope = p  
    q, q̇ = x[1:3], x[4:6]
    q1, q2, q3 = q
    heelStrikeCheck = slope + (q1 + q2)/2
    return heelStrikeCheck
end

function idGuard(x,t)
    return 1
end