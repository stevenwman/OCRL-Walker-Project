# 5-link biped model in 2D, adapted from https://web.eecs.umich.edu/~grizzle/biped_book_web/
# dynamics equation: M(q)q̈ + N(q,q') = B*u + J(q,q̇)ᵀλ

params = (m1 = 1,  m2 = 1,  m3 = 1,  m4 = 1,  m5 = 1,  m6 = 1,
		 l12 = 1, l23 = 1, l34 = 1, l45 = 1, l36 = 1, g = 9.81)


function biped5link_kinematics(q, params)
	m1, m2, m3, m4, m5, m6, l12, l23, l34, l45, l36, g = params
	x, y, t1, t2, t3, t4, t5 = q

	r3 = [x y]
	r2 = r3 + l23*[cos(t2) sin(t2)]
	r4 = r3 + l34*[cos(t4) sin(t4)]
	r1 = r2 + l12*[cos(t1) sin(t1)]
	r5 = r4 + l45*[cos(t5) sin(t5)]
	r6 = r3 + l36*[cos(t3) sin(t3)]

	return [r1; r2; r3; r4; r5; r6]
end

function M_matrix(q, params)
	m1, m2, m3, m4, m5, m6, l12, l23, l34, l45, l36, g = params
	x, y, t1, t2, t3, t4, t5 = q

	# [  - m1 - m2 - m3 - m4 - m5 - m6,                                 
	# 0,                                                  
	# l12*m1*sin(t1),                                                                   
	# l23*m1*sin(t2) + l23*m2*sin(t2),                                  
	# l36*m6*sin(t3),                                                                   
	# l34*m4*sin(t4) + l34*m5*sin(t4),                                                  
	# l45*m5*sin(t5)]

	m11 = - m1 - m2 - m3 - m4 - m5 - m6
	m12 = 0 
	m13 = l12*m1*sin(t1)
	m14 = l23*m1*sin(t2) + l23*m2*sin(t2)
	m15 = l36*m6*sin(t3)
	m16 = l34*m4*sin(t4) + l34*m5*sin(t4)
	m17 = l45*m5*sin(t5)

	# [                              0,     
	# - m1 - m2 - m3 - m4 - m5 - m6,                                                 
	# -l12*m1*cos(t1),                                                                 
	# - l23*m1*cos(t2) - l23*m2*cos(t2),                                 
	# -l36*m6*cos(t3),                                                                 
	# - l34*m4*cos(t4) - l34*m5*cos(t4),                                                 
	# -l45*m5*cos(t5)]

	m22 = - m1 - m2 - m3 - m4 - m5 - m6
	m23 = -l12*m1*cos(t1)
	m24 = - l23*m1*cos(t2) - l23*m2*cos(t2)
	m25 = -l36*m6*cos(t3)
	m26 = - l34*m4*cos(t4) - l34*m5*cos(t4)
	m27 = -l45*m5*cos(t5)

	# [                 l12*m1*sin(t1),                   
	# -l12*m1*cos(t1),                 
	# -(m1*(2*l12^2*cos(t1)^2 + 2*l12^2*sin(t1)^2))/2,                                   
	# -(m1*(2*l12*l23*sin(t1)*sin(t2) + 2*l12*l23*cos(t1)*cos(t2)))/2,                                               
	# 0,                                                                                                 
	# 0,                                                               
	# 0]

	m33 = -(m1*(2*l12^2*cos(t1)^2 + 2*l12^2*sin(t1)^2))/2
	m34 = -(m1*(2*l12*l23*sin(t1)*sin(t2) + 2*l12*l23*cos(t1)*cos(t2)))/2
	m35 = 0
	m36 = 0
	m37 = 0

	# [l23*m1*sin(t2) + l23*m2*sin(t2), 
	# - l23*m1*cos(t2) - l23*m2*cos(t2), 
	# -(m1*(2*l12*l23*sin(t1)*sin(t2) + 2*l12*l23*cos(t1)*cos(t2)))/2, 
	# - (m1*(2*l23^2*cos(t2)^2 + 2*l23^2*sin(t2)^2))/2 - (m2*(2*l23^2*cos(t2)^2 + 2*l23^2*sin(t2)^2))/2,                                               
	# 0,                                                                                                 
	# 0,                                                               
	# 0]

	m44 = - (m1*(2*l23^2*cos(t2)^2 + 2*l23^2*sin(t2)^2))/2 - (m2*(2*l23^2*cos(t2)^2 + 2*l23^2*sin(t2)^2))/2
	m45 = 0
	m46 = 0
	m47 = 0

	# [                 l36*m6*sin(t3),                   
	# -l36*m6*cos(t3),                                                               
	# 0,                                                                                                 
	# 0, 
	# -(m6*(2*l36^2*cos(t3)^2 + 2*l36^2*sin(t3)^2))/2,                                                                                                 
	# 0,                                                               
	# 0]

	m55 = -(m6*(2*l36^2*cos(t3)^2 + 2*l36^2*sin(t3)^2))/2
	m56 = 0
	m57 = 0

	# [l34*m4*sin(t4) + l34*m5*sin(t4), 
	# - l34*m4*cos(t4) - l34*m5*cos(t4),                                                               
	# 0,                                                                                                 
	# 0,                                               
	# 0, 
	# - (m4*(2*l34^2*cos(t4)^2 + 2*l34^2*sin(t4)^2))/2 - (m5*(2*l34^2*cos(t4)^2 + 2*l34^2*sin(t4)^2))/2, 
	# -(m5*(2*l34*l45*sin(t4)*sin(t5) + 2*l34*l45*cos(t4)*cos(t5)))/2]

	m66 = - (m4*(2*l34^2*cos(t4)^2 + 2*l34^2*sin(t4)^2))/2 - (m5*(2*l34^2*cos(t4)^2 + 2*l34^2*sin(t4)^2))/2
	m67 = -(m5*(2*l34*l45*sin(t4)*sin(t5) + 2*l34*l45*cos(t4)*cos(t5)))/2

	# [                 l45*m5*sin(t5),                   
	# -l45*m5*cos(t5),                                                               
	# 0,                                                                                                 
	# 0,                                               
	# 0,                                   
	# -(m5*(2*l34*l45*sin(t4)*sin(t5) + 2*l34*l45*cos(t4)*cos(t5)))/2,                 
	# -(m5*(2*l45^2*cos(t5)^2 + 2*l45^2*sin(t5)^2))/2]

	m77 = -(m5*(2*l45^2*cos(t5)^2 + 2*l45^2*sin(t5)^2))/2

	M = [m11 m12 m13 m14 m15 m16 m17;
			 m12 m22 m23 m24 m25 m26 m27;
			 m13 m23 m33 m34 m35 m36 m37;
			 m14 m24 m34 m44 m45 m46 m47;
			 m15 m25 m35 m45 m55 m56 m57;
			 m16 m26 m36 m46 m56 m66 m67;
			 m17 m27 m37 m47 m57 m67 m77]

	return M
end

function N_matrix(q, q̇, params)
	m1, m2, m3, m4, m5, m6, l12, l23, l34, l45, l36, g = params
	x, y, t1, t2, t3, t4, t5 = q
	xd, yd, t1d, t2d, t3d, t4d, t5d = q̇
	 
	n1 = - (m1*(2*l12*cos(t1)*t1d^2 + 2*l23*cos(t2)*t2d^2))/2 - (m5*(2*l34*cos(t4)*t4d^2 + 2*l45*cos(t5)*t5d^2))/2 - l23*m2*t2d^2*cos(t2) - l34*m4*t4d^2*cos(t4) - l36*m6*t3d^2*cos(t3)
	n2 = g*m1 + g*m2 + g*m3 + g*m4 + g*m5 + g*m6 - (m1*(2*l12*sin(t1)*t1d^2 + 2*l23*sin(t2)*t2d^2))/2 - (m5*(2*l34*sin(t4)*t4d^2 + 2*l45*sin(t5)*t5d^2))/2 - l23*m2*t2d^2*sin(t2) - l34*m4*t4d^2*sin(t4) - l36*m6*t3d^2*sin(t3)
	n3 = (m1*(2*l12*t1d*sin(t1)*(yd + l12*t1d*cos(t1) + l23*t2d*cos(t2)) - 2*l12*t1d*cos(t1)*(l12*t1d*sin(t1) - xd + l23*t2d*sin(t2))))/2 - (m1*(2*l12*cos(t1)*(l12*sin(t1)*t1d^2 + l23*sin(t2)*t2d^2) - 2*l12*sin(t1)*(l12*cos(t1)*t1d^2 + l23*cos(t2)*t2d^2) + 2*l12*t1d*sin(t1)*(yd + l12*t1d*cos(t1) + l23*t2d*cos(t2)) - 2*l12*t1d*cos(t1)*(l12*t1d*sin(t1) - xd + l23*t2d*sin(t2))))/2 + g*l12*m1*cos(t1)
	n4 = (m1*(2*l23*t2d*sin(t2)*(yd + l12*t1d*cos(t1) + l23*t2d*cos(t2)) - 2*l23*t2d*cos(t2)*(l12*t1d*sin(t1) - xd + l23*t2d*sin(t2))))/2 - (m1*(2*l23*cos(t2)*(l12*sin(t1)*t1d^2 + l23*sin(t2)*t2d^2) - 2*l23*sin(t2)*(l12*cos(t1)*t1d^2 + l23*cos(t2)*t2d^2) + 2*l23*t2d*sin(t2)*(yd + l12*t1d*cos(t1) + l23*t2d*cos(t2)) - 2*l23*t2d*cos(t2)*(l12*t1d*sin(t1) - xd + l23*t2d*sin(t2))))/2 + g*l23*m1*cos(t2) + g*l23*m2*cos(t2)
	n5 = g*l36*m6*cos(t3)
	n6 = (m5*(2*l34*t4d*sin(t4)*(yd + l34*t4d*cos(t4) + l45*t5d*cos(t5)) - 2*l34*t4d*cos(t4)*(l34*t4d*sin(t4) - xd + l45*t5d*sin(t5))))/2 - (m5*(2*l34*cos(t4)*(l34*sin(t4)*t4d^2 + l45*sin(t5)*t5d^2) - 2*l34*sin(t4)*(l34*cos(t4)*t4d^2 + l45*cos(t5)*t5d^2) + 2*l34*t4d*sin(t4)*(yd + l34*t4d*cos(t4) + l45*t5d*cos(t5)) - 2*l34*t4d*cos(t4)*(l34*t4d*sin(t4) - xd + l45*t5d*sin(t5))))/2 + g*l34*m4*cos(t4) + g*l34*m5*cos(t4)
	n7 = (m5*(2*l45*t5d*sin(t5)*(yd + l34*t4d*cos(t4) + l45*t5d*cos(t5)) - 2*l45*t5d*cos(t5)*(l34*t4d*sin(t4) - xd + l45*t5d*sin(t5))))/2 - (m5*(2*l45*cos(t5)*(l34*sin(t4)*t4d^2 + l45*sin(t5)*t5d^2) - 2*l45*sin(t5)*(l34*cos(t4)*t4d^2 + l45*cos(t5)*t5d^2) + 2*l45*t5d*sin(t5)*(yd + l34*t4d*cos(t4) + l45*t5d*cos(t5)) - 2*l45*t5d*cos(t5)*(l34*t4d*sin(t4) - xd + l45*t5d*sin(t5))))/2 + g*l45*m5*cos(t5)

	# flip sign because how output of "equationsToMatrix" is defined
	N = - [n1; n2; n3; n4; n5; n6; n7]
	return N
end

function B_matrix()
	# q = [x y t1 t2 t3 t4 t5]
	B = [0 0 0 0 0 0 0;
			 0 0 0 0 0 0 0;
			 0 0 1 0 0 0 0;
			 0 0 0 1 0 0 0;
			 0 0 0 0 1 0 0;
			 0 0 0 0 0 1 0;
			 0 0 0 0 0 0 1]
	return B
end

function left_foot_constraint(q, params, fpos)
	f1posX, f1posY, f2posX, f2posY = fpos

	r = biped5link_kinematics(q, params)
	r1 = r[1:2]

	c1 = r1[1] - f1posX
	c2 = r1[2] - f1posY
	c3 = 0
	c4 = 0

	C = [c1; c2; c3; c4]
	return C
end

function right_foot_constraint(q, params, fpos)
	f1posX, f1posY, f2posX, f2posY = fpos

	r = biped5link_kinematics(q, params)
	r2 = r[3:4]

	c1 = 0
	c2 = 0
	c3 = r2[1] - f2posX
	c4 = r2[2] - f2posY

	C = [c1; c2; c3; c4]
	return C
end

function both_foot_constraint(q, params, fpos)
	f1posX, f1posY, f2posX, f2posY = fpos

	r = biped5link_kinematics(q, params)
	r1 = r[1:2]
	r2 = r[3:4]

	c1 = r1[1] - f1posX
	c2 = r1[2] - f1posY
	c3 = r2[1] - f2posX
	c4 = r2[2] - f2posY

	C = [c1; c2; c3; c4]
	return C
end

function J_matrix(q, params, constraint, fpos)
	J  = ForwardDiff.jacobian(_q -> constraint(_q, params, fpos), q)
	return J
end

function kkt_conditions(q, q̇, u, params, h)
	# rigth hand side of the KKT system:
	# [  M(q)  -J(q)ᵀ*h] [q̇ₖ₊₁] = [M(q)*q̇ₖ + h*B*u - h*N]
	# [J(q)*h         0] [   λ] = [               -C(q)]

	M = M_matrix(q, params)
	N = N_matrix(q, q̇, params)
	B = B_matrix()

	kkt_conditions = [M*q̇ + h*B*u - h*N;
						 				-C(q)]
	return kkt_conditions
end

function kkt_jacobian(q, params, constraint, fpos, h)
	# returns the left hand side of the KKT system:
	# [  M(q)  -J(q)ᵀ*h] [q̇ₖ₊₁] = [M(q)*q̇ₖ + h*B*u - h*N]
	# [J(q)*h         0] [   λ] = [               -C(q)]

	M = M_matrix(q, params)
	J = J_matrix(q, params, constraint, fpos)

	kkt_jac = [M   -J'*h;
						 J*h zeros(size(J, 1), size(J, 1))]
	return kkt_jac
end

function kkt_newton_step(q, q̇, u, params, fpos, constraint, h)
	# solve the KKT system for the newton step
	kkt_jac = kkt_jacobian(q, params, constraint, fpos, h)
	kkt_cond = kkt_conditions(q, q̇, u, params, h)

	newton_step = kkt_jac \ kkt_cond
	return newton_step
end

function forward_dynamics(q, q̇, u, params, fpos, constraint, h)
	# solve the KKT system for the newton step
	newton_step = kkt_newton_step(q, q̇, u, params, fpos, constraint, h)

	# update q and q̇
	q̇ₖ₊₁ = newton_step[1:7]
	λ = newton_step[8:end]

	qₖ₊₁ = q + q̇ₖ₊₁*h
	return qₖ₊₁, q̇ₖ₊₁, λ
end

function simulate(q, q̇, u, params, fpos, constraint, h, T)
	q_hist = zeros(T, 7)
	q̇_hist = zeros(T, 7)
	λ_hist = zeros(T, 4)

	q_hist[1, :] = q
	q̇_hist[1, :] = q̇

	for t = 2:T
		q, q̇, λ = forward_dynamics(q, q̇, u, params, fpos, constraint, h)
		q_hist[t, :] = q
		q̇_hist[t, :] = q̇
		λ_hist[t, :] = λ
	end

	return q_hist, q̇_hist, λ_hist
end