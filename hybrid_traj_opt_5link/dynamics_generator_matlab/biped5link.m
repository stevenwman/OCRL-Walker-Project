clear; clc; close all

syms t g
syms x y t1 t2 t3 t4 t5
syms xd yd t1d t2d t3d t4d t5d
syms xdd ydd t1dd t2dd t3dd t4dd t5dd
syms xt(t) yt(t) t1t(t) t2t(t) t3t(t) t4t(t) t5t(t)
syms f1x f1y f2x f2y
syms u1 u2 u3 u4
syms m1 m2 m3 m4 m5 m6
syms l12 l23 l34 l45 l36

% define individual symbols
q = [x y t1 t2 t3 t4 t5];
qd = [xd yd t1d t2d t3d t4d t5d];
qdd = [xdd ydd t1dd t2dd t3dd t4dd t5dd];

% define derived symbols
qt = [xt yt t1t t2t t3t t4t t5t];
qt = qt(t);
qdt = diff(qt,t);
qddt = diff(qt,t,t);

% params
m = [m1 m2 m3 m4 m5 m6];
u = [u1 u2 u3 u4];

% reminder: q = [x y t1 t2 t3 t4 t5];
r = {};
r{3} = @(q) [q(1); q(2)];
r{2} = @(q) (r{3}(q) + l23*[cos(q(4)); sin(q(4))]);
r{4} = @(q) (r{3}(q) + l34*[cos(q(6)); sin(q(6))]);
r{1} = @(q) (r{2}(q) + l12*[cos(q(3)); sin(q(3))]);
r{5} = @(q) (r{4}(q) + l45*[cos(q(7)); sin(q(7))]);
r{6} = @(q) (r{3}(q) + l36*[cos(q(5)); sin(q(5))]);
rd = cell(size(r));

T = 0;  % kinetic energy
V = 0;  % potential energy
for i=1:length(r)
    r{i} = r{i}(qt);
    rd{i} = diff(r{i},t);
    T = T + 1/2*(m(i)*sum(rd{i}.^2));
    V = V + g*(m(i)*r{i}(2));
end

L = T - V;
% gen_forces = forcing(q,t);
act_angles = [t1t-t2t t2t-t4t t4t-t3t t5t-t4t]
gen_forces = act_angles(t) * u.' 

% idk why it only works when highest order derivative starts first
ind_vars = flip([q qd qdd]);
der_vars = flip([qt qdt qddt]);

eqns = [];
for i=1:length(q)
    eqns = [eqns; 0 == diff(diff(L,diff(qt(i),t)),t) - diff(L,qt(i)) - diff(gen_forces, qt(i))];
    eqns(i) = subs(eqns(i), der_vars, ind_vars);
end

% based on M*qdd + B*qd + G = Xi, N = -B*dq - G + Xi, which we can eval through [q; qd]
[M, N] = equationsToMatrix(eqns, qdd)
state_vars = [q qd];

[B, ~] = equationsToMatrix(N, u)