clear; clc; close all

syms t
syms g
syms x y t1 t2 t3 t4 t5
syms xd yd t1d t2d t3d t4d t5d
syms xdd ydd t1dd t2dd t3dd t4dd t5dd
syms xt(t) yt(t) t1t(t) t2t(t) t3t(t) t4t(t) t5t(t)

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
% l12 = 0.5;
% l23 = l12;
% l34 = l12;
% l45 = l12;
% l36 = l12;

m = ones(6,1);
m = [m1 m2 m3 m4 m5 m6];

% g = 9.81;

% reminder: q = [x y t1 t2 t3 t4 t5];
r = {};
r{3} = @(q) [q(1); q(2)];
r{2} = @(q) (r{3}(q) + l23*[cos(q(4)); sin(q(4))]);
r{4} = @(q) (r{3}(q) + l34*[cos(q(6)); sin(q(6))]);
r{1} = @(q) (r{2}(q) + l12*[cos(q(3)); sin(q(3))]);
r{5} = @(q) (r{4}(q) + l45*[cos(q(7)); sin(q(7))]);
r{6} = @(q) (r{1}(q) + l36*[cos(q(5)); sin(q(5))]);
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
gen_forces = forcing(q,t);
gen_forces = 0;

% idk why it only works when highest order derivative starts first
ind_vars = flip([q qd qdd]);
der_vars = flip([qt qdt qddt]);

eqns = [];
for i=1:length(q)
    eqns = [eqns; 0 == diff(diff(L,diff(qt(i),t)),t) - diff(L,qt(i)) - diff(gen_forces, qt(i))];
    eqns(i) = subs(eqns(i), der_vars, ind_vars);
end

% constraints
cons = sym.empty(4,0);
cons(1) = (0 == rd{1}(1));
cons(2) = (0 == rd{1}(2));
cons(3) = (0);
cons(4) = (0);
cons = cons.';

for i=1:length(cons)
    cons(i) = subs(cons(i), der_vars, ind_vars);
end

[A, b] = equationsToMatrix(cons, qd);
A_der = subs(A, ind_vars, der_vars); % maintain diff()-able terms to diff() again
Ad_der = diff(A_der, t);
Ad = subs(Ad_der, der_vars, ind_vars);

% based on M*qdd + B*qd + G = Xi, N = -B*dq - G + Xi, which we can eval through [q; qd]
[M, N] = equationsToMatrix(eqns, qdd);
state_vars = [q qd];

% solve [M A^T] [q_dd] = [-N     ] equation
%       [A 0  ] [lamb] = [A_d*q_d] 

LHS_mat = [M A.'; A zeros(size(A, 1), size(A, 1))]
RHS_mat = [N; -Ad*qd.']

%For LHS_mat
size(LHS_mat)
LHS_mod = string(LHS_mat);
rows = size(LHS_mod,1);
% Reshape the matrix into a 1D array
LHS_mod = reshape(LHS_mod, 1, []);
% Insert semicolon after every row
semicol = ';';
for i = rows:rows:numel(LHS_mod)
    LHS_mod = [LHS_mod(1:i), {semicol}, LHS_mod(i+2:end)];
end

writematrix(LHS_mod,'LHS_mat.txt','Delimiter','tab');
type 'LHS_mat.txt'

% For RHS_mat
size(RHS_mat)
RHS_mod = string(RHS_mat);
rows = size(RHS_mod,1);
% Reshape the matrix into a 1D array
RHS_mod = reshape(RHS_mod, 1, []);
% Insert semicolon after every row
semicol = ';';
for i = rows:rows:numel(RHS_mod)
    RHS_mod = [RHS_mod(1:i), {semicol}, RHS_mod(i+2:end)];
end

writematrix(RHS_mod,'RHS_mat.txt','Delimiter','tab');
type 'RHS_mat.txt'

function [gen_forces] = forcing(q, time)
    Amp = pi/3;
    omega = 2;
    Kp = 15;
    Kd = 2;
    ref = @(t) Amp*(sin(omega * t));
    th2 = q(4);
    th4 = q(6);
    err = (th2 - th4) - ref(time);
    tau24 = - Kp * err - Kd * diff(err, time);
    forcing = tau24 * (th2 - th4);
    gen_forces = forcing; 
end