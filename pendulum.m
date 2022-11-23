% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];

X = zonotope(interval([-.1; -.00001],[.1; .00001]));
U = 10 * zonotope(interval([0.; -1.],[0. ; 1.]));
k = 0.34364;
W = k*U;

close all;
f = figure;
hold on;
% legend([f, b, i],'F_4(S, W=0.2*U)','B_9(T, U)', 'Initial=Target');



k1 = 4;
k2 = 9;
fX = k_step_forward(X, W, A, dt, k1);
bX = k_step_backward(X, U, A, dt, k2);

% f = plot(fX, [1,2], 'r', 'linewidth', 2);
% b = plot(bX, [1,2], 'g', 'linewidth', 2); 
% 
% 
% xlim([-2 2]);
% ylim([-2 2]);
% legend([f, b, i],'F_4(S, W=0.2*U)','B_9(T, U)', 'Initial=Target');



% --Goal: retrieve k = 0.2 -> w = 1 --Draw W
n_points = 200;
Theta = linspace(0, 2*pi, n_points);
c = zeros(2, n_points);
d = ones(n_points, 1);

for i = 1:length(Theta)
    theta = Theta(i);
    c(1, i) = cos(theta);
    c(2, i) = sin(theta);
    sumai = [0 0; 0 0];
    l = c(:, i);
    for k = k2:k1+k2-1
        sumai = sumai + A^k;
    end
    invsumai = inv(sumai);
    rho_s = supportFunc(X, l);
    rho_aiu = 0;
    for k = 0:k2-1
        rho_aiu = rho_aiu + supportFunc(-A^k*U*dt, l);
    end 
    rho_ak1k2s = supportFunc(A^(k1+k2)*X, l);
    rhs = rho_s + rho_aiu - rho_ak1k2s;
    d(i, 1) = rhs;
end

w = (1/dt)*invsumai*mptPolytope(c', d);

%fprintf("%f\n", supportFunc(poly, [0; 1]));

x = ones(20);
y = linspace(-20, 20, 20);
initial_set = plot(X, [1,2], 'b', 'linewidth', 2);

plot(0*ones(20), linspace(-20, 20, 20));
plot(linspace(-20, 20, 20), 3.43*ones(20));

xlim([-5 5]);
ylim([-5 5]);
wplot = plot(w & U, [1,2], 'k', 'linewidth', 2);
plot(w , [1,2], 'k', 'linewidth', 2);
f = plot(fX, [1,2], 'r', 'linewidth', 2);
b = plot(bX, [1,2], 'g', 'linewidth', 2);
legend([f, b, initial_set, wplot],'F_4(S, W=0.343*U)','B_9(T, U = [-10, 10])', 'Initial=Target', 'W');

function res = k_step_forward(x, u, A, dt, k)
    for i = 1:k
        x = A*x + dt*u;
    end
    res = x;
end

function res = k_step_backward(x, u, A, dt, k)
    for i = 1:k
        x = inv(A) *  (x + -dt*u);
    end
    res = x;
end
