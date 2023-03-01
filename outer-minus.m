% PARAMS
dt = 0.01;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
max_u = 6.;
max_w = .1;
% END PARAMS

close all;
figure; 
hold on;

% Computation

U = dt*max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
W = dt*max_w * zonotope(interval([-1.; -1.],[1. ; 1.]));
T =  zonotope(interval([-.1; -.1],[.1 ; .1]));
h = 1e-3 * [1; 1; 1; 1];
H = [1 0; 0 1; -1 0; 0 -1];


B1 = k_step_backward(T, U, A, 10);
% Br1 = k_robust_backward(T, U, W, A, 10);

G = B1.generators;
c = B1.center;
n = size(G, 2);

[x,fval] = optimize_alpha_c(H, h, G);
disp(x)


gr =(1 - x(1: n))' .* G;
cr = c - x(n+1:n+2);

rb = myMinkDiff(B1, H, h);
rmat = B1.minkDiff(W,'outer');
% PLOT
plot(T, [1,2], 'g', 'linewidth', 1., 'FaceColor', 'none');
plot(B1, [1,2], 'b', 'linewidth', 1.,'FaceColor', 'none');
plot(rb, [1,2], 'r', 'linewidth', 1., 'FaceColor', 'none');
plot(rmat, [1,2], 'm', 'linewidth', 1., 'FaceColor', 'none');


function res = myMinkDiff(Z, H, h)
    G = Z.generators;
    c = Z.center;
    n = size(G, 2);
    [x, fval] = optimize_alpha_c(H, h, G);
    gr =(1 - x(1: n))' .* G;
    cr = c - x(n+1:n+2);
    res = zonotope(cr, gr);
end

function [x,fval] = optimize_alpha_c(H, h, G)
n = size(G, 2);
d = sqrt(sum(G.^2, 1));
% Define the objective function
fun = @(x) -d*log(x(1:n));
% Define the initial guess
x0 = 0.9*ones(n + 2, 1);
% Define the bounds on the optimization variables
lb = -1*ones(1, n + 2);
ub = ones(1, n + 2);
% Define the options for the optimization algorithm
options = optimoptions('fmincon', 'Display', 'iter');
% Call the fmincon function to optimize the objective function subject to the constraints
[x, fval] = fmincon(fun, x0, [], [], [], [], lb, ub,@(x) nonlinear_constraints(x, H, h, G, n), options);
% Define the nonlinear constraints
function [c, ceq] = nonlinear_constraints(x, H, h, G, n)
    % Extract alpha and c from the optimization variables
    alph = x(1:n);
    C = x(n+1:n+2);
    % Define the nonlinear inequality constraint
    c = H*C + abs(H*G)*alph - h;
    % Define the nonlinear equality constraint
    ceq = [];
end
end
function res = k_step_backward(X, U, A, k)
    for i = 1:k
        X = inv(A) *  (X + -1*U);
        %plot(X, [1,2], 'b', 'linewidth', 1.,'FaceColor', 'none');
    end
    res = X;
end

function res = k_robust_backward(X, U, W, A, k)
    for i = 1:k
        X = inv(A) *  (minkDiff(X, W) + -1*U);
        plot(X, [1,2], 'm', 'linewidth', 1.,'FaceColor', 'none');
    end
    res = X;
end