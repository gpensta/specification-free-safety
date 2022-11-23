% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
k1 = 4;
k2 = 9;
max_u = 18;
% END PARAMS

X = zonotope(interval([-.1; -.00001],[.1; .00001]));
U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));

max_w = get_max_w(X, U, A, dt, k1, k2);

W = max_w * zonotope(interval([0.; -1.],[0. ; 1.]));
fX = k_step_forward(X, W, A, dt, k1);
bX = k_step_backward(X, U, A, dt, k2);

close all;
f = figure;
hold on;
xlim([-5 5]);
ylim([-5 5]);
initial_set = plot(X, [1,2], 'b', 'linewidth', 2);
wplot = plot(w & U, [1,2], 'k', 'linewidth', 2);
plot(w , [1,2], 'k', 'linewidth', 2);
f = plot(fX, [1,2], 'r', 'linewidth', 2);
b = plot(bX, [1,2], 'g', 'linewidth', 2);
legend([f, b, initial_set, wplot], sprintf('F_%d(S, W=%.2f*U)', k1, max_w/max_u),sprintf('B_%d(T, U = [-%.0f, %.0f])', k2, max_u, max_u), 'Initial=Target', 'W');


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

function res = get_max_w(X, U, A, dt, k1, k2)
    n_points = 200;
    Theta = linspace(0, 2*pi, n_points);
    c = zeros(2, n_points);
    d = ones(n_points, 1);
    for i = 1:length(Theta)
        theta = Theta(i);
        c(1, i) = cos(theta);
        c(2, i) = sin(theta); % n_points direction on the unit circle.
        sumai = [0 0; 0 0];
        l = c(:, i);
        for k = k2:k1+k2-1
            sumai = sumai + A^k;
        end
        invsumai = inv(sumai); % lhs
        % rhs
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
    res =supportFunc(w & U, [0; 1]);
end 