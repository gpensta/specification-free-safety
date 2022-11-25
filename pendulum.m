% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
k1 = 4;
k2 = 9;
max_u = 9;
% END PARAMS

T = zonotope(interval([-.15; -.1],[.15; .1]));
X = zonotope(interval([.4; -.1],[.6; .1]));
U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
U_0 = 0 * zonotope(interval([0.; -6.67],[0. ; -6.67]));
W = get_max_w(X, T, U, A, dt, k1, k2);
if isempty(W.vertices)
    disp("W is empty.");
else

    % W = max_w * zonotope(interval([0.; -1.],[0. ; 1.]));
    fX = k_step_forward(X, W, A, dt, k1);
    bX = k_step_backward(T, U, A, dt, k2);

    close all;
    f = figure;
    hold on;
    xlim([-5 5]);
    ylim([-5 5]);
    initial_set = plot(X, [1,2], 'm', 'linewidth', 2);
    target_set = plot(T, [1,2], 'b', 'linewidth', 2);
    w = plot(W , [1,2], 'k', 'linewidth', 2);
    f = plot(fX, [1,2], 'r', 'linewidth', 2);
    b = plot(bX, [1,2], 'g', 'linewidth', 2);
    legend([f, b, initial_set, target_set, w], sprintf('F_%d(S, W = [%.2f, %.2f])', k1, w.YData(1), w.YData(2)),sprintf('B_%d(T, U = [-%.2f, %.2f])', k2, supportFunc(U, [0; -1]), supportFunc(U, [0; 1])), 'Initial', 'Target', 'W');
end

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

function w = get_max_w(X, T, U, A, dt, k1, k2)
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
        rho_s = supportFunc(T, l);
        rho_aiu = 0;
        for k = 0:k2-1
            rho_aiu = rho_aiu + supportFunc(-A^k*U*dt, l);
        end 
        rho_ak1k2s = supportFunc(A^(k1+k2)*X, l);
        rhs = rho_s + rho_aiu - rho_ak1k2s;
        d(i, 1) = rhs;
    end

    w_box = (1/dt)*invsumai*mptPolytope(c', d);
    w = w_box & U;
end 