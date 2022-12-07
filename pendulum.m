% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
% A = [1 0; 0 1];
k1 = 4;
k2 = 9;
% max_u = 14;
% END PARAMS
profile on

for max_u = 7:1:7
    T = zonotope(interval([-.15; -.1],[.15; .1]));
    X = zonotope(interval([-.4; .05],[-.3; .1]));
    U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
    
    [W, w_box] = get_max_w(X, T, U, A, dt, k1, k2);
    if isempty(W.vertices)
        disp("W is empty.");
    else
        close all;
        f = figure('Position', [10 10 1200 780]);
        subplot(1, 2, 1);
        hold on;
        xlim([-1 1]);
        ylim([-8 8]);
        fX = k_step_forward(X, W, A, dt, k1);
        [bX, kr] = k_step_backward(T, U, A, dt, k2);    
        initial_set = plot(X, [1,2], 'm', 'linewidth', .5);
        target_set = plot(T, [1,2], 'g', 'linewidth', .5);
        plot(w_box, [1,2], 'm', 'linewidth', .5);
        draw_streamlines(A, dt, max_u, 7);
        subplot(1, 2, 2);
        hold on;
        title("W and U");
        u = plot(U, [1,2], 'b', 'linewidth', 10.);
        w =plot(W, [1,2], 'r', 'linewidth', 7.);
        legend([w, u], 'w', 'u');
        ylim([-20 20]);
        xlim([-.1 .1]);
%       saveas(f, sprintf("C:\\Users\\kiwin\\Pictures\\article\\pres\\%d.png", floor(max_u)));
end
end
profile viewer

function res = k_step_forward(x, u, A, dt, k)
    for i = 1:k
        x = A*x + dt*u;
        f = plot(x, [1,2], 'r', 'linewidth', .5);
    end
    res = x;
end

function [res, kr] = k_step_backward(x, u, A, dt, k)
    for i = 1:k
        old_x = x;
        x = inv(A) *  (old_x + -dt*u);
        b = plot(x, [1,2], 'b', 'linewidth', 0.5);
        
        
        inter_old = old_x & zonotope(interval([-1.; -8], [1.; 8]));
        inter = x & zonotope(interval([-1.; -8], [1.; 8]));
        
        if inter_old.volume / inter.volume > 0.995
            kr = i;
            break
        end
    end
    res = x;
end

function [w, w_box] = get_max_w(X, T, U, A, dt, k1, k2)
    % Theta = linspace(0, 2*pi, n_points);
    Theta = [0,pi, pi*0.5, pi*1.5];
    c = zeros(2, length(Theta));
    d = ones(length(Theta), 1);
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


function res = clip(x, inf , sup)
    if x < inf
        x = inf;
    elseif x > sup
        x = sup;
    end
    res = x;
end

function res = pd_control(x, sup_u)
    g = -10.;
    l = 1.;
    m = 1.;
    beta = 1.;
    a = -g/l;
    b = -beta/m;
    u = -(a+1)*x(1,1) - (b+2)*x(2, 1);
    res = clip(u, -sup_u, sup_u);
end

function res = move(x, u, A, dt)
        res = A*x + [0; dt]*u;
end

function draw_streamlines(A, dt, max_u, num)
    n = 50;
    [X, Y] = meshgrid(linspace(-1, 1, n), linspace(-8, 8, n));
    [SX, SY] = meshgrid(linspace(-1, 1, num), linspace(-8, 8, num));
    U = ones(n);
    V = ones(n);
    for i = 1:1:n
        for j = 1:1:n
            x1 = X(i, j);
            x2 = Y(i, j);
            x = [x1; x2];
            u = pd_control(x, max_u);
            x = move(x, u, A, dt);
            U(i, j) = x(1, 1) - x1;
            V(i, j) = x(2, 1) - x2;
        end
    end
    slines = streamline(X, Y, U, V, SX, SY);
    set(slines, 'Color', 'k');
end

