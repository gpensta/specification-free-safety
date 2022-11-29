% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
k1 = 4;
k2 = 20;
% max_u = 14;
% END PARAMS


rollout(A, dt, 6, 25, 1000) 
% supervised_rollout(A, dt, 4, 20, 6, 25, 1000)

function res = supervised_rollout(A, dt, k1, k2, max_u, horizon, n)
    max_x = .9;
    resets = 2 * max_x .*rand(1, n) - max_x;
    T = zonotope(interval([-.15; -.1],[.15; .1]));
    U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
    for i = 1:n
        % disp(i);
        x = [resets(i); 0];
        for t = 1:horizon
             if mod(t, 4) == 1
                X = zonotope(interval([x(1,1); -.01],[x(1,1); .01]));
                W = get_max_w(X, T, U, A, dt, k1, k2);
             end
             u = stochastic_control(x, max_u, .4);
             if not(isempty(W.vertices))
                u = supervision(u, W);
                % disp(W.vertices);
             end
            x = move(x, u, A, dt); 
        end
        if norm(x) < 1.
            resets(i) = 1;
        else
            resets(i) = 0;
        end 
    end
    res = sum(resets)/n; 
end

function res = rollout(A, dt, max_u, T, n)
    max_x = .9;
    resets = 2 * max_x .*rand(1, n) - max_x;
    for i = 1:n
        x = [resets(i); rand()];
        for t = 1:T
             u = stochastic_control(x, max_u, .4);
            %u = pd_control(x, max_u);
            x = move(x, u, A, dt); 
        end
        if norm(x) < 1.
            resets(i) = 1; 
        else
            resets(i) = 0;
        end 
    end
    % res = sum(resets)/n;
    res = sum(resets) / n;
end

function res = supervision(u, w)
    inter = w.vertices;
    sup = inter(2, 2);
    low = inter(2, 1);
    if u < low
        u = low;
    elseif u > sup
        u = sup;
    end
    res = u;
end

function res = k_step_forward(x, u, A, dt, k)
    for i = 1:k
        x = A*x + dt*u;
        f = plot(x, [1,2], 'r', 'linewidth', .5);
    end
    res = x;
end
function res = k_step_backward(x, u, A, dt, k)
    for i = 1:k
        x = inv(A) *  (x + -dt*u);
        if mod(i, 2) == 0
            b = plot(x, [1,2], 'b', 'linewidth', 0.5);
        end
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

function res = stochastic_control(x, max_u, eta)
    g = -10.;
    l = 1.;
    m = 1.;
    beta = 1.;
    a = -g/l;
    b = -beta/m;
    u = -(a+1)*x(1,1) - (b+2)*x(2, 1);
    res = clip(u + eta * (2 * max_u *rand() - max_u), -max_u, max_u);
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
