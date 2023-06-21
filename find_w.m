function [W] = find_w(S, U)
    [A, B] = get_lin_model();
    n_points = 10; m = 1.;  g = -9.81; l = 1.; beta = 1 ; max_u  = 10.; dt = .01; 
    a = (m * g * l * cos(0)) / (m * l * l);
    b = l  / (m * l * l);
    temp = S.interval(); x1 = temp(1); x2 = temp(2);
    temp = U.interval(); u1 = temp(1); u2 = temp(2);
    xm = linspace(x1.inf, x1.sup, n_points); 
    ym = linspace(x2.inf, x2.sup, n_points);
    u = linspace(u2.inf, u2.sup, n_points);
    [X1, X2, Um] = meshgrid(xm, ym, u);
    erreur = X2 + dt * (9.81 * sin(X1) - X2 + Um) - (-dt * a * X1 + (1 - dt * beta / (m * l * l)) * X2 + dt * Um);
    erreur_vector = erreur(:);
    min_erreur = min(erreur_vector);
    max_erreur = max(erreur_vector);
    W = zonotope(interval([0.; min_erreur],[0. ; max_erreur]));
end

