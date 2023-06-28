function [A, B] = get_lin_model()
    m = 1.;  g = -9.81; l = 1.; dt = .05; beta = 1;
    a = (m * g * l * cos(0)) / (m * l * l);
    b = l  / (m * l * l);
    A = [1,    dt;
        -dt * a,  1 - dt * beta / (m * l * l)]; 
    B = [0 0;
         0 dt * b];
end

