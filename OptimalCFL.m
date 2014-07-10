function OptimalCFL(V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
    while 1
        vars = CFL(V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
        
    end
end

function vars_to_mat(vars, V, P, T)
    f = zeros(V, V, P, T);
    x = zeros(V, V, P);
    
    for i = 1:V
        for j = 1:V
            for s = 1:P
                for t = 1:T
                    f(i, j, s, t) = vars(
end