function Collection()
    pairing_avg = 1.3; 
    S = 2;
    T = 2;

    costs = zeros(0, 3);

    for iteration = 1:2
        [cost_exhaustive, cost_routing, cost_atoms] = GenerateGraph(iteration, pairing_avg, S, T);
        cost_exhaustive
        cost_routing
        cost_atoms
        costs(iteration, 1) = cost_exhaustive;
        costs(iteration, 2) = cost_routing;
        costs(iteration, 3) = cost_atoms;

        fprintf('\n============================\n');
    end     

    costs        
end
