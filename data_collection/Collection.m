function Collection(S, T, pairing_avg, iterations)
    costs = zeros(0, 4);

    for iteration = 1:iterations
        [cost_exhaustive, cost_exhaustive_no_expansion, cost_routing, cost_atoms] = GenerateGraph(iteration, pairing_avg, S, T);

        %cost_exhaustive
        %cost_exhaustive_no_expansion
        %cost_routing
        %cost_atoms

        costs(iteration, 1) = cost_exhaustive;
        costs(iteration, 2) = cost_exhaustive_no_expansion;
        costs(iteration, 3) = cost_routing;
        costs(iteration, 4) = cost_atoms;

        fprintf('Iteration: %d\n', iteration);
    end     

    costs        
end
