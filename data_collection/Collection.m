function [costs,complexity_mat]=Collection(S, T, pairing_avg, iterations, scheme_name)
    clearvars -except S T pairing_avg iterations scheme_name
    clc

    costs = zeros(0, 4);
    classifications = cell(0, 0);
    path_depths = zeros(0, 1); 

    scheme_suffix = 'DIRECTED';

    complexity_mat = zeros(iterations, 4); 
    
    for iteration = 1:iterations
         [cost_exhaustive, shortest_path_count, path_count, cost_exhaustive_no_expansion, cost_routing, cost_atoms, ST_classification, RT, complexity]=...
 GenerateGraph(iteration, pairing_avg, S, T, scheme_name, scheme_suffix);
        %cost_atoms = GenerateGraph(iteration, pairing_avg, S, T, 'NSFNET');

        complexity_mat(iteration, :) = complexity;
        
        %cost_exhaustive
        %cost_exhaustive_no_expansion
        %cost_routing
        %cost_atoms

        %path_depths(iteration, 1) = shortest_path_count;

        costs(iteration, 1) = cost_exhaustive;
        costs(iteration, 2) = cost_exhaustive_no_expansion;
        costs(iteration, 3) = cost_routing;
        costs(iteration, 4) = cost_atoms;

        classifications(iteration, 1) = cellstr(ST_classification);

        fprintf('Iteration: %d\n', iteration);
    end     

    mkdir('costs');
    filename = strcat('costs/costs-', scheme_name, '-', scheme_suffix,...
        '-', num2str(S), '-', num2str(T), '-', num2str(pairing_avg), '-', num2str(iterations));

    diary filename;
    costs
    classifications
    filename_mat = strcat(filename, '.mat');
    save(filename_mat, 'costs', 'classifications')
    diary off;

    exh = costs(:, 1);
    exh_no_exp = costs(:, 2);
    routing = costs(:, 3);
    atoms = costs(:, 4);
    

    DIFF_EXH_ATOMS = sum(exh-atoms);
end
