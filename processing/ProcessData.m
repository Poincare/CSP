function ProcessData(cost_file)
    load(cost_file, '-mat');

    costs

    exh_costs = costs(:, 1);
    exh_no_exp_costs = costs(:, 2);
    routing_costs = costs(:, 3);
    atoms_costs = costs(:, 4);

    routing_comparison = zeros(0, 4);
    atoms_comparison = zeros(0, 3);
    exh_no_exp_comparison = zeros(0, 2);
    exh_comparison = zeros(0, 1);

    j = 1;
    for i = 1:length(exh_costs)
        if routing_costs(i) ~= -1
            routing_comparison(j, 1) = exh_costs(i);
            routing_comparison(j, 2) = exh_no_exp_costs(i);
            routing_comparison(j, 3) = atoms_costs(i);
            routing_comparison(j, 4) = routing_costs(i);
            j = j + 1;

        elseif atoms_costs(i) ~= -1
            atoms_comparison(j, 1) = exh_costs(i);
            atoms_comparison(j, 2) = exh_no_exp_costs(i);
            atoms_comparison(j, 3) = atoms_costs(i);
        elseif exh_no_exp_costs(i) ~= -1
            exh_no_exp_comparison(j, 1) = exh_costs(i);
            exh_no_exp_comparison(j, 2) = exh_no_exp_costs(i);
        elseif exh_costs(i) ~= - 1;
            exh_comparison(j, 1) = exh_costs(i);
        end

    end

    fprintf('Routing feasible:%d\n', length(routing_costs(routing_costs ~= -1)));
    fprintf('2004 feasible: %d\n', length(atoms_costs(atoms_costs ~= -1)));
    fprintf('Exh_no_exp feasible: %d\n', length(exh_no_exp_costs(exh_no_exp_costs ~= -1)));
    fprintf('Exh feasible: %d\n', length(exh_costs(exh_costs ~= -1)));
    
    routing_comparison
    atoms_comparison
    exh_no_exp_comparison
    exh_comparison

    
end
