%E: number of edges
%P: number of flows
%T: number of terminal nodes
%EE: adjacency matrix
function OptimalCFL(V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
    global z;
    global fea_z;
    global fea_x;
    global fea_f;
    global fea_cost;
    global fea_idx;
    global min_cost;
    global NO_CHANGE_P


    %sets f values for path switching
    %on: the terminal inputs that are on for this iteration
    function pick_paths(vars, prob_matrix, TS, V, P, T, p, on, t, mq)
        if p > P
            fprintf('Path: ');

            on
            %for all of the edges, set them to 0 and NO_CHANGE_P in prob_matrix
            for s = 1:P
                inputs = TS{s};
                for i = 1:length(inputs)
                    index = linear_index_f(inputs(i), TT(t), s, t, V, P, T);
                    vars(index) = 0;
                    prob_matrix(index, :) = 0;
                    prob_matrix(index, 1) = NO_CHANGE_P; 
                end 
            end
 
            %for all of the on edges, we set the probability so that they
            %cannot change throughout the iterations of CFL
            for i = 1:length(on)
                fprintf('on(i): %d, t: %d, i: %d, t: %d\n', on(i), t, i, t);
                index = linear_index_f(on(i), TT(t), i, t, V, P, T);
                vars(index) = 1;
                prob_matrix(index, :) = 0;
                prob_matrix(index, 2) = NO_CHANGE_P;
            end
            
            z = zeros(V, V);
            [vars, p_cell_mat] = CFL(vars, clause_mat, prob_matrix, V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT);
            p_filename = strcat('p_cell_mats/', 'p_cell_mat.mat');
            save(p_filename, 'p_cell_mat');
            fprintf('SAVED p_cell_mat\n');
            
            if length(vars) == 0
                return
            end

            cost = sum(sum(z(z ~= -1)))
            [f, x] = vars_to_mat(vars, V, P, T);
            fea_x(:, :, :, fea_idx) = x;
            fea_f(:, :, :, :, fea_idx) = f;
            fea_z(:, :, fea_idx) = z;
            fea_cost(fea_idx) = cost
            fea_cost
            fea_idx = fea_idx + 1;
            disp_vars(vars, V, P, T, E, edge);
            return
        end 

        inputs = TS{p};

        for i = 1:length(inputs)
            %prevent overlapping terminal edges
            if ~(any(on == inputs(i)))
                ons = on;
                ons(p) = inputs(i);
                pick_paths(vars, prob_matrix, TS, V, P, T, p + 1, ons, t, mq+1);
            end
        end
    end

    S=P;
    vars = zeros(1, V*V*P*T + V*V*P);
    
    %number of variables
    N = length(vars);
    
    %possible values of variables
    D = 2;
    
    %if this value is found in the probability matrix,
    %the value of the variable cannot be changed
    %global NO_CHANGE_P
    NO_CHANGE_P = 2;

    fea_idx = 1;
    %infinity initially
    min_cost = Inf;
    fea_x = zeros(V, V, P, V);
    fea_f = zeros(V, V, P, T, V);
    fea_z = zeros(V, V, V);
    fea_cost = zeros(1, V);
    
    %vector of variable values
    %f: first V* V * P * T entries
    %x: next V * V * P entries
    %vars(:) = -1;
    vars(1:V*V*P*T) = -1;
    vars = remove_nonexistent_edges(vars, V, P, T, EE);
    %vars = explore_vars_x(vars, EE, SS, V, P, T, E, TT, ST, edge);
    vars = explore_vars_f(vars, EE, SS, V, P, T, E, TT, ST, edge);
    vars_after = length(vars(vars >= 0));

    %clause matrix
    %i,j is marked as 1 if variable i participates in clause j
    clause_mat = generate_clause_mat(IV, OV, V, P, T, ST);
 
    %all x variables participate in clause 5
    clause_mat(V*V*P*T+1:N, 4*V+1) = 1;

    %all f variables participate in the cost clause
    clause_mat(1:V*V*P*T, 4*V+2) = 1;      

    N
    %probability matrix
    p = zeros(N, 2);
    p(1:N, 1:D) = 1/D;
  
    [vars, p, clause_mat] = set_source_links(vars, p, clause_mat, V, P, T, EE, SS);

    iter_counter = 0;
    last_time = cputime;

    [vars, p_cell_mat, iter_counter] = CFL(vars, clause_mat, p, V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
    save('p_cell_mat_fun.mat', '-mat');
    iter_counter
    return

    for t = 1:T
        terminals = cell(1, P);
        perm_eyes = cell(1, P);
        for s = 1:S
            if ST(s, t) == 1
                terminal_inputs = get_terminal_inputs(vars, s, t, V, P, T, E, IV, SS, TT);
                terminals{s} = terminal_inputs;
                width = length(terminal_inputs);
                perm_eyes{s} = eye(width) * reshape(terminal_inputs, width, 1);
            end
        end

        z = zeros(V, V);
        vars = CFL(vars, clause_mat, p, V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT); 
        [f, x] = vars_to_mat(vars, V, P, T);
        fea_f(:, :, :, :, fea_idx) = f;
        fea_x(:, :, :, fea_idx) = x;
        fea_z(:, :, fea_idx) = z;
        fea_cost(fea_idx) = sum(sum(z));

        fea_idx = 2;
        pick_paths(vars, p, terminals,  V, P, T, 1, zeros(1, P), t, 1); 
    end

    fea_cost
    min_cost =  min(fea_cost(fea_cost ~= 0));
    min_id = find(fea_cost == min_cost);
    fea_z(:, :, min_id)
    fea_f(:, :, :, :, min_id)
    
    min_cost

    return

    %while 1
    %    z = zeros(V, V);
    %    vars = CFL(vars, V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT);
    %    %implies no satisfying soln. found
    %    if length(vars) == 0
    %        break;
    %    end

    %    [f, x] = vars_to_mat(vars, V, P, T);
    %    fea_x(:, :, :, fea_idx) = x;
    %    fea_f(:, :, :, :, fea_idx)  = f;
    %    fea_z(:, :, fea_idx) = z;
    %    min_cost = sum(sum(z(z ~= -1)));
    %    fea_cost(:, fea_idx) = min_cost;

    %    min_cost       
    %end
end

function [vars, p]=ReduceSearch(vars, p, v, V, P, T, OV, SS)
    i = v;
    ov = OV{v};

    if sum(ov) == 0
        return
    end

    for ovi = 1:length(ov)  
        j = ov(ovi);

        found_one = 0;
        for s = 1:P
            for t = 1:T
                if get_f_from_vars(vars, i, j, s, t, V, P, T) == 1
                    found_one = 1;
                    break;
                end
            end

            if found_one
                vars = set_x_in_vars(vars, 1, i, j, s, V, P, T);
                p(linear_index_x(i, j, s, V, P, T), 1) = 0;
                p(linear_index_x(i, j, s, V, P, T), 2) = 1;
            end
        end
    end
end

%this function takes the variable group defined by nodes i, j and terminal
%t and the random variable value e in order to pack those values into vars
function vars=unpackVariableGroup(vars, i, j, t, e, D, V, P, T)
    for p = 1:P
        %zero all the variables in this variable group
        vars = set_f_in_vars(vars, 0, i, j, p, t, V, P, T);
    end
    
    %set the variable to which "e" refers as 1
    p = e + 1;
    %if p = 0, e = -1, implying that all the variables in this group
    %should be zero
    if p ~= 0
        vars = set_f_in_vars(vars, 1, i, j, p, t, V, P, T);
    end
end

%randomizes a variable group's value given a variable group probability
%table
function randomizeVariableGroup(i, j, t, D)
    global p_variable_group;
    global variable_group_vals;
    
    r = rand;
    low = 0;
    for e = 0:D-1
        prob = p_variable_group(i, j, t, e + 1);
        if (r >= low) && (r <= low + prob)
            variable_group_vals(i, j, t) = e;
            return
        end
        
        low = low + prob;
    end
end

function satisfied=checkFlowConservation(node, var_index, vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied)
    if clause_mat(var_index, node) == 1 && satisfied
        if flow_conservation(node, vars, V,P,T,OV,IV,sigma,ST,E,edge) ~= 1
            %fprintf('Flow conservation failed.\n\n');
            satisfied = 0;
            return;
        end  
    end
end

function satisfied=checkFlowX(node,var_index, vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied)
    if clause_mat(var_index, V+node) == 1 && satisfied
        %checking clause 3
        if flow_x(node, vars, OV, V, P, T, E, ST, edge) ~= 1
            %fprintf('Flow x failed. i: %d, c: %d, mat: %d, var_value: %d\n\n', i, c, clause_mat(i,c), vars(i));                        satisfied = 0;
        end
    end    
end

function satisfied=checkX(node,var_index, vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,TT,satisfied)
    if clause_mat(var_index, 2*V+node) == 1 && satisfied
        if checkx(node, vars, V, P, T, ST, IV, TT) ~= 1
            %fprintf('Checkx failed.\n')
            satisfied = 0;
        end
    end    
end

function satisfied=checkFlowLimit(node,var_index, vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied)
    if clause_mat(var_index, 3*V+node) == 1 && satisfied
        if flow_limit(node, vars, V, P, T, E, edge, ST, OV) ~= 1
            %fprintf('Flow limit failed.\n')
            satisfied = 0;
        end
    end    
end

function satisfied=checkBeta(node,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied)
    if clause_mat(var_index, 4*V+1) == 1 && satisfied
        %checking clause 5
        if checkbeta(vars, V, P, T, IV, E, edge) ~= 1
            %fprintf('Checkbeta failed.\n')
            satisfied = 0;
        end
    end    
end

function satisfied=checkCostClause(node,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied)
    if clause_mat(var_index, 4*V+2) == 1 && satisfied
        if cost_clause(vars, V, P, T, E, edge, ST, OV) ~= 1
            satisfied = 0;
            cost_clause_failure = cost_clause_failure + 1;
            if cost_clause_failure > cost_clause_thresh
                vars = [];
                return
            end
        end    
    end    
end

function satisfied=checkClauses(i, j, var_index, vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,TT)
    satisfied = 1;
    satisfied = checkFlowConservation(i,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkFlowConservation(j,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkFlowX(i,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkFlowX(j,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkX(i,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,TT,satisfied);
    satisfied = checkX(j,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,TT,satisfied);
    satisfied = checkFlowLimit(i, var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkFlowLimit(j, var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkBeta(i, var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkBeta(j, var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkCostClause(i,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
    satisfied = checkCostClause(j,var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,satisfied);
end

function [vars,p_cell_mat,iter_counter]=CFL(vars, clause_mat, p, V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
    global NO_CHANGE_P

    global z
    z = zeros(V, V);
    
    %meant for Octave (not sure about matlab syntax)
    more off

    %global beta matrix (will not work when algorithm is
    %run as distributed)
    global beta_mat
    beta_mat = zeros(V, V, V);

    %number of variables
    N = length(vars);
 
    %possible values of variables
    D = 2;
    
    %number of clauses
    %(17): flow in {0, 1} - every f participates
    C = 4*V+3;

    %vector of variable values
    %f: first V* V * P * T entries
    %x: next V * V * P entries
    %vars(:) = -1;
    vars(1:V*V*P*T) = -1;
    vars = remove_nonexistent_edges(vars, V, P, T, EE);
    %vars = explore_vars_x(vars, EE, SS, V, P, T, E, TT, ST, edge);
    vars = explore_vars_f(vars, EE, SS, V, P, T, E, TT, ST, edge);
    vars_after = length(vars(vars >= 0));

    %clause matrix
    %i,j is marked as 1 if variable i participates in clause j
    clause_mat = generate_clause_mat(IV, OV, V, P, T, ST);
 
    %all x variables participate in clause 5
    clause_mat(V*V*P*T+1:N, 4*V+1) = 1;

    %all f variables participate in the cost clause
    clause_mat(1:V*V*P*T, 4*V+2) = 1;      

    global D_variable_group
    global p_variable_group
    global N_variable_group
    global variable_group_vals
    
    init_variable_groups(V, P, T, ST);
    D_variable_group
    
    [vars, p, clause_mat] = set_source_links(vars, p, clause_mat, V, P, T, EE, SS);

    iter_counter = 0;
    last_time = cputime;

    %number of times the cost clause is unsatisfied
    %after the threshold, we give up 
    cost_clause_failure = 0;
    cost_clause_thresh = 100;

    p_cell_mat = cell(0, 1);

    %design parameter
    b = 0.3;
    
    N_xvars = V*V*P;
    
    while 1
        done = 1;
        globally_satisfied = 1;

        for i = 1:V
            for j = 1:V
                for s = 1:P
                    fixed = 0;
                    var_index = linear_index_x(i, j, s, V, P, T);
                    if vars(var_index) == -1
                        continue;
                    end

                    p
                    if any(p(var_index, :) == NO_CHANGE_P) ~= 0 
                        fixed = 1;
                    end

                    satisfied = checkClauses(i, j, var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,TT);
                    if satisfied == 1
                        p(i, :) = 0;
                        p(i, vars(var_index) + 1) = 1;
                        continue;
                    %otherwise, we have to interpolate the distribution
                    else
                        globally_satisfied = 0;
                        if ~fixed
                            t = p(i,vars(i)+1);
                            for val = 1:D
                                p(var_index, val) = (1-b)*(p(var_index,val)) + (b/(D-1));
                            end

                            %i
                            %vars(i) 
                            %(1-b)*t

                            p(i, vars(i) + 1) = (1-b)*t;
                        end
                    end
                end
            end
        end
        
        for v = 1:V
            for i = 1:V
                for j = 1:V
                    for t = 1:T
                        D = D_variable_group(i, j, t);
                        
                        %randomize the variable group's value
                        randomizeVariableGroup(i, j, t, D);
                        
                        e = variable_group_vals(i, j, t);
                        
                        %unpack the variable group into vars
                        unpackVariableGroup(vars, i, j, t, e, D, V, P, T);
                        
                        satisfied = 1;
                        %check feasibility
                        for p = 1:D-1
                            var_index = linear_index_f(i, j, p, t, V, P, T);
                            satisfied = checkClauses(i, j, var_index,vars,clause_mat,V,P,T,OV,IV,sigma,ST,E,edge,TT);
                        end
                        
                        %update the probability distribution for this
                        %variable group.
                        
                        if satisfied
                            p_variable_group(i, j, t, :) = 0;
                            p_variable_group(i, j, t, e + 1) = 1;
                            continue;                            
                        else
                            globally_satisfied = 0;
                            prev = p_variable_group(i,j,t,e+1);
                            for k = 0:D-1
                                p_variable_group(i,j,t,k+1) = (1-b)*(p_variable_group(i,j,t,k+1)) + (b/(D-1));
                            end

                            p_variable_group(i, j, t, e+1) = (1-b)*prev;
                        end
                    end
                end
            end        
        end
        
        if globally_satisfied
            p_size = size(p);

            p_cell_mat = [p_cell_mat; p];
            save('p_cell_mat.mat', 'p_cell_mat');
            break 
        end
    end
    
%     while 1
%         done = 1;
%         for v = 1:V
%             for i = 1:N
%                 if rem(iter_counter, 100) == 0
%                     p_cell_mat = [p_cell_mat; p];
%                 end
%                 
%                 %V
%                 %P
%                 %T
%                 %index_1_1 = linear_index_f(3, 4, 1, 1, V, P, T)
%                 %index_1_2 = linear_index_f(3, 4, 1, 2, V, P, T)
%                 %index_2_1 = linear_index_f(3, 4, 2, 1, V, P, T)
%                 %index_2_2 = linear_index_f(3, 4, 2, 2, V, P, T)
%                 
%                 %p_3_4_1_1 = p(linear_index_f(3, 4, 1, 1, V, P, T), :)
%                 %p_3_4_1_2 = p(linear_index_f(3, 4, 1, 2, V, P, T), :)
%                 %p_3_4_2_1 = p(linear_index_f(3, 4, 2, 1, V, P, T), :)
%                 %p_3_4_2_2 = p(linear_index_f(3, 4, 2, 2, V, P, T), :)
%                 
%                 fixed = 0;
% 
%                 %check if the variable is supposed to be considered
%                 if vars(i) == -1
%                     continue;    
%                 end
% 
%                 if any(p(i, :) == NO_CHANGE_P) ~= 0 
%                     fixed = 1;
%                 end
% 
%                 %design variable
%                 b = 0.1;
% 
%                 %evaluate the clauses and see if they are
%                 %satisfied
%                 satisfied = 1;
%                 %disp_vars(vars, V, P, T, E, edge)
% 
% 
%                 [vars, p] = ReduceSearch(vars, p, v, V, P, T, OV, SS);
%                 
%                 %for v = 1:v
%                 %     ov = OV{v};
%                 %     if sum(ov) ~= 0
%                 %         for ovi = 1:length(ov)
%                 %             j = ov(ovi);
%                 %             j
%                 %             for s = 1:P
%                 %                 fprintf('x, i: %d, j;%d, s: %d, t: %d, val: %d\n',    v, j, s, get_x_from_vars(vars, v, j, s, V, P, T));
%                 %                 for t = 1:T
%                 %                     fprintf('f, i: %d, j: %d, p: %d, t: %d, val: %d\n', v, j, s, t, get_f_from_vars(vars, v, j, s, t, V, P, T));
%                 %                 end
%                 %             end
%                 %         end
%                 %     end
% 
%                 % end
% 
%                 for c = 1:V
%                     if clause_mat(i, c) == 1 && satisfied
%                         if flow_conservation(c, vars, V,P,T,OV,IV,sigma,ST,E,edge) ~= 1
%                             %fprintf('Flow conservation failed.\n\n');
%                             satisfied = 0;
%                             break;
%                         end  
%                     end
%                 end
% 
%                 for c = 1:V
%                     if clause_mat(i, V+c) == 1 && satisfied
%                         %checking clause 3
%                         if flow_x(c, vars, OV, V, P, T, E, ST, edge) ~= 1
%                             %fprintf('Flow x failed. i: %d, c: %d, mat: %d, var_value: %d\n\n', i, c, clause_mat(i,c), vars(i));                        satisfied = 0;
%                         end
%                     end
%                 end
% 
%                 
%                for c = 1:V
%                    if clause_mat(i, 2*V+c) == 1 && satisfied
%                         if checkx(c, vars, V, P, T, ST, IV, TT) ~= 1
%                             %fprintf('Checkx failed.\n')
%                             satisfied = 0;
%                         end
%                     end
%                end
%                
%                for c = 1:V
%                    if clause_mat(i, 3*V+c) == 1 && satisfied
%                         if flow_limit(c, vars, V, P, T, E, edge, ST, OV) ~= 1
%                             %fprintf('Flow limit failed.\n')
%                             satisfied = 0;
%                         end
%                     end
%                end
%                
%                 %if clause_mat(i, 2) == 1 && satisfied 
%                 %    %checking clause2
%                 %    if flow_conservation(vars,V,P,T,OV,IV,sigma,ST) ~= 1
%                 %        %fprintf('Flow conservation failed.\n')
%                 %        satisfied = 0;
%                 %    end
%                 %end
%                 
%                 if clause_mat(i, 4*V+1) == 1 && satisfied
%                     %checking clause 5
%                     if checkbeta(vars, V, P, T, IV, E, edge) ~= 1
%                         %fprintf('Checkbeta failed.\n')
%                         satisfied = 0;
%                     end
%                 end
% 
%                 %check cost clause
%                 if clause_mat(i, 4*V+2) == 1 && satisfied
%                     if cost_clause(vars, V, P, T, E, edge, ST, OV) ~= 1
%                         satisfied = 0;
%                         cost_clause_failure = cost_clause_failure + 1;
%                         if cost_clause_failure > cost_clause_thresh
%                             vars = [];
%                             return
%                         end
%                     end    
%                 end
% 
% 
%                 %fprintf('-----------\n\n\n')
% 
%                 %fprintf('Before: \n');
%                 %p           
%                 %fprintf('Satisfied: %d, i: %d, vars(i): %d\n', satisfied, i, vars(i)); 
%                 
%                 %if it is satisfied, we can clear the probabilities
%                 if satisfied == 1
%                     p(i, :) = 0;
%                     p(i, vars(i) + 1) = 1;
%                     continue;
%                 %otherwise, we have to interpolate the distribution
%                 else
%                     if ~fixed
%                         t = p(i,vars(i)+1);
%                         for j = 1:D
%                             p(i, j) = (1-b)*(p(i,j)) + (b/(D-1));
%                         end
%                   
%                         %i
%                         %vars(i) 
%                         %(1-b)*t
%      
%                         p(i, vars(i) + 1) = (1-b)*t;
%                     end
%                 end
%                 
%                 r = rand;
%                 %realize random bernoulli variable
%                 %remember that p(i, 1) is actually the probability of j=0
%                 if r <= p(i, 1)
%                     vars(i) = 0;
%                 else
%                     vars(i) = 1;
%                 end                   
% 
%                 if rem(iter_counter, 10000) == 0
%                     fprintf('Iter: %d\n', iter_counter);
%                     fprintf('%.10f\n', (cputime - last_time)*1000)
%                     save_vars(vars);
%                     last_time = cputime;
%                 end
%         
%                 iter_counter = iter_counter + 1;
%             end
%         end
% 
%         %check if we are completely done
%         for r = 1:N
%             if vars(r) == -1
%                 continue
%             end
%             
%             for c = 1:D
%                 if p(r, c) ~= 0 && p(r, c) ~= 1
%                     done = 0;
%                 end
%             end
%         end
% 
%         if done
%             p_size = size(p)
% 
%             p_cell_mat = [p_cell_mat; p];
%             save('p_cell_mat.mat', 'p_cell_mat');
%            break 
%         end
% 
% 
%     end

    %fprintf('ANSWER:: \n')
    %disp_vars(vars, V, P, T, E, edge)
end

%get of input nodes for terminal edges given a terminal and a flow
function res = get_terminal_inputs(vars, s, t, V, P, T, E, IV, SS, TT)
    %input nodes to terminal
    iv = IV{TT(t)};
    res = [];
    for i = 1:length(iv)
        f_val = get_f_from_vars(vars, iv(i), TT(t), s, t, V, P, T);
        
        %if f_val is -1, then by the dimensionality reduction, the edge
        %from iv(i) to t does not fall on any path from s to t, meaning it
        %is not a terminal input node
        if f_val ~= -1
            res = [res, iv(i)]; 
        end 
    end
end

function init_variable_groups(V, P, T, ST)
    global D_variable_group
    global N_variable_group
    global p_variable_group
    global variable_group_vals
    
    max_D = 0;
    D_variable_group = zeros(V, V, T);
    for i = 1:V
        for j = 1:V
            for t = 1:T
                D = sum(ST(:, t)) + 1;
                D_variable_group(i, j, t) = D;
                N_variable_group = N_variable_group + 1;
                
                if D > max_D
                    max_D = D;
                end
            end
        end
    end
    
    p_variable_group = zeros(V, V, T, max_D);
    for i = 1:V
        for j = 1:V
            for t = 1:T
                D = D_variable_group(i, j, t);
                for d = 1:D
                    p_variable_group(i, j, t, d) = 1/D;
                end
            end
        end
    end
    
    variable_group_vals = zeros(V, V, T);
    
end

%cost constraint - changes w/ every satisfying solution found
function res=cost_clause(vars, V, P, T, E, edge, ST, OV)
    global min_cost
    global z

    %used to the set the z matrix correctly
    for v = 1:V
        flow_limit(v, vars, V, P, T, E, edge, ST, OV);
    end

    res = 0;

    z = compute_z(vars, V, P, T, E, edge, ST, OV);
    [f, ~] = vars_to_mat(vars, V, P, T);
    cost = sum(sum(z(z ~= -1)));
    %f
    %z
    %fprintf('Cost: %d, min_cost: %d\n', cost, min_cost);

    if cost < min_cost
        res = 1; 
    end   
end
 
function [f, x]=vars_to_mat(vars, V, P, T)
    f = zeros(V, V, P, T);
    x = zeros(V, V, P);
    
    for i = 1:V
        for j = 1:V
            for s = 1:P
                for t = 1:T
                    f(i, j, s, t) = get_f_from_vars(vars, i, j, s, t, V, P, T);
                end
            end
        end
    end
    
    for i = 1:V
        for j = 1:V
            for s = 1:P
                x(i, j, s) = get_x_from_vars(vars, i, j, s, V, P, T);
            end
        end
    end
    
end

%generates the clause-variable participation matrix
function clause_mat=generate_clause_mat(IV, OV, V, P, T, ST)
    nclauses = V+V+V+V+1+1;
    nvars = V*V*P*T + V*V*P;
    clause_mat = zeros(nvars, nclauses);

    for v = 1:V
        iv = IV{v};
        ov = OV{v};
        for s = 1:P
            for i = 1:length(ov)
                if ov(i) ~= 0
                    clause_mat(linear_index_x(v, ov(i), s, V, P, T), V+v) = 1;
                end
            end
            
            for t = 1:T
                if ST(s, t) == 1
                    %constraint (28)
                    for i = 1:length(iv)
                        if iv(i) ~= 0
                            clause_mat(linear_index_f(iv(i), v, s, t, V, P, T), v) = 1;
                        end
                    end
                    for i = 1:length(ov)
                        if ov(i) ~= 0
                            clause_mat(linear_index_f(v, ov(i), s, t, V, P, T), v) = 1; 
                        end
                    end
                    
                    %constraint (29)
                    for i = 1:length(ov)
                        if ov(i) ~= 0
                            clause_mat(linear_index_f(v, ov(i), s, t, V, P, T), V+v) = 1;
                        end
                    end
                    
                    %constraint (32)
                    ivt = IV{t};
                    if length(ivt(ivt == v)) ~= 0
                        clause_mat(linear_index_x(v, t, s, V, P, T), 2*V+v) = 1;
                    end
                    
                    %constraint (17)
                    for i = 1:length(ov)
                        clause_mat(linear_index_f(v, ov(i), s, t, V, P, T), 3*V+v) = 1;
                    end
                end 
            end
        end      
    end
end

function bitstr=vars_to_bitstr(vars)
    bitstr = '';
    for v = 1:length(vars)
        if vars(v) >= 0
            bitstr = strcat(bitstr, num2str(vars(v)));
        end
    end
end

function save_vars(vars)
    fileId = fopen('vars.out', 'a+');

    bitstr = vars_to_bitstr(vars);
    bitstr = strcat(bitstr, '\n');

    fprintf(fileId, bitstr); 
end

function [vars, cont]=explore_vars_x_iter(vars, v, s, EE, SS, V, P, T, E, edge_imag)
	cont = 0;
	edges = EE(v, :);

    if sum(edges) == 0
        if v == SS(s)
            %fprintf('Found source (x). v: %d, s: %d\n', v, s)
            cont = 1;
            return
        end
    end
        
    for e = 1:length(edges)
        if edges(e) == 1
            [vars,cont_x] = explore_vars_x_iter(vars, e, s, EE, SS, V, P, T, E, edge_imag);
            %this means we reached the terminal node when explored
            if cont_x
                %fprintf('Travelling down (x). i: %d, j: %d, s: %d\n', e, v, s)
                vars = set_x_in_vars(vars, 0, e, v, s, V, P, T);
                %should set f to 0
                cont = cont_x;
            end
        end
    end
end


function vars=explore_vars_x(vars, EE, SS, V, P, T, E, TT, ST, edge_imag)
        for ti = 1:length(TT)
        for si = 1:length(SS)
            [vars, cont] = explore_vars_x_iter(vars, TT(ti), si, transpose(EE), SS, V, P, T, E, edge_imag);
        end
    end
end

function [vars, cont]=explore_vars_f_iter(vars, v, s, t, EE, SS, V, P, T, E, edge_imag)
    cont = 0;
    edges = EE(v, :);

    if sum(edges) == 0
        if v == SS(s)
            %fprintf('Found source. v: %d, s: %d, t:%d\n', v, s, t)
            cont = 1;
            return
        end
    end

    for e = 1:length(edges)
        if edges(e) == 1
            [vars,cont_x] = explore_vars_f_iter(vars, e, s, t, EE, SS, V, P, T, E, edge_imag);
            %this means we reached the terminal node when explored
            if cont_x
                %fprintf('Travelling down. i: %d, j: %d, s: %d, t: %d\n', v, e, s, t)
                vars = set_f_in_vars(vars, 0, e, v, s, t, V, P, T);
                %should set f to 0
                cont = cont_x;
            end
        end
    end
end

function vars=explore_vars_f(vars, EE, SS, V, P, T, E, TT, ST, edge_imag)
        
    tt_length = length(TT);
    %fprintf('Exploring f variables...\n')
    
    for ti = 1:length(TT)
        s_vec = ST(:, ti);
        for si = 1:length(s_vec);
            if s_vec(si) == 1
                t = TT(ti);
                [vars, cont] = explore_vars_f_iter(vars, t, si, ti, transpose(EE), SS, V, P, T, E, edge_imag);
                vars;
            end
        end
    end
    
end

function stack=explore(v, s, EE, SS, ~)
    edges = EE(v, :);
    stack = [v];

    if sum(edges) == 0
        if v == SS(s)
            return
        else
            %set all the right things to -1
            return
        end
    end

    for e = 1:length(edges)
        %if an edge exists, we follow it
        if edges(e) == 1
            explored = explore(e, s, EE, SS);
            stack = [stack, explored];
        end
    end
end

%sets possible f's to zero in vars
function vars=explore_paths_f(vars, V, P, T, EE, SS, TT, ST)
    mat = zeros(V, V);


    for ti = 1:length(TT)
        s_vec = ST(:, ti);
        for si = 1:length(s_vec)
            if s_vec(si) == 1
                t = TT(ti);
                stack = explore(t, si, transpose(EE), SS, V);
                %mark the right f values as 0's
                for i = length(stack):-1:2
                    vars = set_f_in_vars(vars,0,stack(i),stack(i-1),si,ti,V,P,T); 
                end
            end
        end
    end
end

function vars=explore_paths_x(vars, V, P, T, EE, SS, TT, ST)
    for ti = 1:length(TT)
        for si = 1:length(SS)
            t = TT(ti);
            stack = explore(t, si, transpose(EE), SS, V);
            
            %mark the right x values as 0's
            for i = length(stack):-1:2
                vars = set_x_in_vars(vars, 0, stack(i), stack(i-1),si,V,P,T);
            end
        end
    end   
end

%Constraint (17)
function checkf=flow_limit(v, vars, V, P, T, E, edge, ST, OV)
    checkf = 1;

    global z;
    zt = zeros(V, V);

    ov = OV{v};

    for t=1:T
        sumf=0;
        for s=1:P
            for i = 1:length(ov)
                fval = get_f_from_vars(vars, v,ov(i),s,t, V,P,T);
                if ST(s,t)==1 && fval ~= -1
                    sumf=sumf + fval;
                end
            end
        end

        if sumf>1
            checkf=0;
        else
            zt(v,t)=sumf;
        end
    end

    if sum(ov) ~= 0
        zt_row = zt(v, :);
        for i = 1:length(ov)
            mq = [zt(v, :), z(v, ov(i))];
            mq;
            z(v,ov(i))=max(mq);
        end
    end

end

function z=compute_z(vars, V, P, T, E, edge, ST, OV)
    z = zeros(V, V);
    zt = zeros(E, T);

    for e = 1:E
        iv=real(edge(e));
        ov=imag(edge(e));
        for t=1:T
            sumf=0;
            for s=1:P
                if ST(s,t)==1
                    val = get_f_from_vars(vars, iv,ov,s,t, V, P, T);
                    if val ~= -1
                       sumf=sumf+val;
                    end
                end
            end
            zt(e,t)=sumf;
            if sumf>1
                checkf=0;
            end
        end
        z(iv,ov)=max(zt(e,:));
    end
end

%clause2: Constraint (28)
function checkfv=flow_conservation(v, vars, V, P, T, OV, IV, sigma, ST, E, edge)
    checkfv = 1;
    for t=1:T
        for s=1:P
            if ST(s,t)==1
                %sum of outgoing flows from v
                sumoutf=0;
                if sum(OV{v}~=0)>=1
                    for ovidx=1:length(OV{v})
                        ov=OV{v}(ovidx);
                        f_val = get_f_from_vars(vars,v,ov,s,t,V,P,T);
                        if f_val ~= -1
                            sumoutf=sumoutf+f_val;
                        end
                    end
                end

                suminf=0;
                if sum(IV{v}~=0)>=1
                    for ividx=1:length(IV{v})
                        iv=IV{v}(ividx);
                        f_val = get_f_from_vars(vars,iv,v,s,t,V,P,T);
                        if f_val ~= -1
                            suminf=suminf+get_f_from_vars(vars,iv,v,s,t,V,P,T);
                        end
                    end
                end
                
                if (sumoutf-suminf~=sigma(v,s,t))
                    %disp_vars(vars, V, P, T, E, edge);
                    %fprintf('flow_conservation, v: %d\n', v); 
                    checkfv=0;
                end
            end
        end
    end
end

%clause 3: Constraint (29)
function checkfx=flow_x(v, vars, OV, V, P, T, E, ST, edge, p)
    checkfx=1;
    for s=1:P
        for t=1:T
            if ST(s,t) == 1
                ov = OV{v};
                for j = 1:length(ov)
                    fval = get_f_from_vars(vars,v,ov(j),s,t,V,P,T);
                    xval = get_x_from_vars(vars,v,ov(j),s,V,P,T);
                    
                    if fval > xval && fval ~= -1 && xval ~= -1
                        %disp_vars(vars, V, P, T, E, edge);
                        %fprintf('flow_x, fval: %d, xval: %d, v: %d, ov(j): %d, s: %d, t: %d\n', fval, xval, v, ov(j), s, t);
                        checkfx=0;
                    end
                end
            end
        end
    end
end


%clause 4: Constraint (32)
function res=checkx(v, vars, V, P, T, ST, IV, TT)
    res = 1;
    for t = 1:T
         for s=1:P
             if (ST(s,t)~=1)
                 xval = get_x_from_vars(vars, v, TT(t), s, V, P, T);
                 if xval ~=0 && xval ~= -1
                     res=0;
                 end
             end
         end 
     end
end

%clause 5: Constraint (33)
%this involves trying multiple beta values
function res=checkbeta(vars, V, P, T, IV, E, edges)
    res = 1;
    global beta_mat

    for e = 1:E
        i = real(edges(e));
        j = imag(edges(e));
        
        max_vec = zeros(P);
        
        for k = 1:length(IV{i})
            for beta = 0:1
                for p = 1:P
                    xval = get_x_from_vars(vars,i, j, p, V, P, T);
                    if xval ~= -1
                        vec_p = beta * get_x_from_vars(vars,i, j, p, V, P, T);
                        beta_mat(k, i, j) = beta;
                        if vec_p > max_vec(p)
                            max_vec(p) = vec_p;
                        end
                    end
                end
            end
        end
        
        for p = 1:P
            xval = get_x_from_vars(vars,i, j, p, V, P, T);
            if xval ~= -1
                if max_vec(p) ~= get_x_from_vars(vars,i, j, p, V, P, T)
                    res = 0;
                end
            end
        end
    end
end


%takes vars and prints it
function disp_vars_p(vars, V, P, T, E, edge)
    res_str = '';
    
    for e = 1:E
        i = real(edge(e));
        j = imag(edge(e));
        
        for p = 1:P
            for t = 1:T
                val = get_f_from_vars(vars, i, j, p, t, V, P, T);
                if val >= 0
                    res_str = strcat(res_str, int2str(val));
                end
            end
        end
    end
    
    for e = 1:E
        i = real(edge(e));
        j = imag(edge(e));
        for p = 1:P
            val = get_x_from_vars(vars, i, j, p, V, P, T);
            if val >= 0
                res_str = strcat(res_str, int2str(val));
            end
        end
    end
    
    disp(res_str)
    
end

%takes vars and prints it in a human-readable form
function disp_vars(vars, V, P, T, E, edge)
    for p = 1:P
        for t = 1:T
            for e = 1:E
                i = real(edge(e));
                j = imag(edge(e));
                val = get_f_from_vars(vars, i, j, p, t, V, P, T);
                if val >= 0
                    fprintf('f; i %d, j: %d, p: %d, t: %d, val: %d\n', i, j, p, t, val);
                end
            end
        end
    end
    
%     for e = 1:E
%         i = real(edge(e));
%         j = imag(edge(e));
%         
%         for p = 1:P
%             for t = 1:T
%                 val = get_f_from_vars(vars, i, j, p, t, V, P, T);
%                 if val >= 0
%                     fprintf('f; i: %d, j: %d, p: %d, t: %d, val: %d\n', i, j, p, t, val)
%                 end
%             end
%         end
%     end
    
    for e = 1:E
        i = real(edge(e));
        j = imag(edge(e));
        for p = 1:P
            val = get_x_from_vars(vars, i, j, p, V, P, T);
            if val >= 0
                fprintf('x; i: %d, j: %d, p: %d, val: %d\n', i, j, p, val);
            end
        end
    end
end


%fix all of the flow x's as 1's
function [res, prob, clause_mat]=set_source_links(vars, prob, clause_mat, V, P, T, EE, SS)
    %loop over the sources, check the edges from each source
    %set those as one in vars
    
    for p1 = 1:P
        s_p = SS(p1);
        for p2 = 1:P
            for j = 1:V
                if EE(s_p, j) == 1 && get_x_from_vars(vars, s_p, j, p2, V, P, T) ~= -1
                    vars = set_x_in_vars(vars, 0, s_p, j, p2, V, P, T);
                    
                    %we also set the initial probability values
                    index = linear_index_x(s_p, j, p2, V, P, T);
                    prob(index, 1) = 1;
                    prob(index, 2) = 0;
                    
                    %these are set in stone so they are not part of the
                    %clauses
                    clause_mat(index, 4) = 0;
                    clause_mat(index, 5) = 0;
                end
            end
        end
        
        for j = 1:V
            if EE(s_p, j) == 1 && get_x_from_vars(vars, s_p, j, p1, V, P, T) ~= -1
                vars = set_x_in_vars(vars, 1, s_p, j, p1, V, P, T);
                prob(linear_index_x(s_p, j, p2, V, P, T), 1) = 0;
                prob(linear_index_x(s_p, j, p2, V, P, T), 2) = 1;
            end
        end
    end
    res = vars;
end

%sets vars so that nonexistent edges are -1
function res=remove_nonexistent_edges(vars, V, P, T, EE)
    for i = 1:V
        for j = 1:V
            %if there is no edge at a given place,
            %then we can set the values of f and x
            %that correspond to it as -1
            if EE(i, j) == 0
                for p = 1:P
                    vars = set_x_in_vars(vars, -1, i, j, p, V, P, T);
                    for t = 1:T
                        vars = set_f_in_vars(vars, -1, i, j, p, t, V,P,T);
                    end
                end
            end
        end
    end
    res = vars;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Path detection                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cont, cell_mat, path_count, paths]=get_paths_explore_iter(s, t, EE, cell_mat, path_count, paths)
    edges = EE(s, :);
    paths = [paths, s];

    %edge exists between the starting point and terminal;
    %that means we are done with this iteration of DFS
    if s == t
        cont = 1;
        cell_mat{path_count} = paths;
        path_count = path_count + 1;
        return
    end

    cont = 0;
    for e = 1:length(edges)
        if edges(e) == 1 
            [cont_x, cell_mat, path_count, ~] = get_paths_explore_iter(e, t, EE, cell_mat, path_count, paths);
            if cont_x == 1
                cont =  1;
            end 
        end
    end
end

%Note: t represents the actual terminal value, not the index of TT

function final_mat=get_paths(t, S, SS, EE, V)
    final_mat = cell(1, S);
    for i = 1:S
        cell_mat = cell(i, V);
        [cont, cell_mat, path_count, paths] = get_paths_explore_iter(...
            SS(i), t, EE, cell_mat, 1, []);
        final_mat{i} = cell_mat;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear indexing                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=linear_index_f(i, j, p, t, V, P, T)
    res = (i-1)*V*P*T + (j-1)*P*T + (p-1)*T + t;
end

function res=linear_index_x(i, j, p, V, P, T)
    res = V*V*P*T + (i-1)*V*P + (j-1)*P + p;
end

function res=set_f_in_vars(vars, val, i, j, p, t, V, P, T)
    vars((i-1)*V*P*T + (j-1)*P*T + (p-1)*T + t) = val;
    res = vars;
end

function res = set_x_in_vars(vars, val, i, j, p, V, P, T)
    vars(V*V*P*T + (i-1)*V*P + (j-1)*P + p) = val;
    res = vars;
end

%get value of f from the long "vars" vector
function res=get_f_from_vars(vars, i, j, p, t, V, P, T)
    res = vars((i-1)*V*P*T + (j-1)*P*T + (p-1)*T + t);
end

%get value of x from the long "vars" vector
function res=get_x_from_vars(vars, i, j, p, V, P, T)
    res = vars(V*V*P*T + (i-1)*V*P + (j-1)*P + p);
end
