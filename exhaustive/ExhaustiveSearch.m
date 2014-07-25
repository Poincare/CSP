%demand_set_expansion: boolean variable that decides whether or the demand set is to be expanded
%shortest_path_depth: how far we had to go down the ranked paths list in order to find a feasible solution
function [min_cost_z, shortest_path_depth, path_count] = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST,...
        demand_set_expansion, ROUTING, ATOMS, clause_mat_l)
    %function z=ExhaustiveSearch()
    
    global m
    global fea_idx
    global fea_x
    global fea_beta
    global fea_z
    global fea_f
    global fea_usum;
    
    global clause_mat
    clause_mat = clause_mat_l;

    shortest_path_depth = -1;
    
    m=1;
    
    fea_idx = 1;
    fea_z = zeros(V, V, V);
    P = S;
    
    path_count = -1;
    
    %expand the demand set if we are told to do
    %so
    if demand_set_expansion == 1
        ST_possib = expand_demand_set(S, T, ST);
    else
        ST_possib = ST;
    end

    min_cost = Inf;
    min_cost_z = zeros(V, V);
    min_cost_st = zeros(S, T);
    
    %input node
    IV=cell(V);
    for v=1:V
        IV{v}=0;
        for iv=1:V
            if EE(iv,v)==1
                if sum(IV{v})==0
                    IV{v}=iv;
                else
                    IV{v}=[IV{v},iv];
                end
            end
        end
    end
    
    %output node
    OV=cell(V);
    for v=1:V
        OV{v}=0;
        for ov=1:V
            if EE(v,ov)==1
                if sum(OV{v})==0
                    OV{v}=ov;
                else
                    OV{v}=[OV{v},ov];
                end
            end
        end
    end
    
    %topologically sort the nodes (used for
    %computing values of x from beta and f)
    global sorted_nodes
    sorted_nodes = TopologicalSort(EE, IV);

    %edge vector
    E = sum(sum(EE));
    edge=zeros(1,E);
    e=1;
    for l=1:V
        for m=1:V
            if EE(l,m)==1
                edge(e)=l+m*1i;
                e=e+1;
            end
        end
    end
    
    for st_i = 1:length(ST_possib)
        %select one ST from the expanded demand set
        if iscell(ST_possib) == 1
            ST = ST_possib{st_i};
        else
            ST = ST_possib;
        end
        
        % sigma
        sigma=zeros(V,S,T);
        for v=1:V
            for s=1:S
                for t=1:T
                    if ST(s,t)==1
                        if v==SS(s)
                            sigma(v,s,t)=1;
                        elseif v==TT(t)
                            sigma(v,s,t)=-1;
                        end
                    end
                end
            end
        end
        
        %st_imag vector
        st_size = sum(sum(ST));
        st_imag = zeros(1, st_size);
        st_index = 1;
        for si = 1:S
            for ti = 1:T
                if ST(si, ti) == 1
                    st_imag(st_index) = si + (ti * 1i);
                    st_index = 1 + st_index;
                end
            end
        end
        
        %initialize variable matrices
        z=-ones(V,V);
        x=-ones(V,V,S);
        beta=-ones(V,V,V);
        f=-ones(V,V,S,T);
        for l=1:V %column of EE matrix
            for m=1:V %row of EE matrix
                if EE(l,m)==1 %if l,m work in the EE matrix(when there is a 1 instead of 0)
                    z(l,m)=0; %show 0 instead of -1
                    for s=1:S
                        x(l,m,s)=0;
                        %for t=1:T
                        %if ST(s,t)==1
                        %f(l,m,s,t)=0;
                        %end
                        %end
                    end
                    for n=1:V
                        if EE(m,n)==1
                            beta(l,m,n)=0;
                        end
                    end
                end
            end
        end
        
        %use Depth First Search to set the correct variable spots for f.
        f = explore_vars_f(f, EE, SS, V, P, T, E, TT, ST, edge);
        
        if ATOMS
            global C
            compute_atom_globals(V, OV, ST, f, S, T);
        end
        
        global path_set_costs
        global path_set_idx
        
        % (path_stack) : (path_cost)
        path_set_costs = cell(1, 2);
        path_set_idx = 1;
        shortest_path_rankings(1, st_imag, cell(0, 0), V, T, TT, S, SS, EE);

        [~, non_empty_size] = size(path_set_costs(~cellfun('isempty', path_set_costs)));
        if non_empty_size == 0
            z = zeros(V, V);
            return
        end
        
        path_set_costs = sortrows(path_set_costs, 2);
        
        [path_count, ~] = size(path_set_costs);
        
        % for k = 1:path_count
        %     path_set = path_set_costs{k, 1};
        %     for l = 1:length(path_set)
        %         path_set{l}
        %     end
        %     fprintf('----\n');
        % end

        global cost_mat
        
        %path_count
        for k = 1:path_count
            fs = f;
            betas = beta;

            path_set = path_set_costs{k, 1};
            %for l = 1:length(path_set)
            %    path_set{l}
            %end


            [fs, betas] = gen_vars_from_path_set(path_set, V, P, T, SS, TT, fs, betas);
            res = check_feasibility(fs, betas, V, S, T, SS, TT, OV, IV, E, EE, edge, ST, sigma, ROUTING, ATOMS);

            if res
                %first feasible solution is the optimal solution
                z = compute_z(fs, edge, E, V, S, T);
                %EE
                %ST
                %min_cost_z
                z = z(1:V-S-T, 1:V-S-T);
                shortest_path_depth = path_count;
                
                cost = GetCost(z, cost_mat);
                if cost < min_cost && (cost > 0)
                   min_cost = cost;
                   min_cost_z = z;
                   min_cost_st = ST;
                end
            else
                z = zeros(V, V);
            end
        end
    end

    return
end

function atom_mat=d_atoms(i, atoms, ti, atom_mat)
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
    global D
 
    %result matrices
    global z x beta f
    
    terminals = TP{i};
 
    %base case
    if ti > length(terminals)
        atom_mat = cell(1, 1);
        atom_mat{1} = atoms;
        return
    end
   
    % if i == 8
    %     terminals_ti = terminals(ti) 
    %     atoms
    %     to_intersect = D{terminals(ti)}
    % end

    atoms_l = intersect(atoms, D{terminals(ti)});
    % if i == 8
    %     atoms_l
    % end

    if length(atoms_l) ~= 0
        atom_mat_l = d_atoms(i, atoms_l, ti + 1, atom_mat);
        if ~isempty(atom_mat_l)
            atom_mat = [atom_mat; atom_mat_l];
        end
    end
    
    whole = 1:S;
    d_comp = setdiff(whole, D{terminals(ti)});
    atoms_r = intersect(atoms, d_comp);
    if length(atoms_r) ~= 0
        atom_mat_r = d_atoms(i, atoms_r, ti + 1, atom_mat);
        if ~isempty(atom_mat_r)
            atom_mat = [atom_mat; atom_mat_r];
        end
    end
end


%expand the demand set (can give feasible solution even when the original demand set
%does not give a feasible solution).
function ST_possib = expand_demand_set(S, T, ST)
    st_zero_width = length(ST(ST == 0));
    
    ST_possib = cell(0, 1);
    
    %no freedom in picking ST supersets
    if st_zero_width == 0
        ST_possib = ST
        return
    end
    
    st_zero_pos = find(ST == 0);
    comb = 2^(st_zero_width) - 1;
    
    for tcomb = 0:comb
        STs = ST;
        
        bin_vec = de2bi(tcomb);
        
        st_i = 1;
        for i = 1:length(bin_vec)
            STs(st_zero_pos(st_i)) = bin_vec(i);
            st_i = st_i + 1;
        end
        
        ST_possib = [ST_possib, STs];
    end
end

%get cost for a path
%REWRITE
function cost=path_stack_cost(path_stack)
    global cost_mat

    cost = 0;
    edges = zeros(1, 2);
    edge_idx = 1;
    
    for i = 1:length(path_stack)
        path = path_stack{i};
        for i = 1:length(path)-1
            m = path(i);
            n = path(i+1);
            
            edges(edge_idx, 1) = m;
            edges(edge_idx, 2) = n;
            edge_idx = edge_idx + 1;
        end
    end
    
    %note: this has to be changed if a different
    %cost function is to be used
    cost_edges = unique(edges, 'rows');
    [num_cost_edges, ~] = size(cost_edges);

    for ce = 1:num_cost_edges
        m = cost_edges(ce, 1);
        n = cost_edges(ce, 2);

        cost = cost + cost_mat(m, n);
    end
end

function shortest_path_rankings(st_index, st_imag, path_stack, V, T, TT, S, SS, EE)
    global path_set_costs
    global path_set_idx

    if st_index > length(st_imag)
        %fprintf('hit.\n');
        cost = path_stack_cost(path_stack);
        path_set_costs{path_set_idx, 1} = path_stack;
        path_set_costs{path_set_idx, 2} = cost;
        path_set_idx = path_set_idx + 1;
        
        return
    end
    
    si = real(st_imag(st_index));
    ti = imag(st_imag(st_index));
    
    paths_c = get_paths(TT(ti), S, SS, EE, V);
    paths = paths_c{si};
    source = SS(si);
    terminal = TT(ti);
    
    path_length = length(paths);
    
    for path_i = 1:length(paths)
        path = paths(path_i);
        if ~isempty(path{1})
            %path{1}
            shortest_path_rankings(st_index + 1, st_imag, [path_stack; path], V, T, TT, S, SS, EE);
        end
    end
    
end

function [f, beta] = gen_vars_from_path_set(path_set, V, P, T, SS, TT, f, beta)
    beta = zeros(V, V, V);
    %f = zeros(V, V, P, T);
    
    for p = 1:length(path_set)
        path_m = path_set{p};

        source = find(SS == path_m(1));
        si = source(1);
        terminal = find(TT == path_m(length(path_m)));
        ti = terminal(1);
        
        %set the right f values for all of the edges
        %along this path
        for i = 1:length(path_m)-1
            m = path_m(i);
            n = path_m(i+1);
            
            f(m, n, si, ti) = 1;
        end
        
        %set beta values for all of the input/output edge
        %pairs along this beta
        for i = 1:length(path_m)-2
            m = path_m(i);
            n = path_m(i + 1);
            o = path_m(i + 2);
            
            beta(m, n, o) = 1;
        end
    end
    %fprintf('----\n')
end

function compute_atom_globals(V, OV, ST, f, S, T)
    global C TP D
    
    TP = cell(V, 1);
    for v = 1:V
        ov = OV{v};
        terminals = [];
        
        for t = 1:T
            for s = 1:S
                if ST(s, t) == 1
                    for j = 1:length(ov)
                        if ov(j) ~= 0 && f(v, ov(j), s, t) ~= -1
                            terminals = unique([terminals, t]);
                        end
                    end
                end
            end
        end
        
        TP{v} = terminals;
    end
    
    %set up D(t) from 2004isit paper
    D = cell(T, 1);
    for t = 1:T
        sources = [];
        
        for s = 1:S
            if ST(s, t) == 1
                sources = [sources, s];
            end
        end
        
        D{t} = sources;
    end

    C = cell(1, V);
    for i = 1:V
        ti_size = length(TP{i});
        atom_mat = d_atoms(i, 1:S, 1, cell(1, 1));
 
        atoms = cell(1, S);
        [rows, cols] = size(atom_mat);
        atoms = cell(1, 2^ti_size);
        j = 1;
        for r = 1:rows
            if length(atom_mat{r}) ~= 0
                atoms{j} = atom_mat{r};
                j = j + 1;
            end
        end
        
        C{i} = atoms;
    end
 
    save('c.mat', 'C', 'TP', 'D'); 
end

function z=compute_z(f, edge, E, V, S, T)
    z = zeros(V, V);
    for e = 1:E
        iv = real(edge(e));
        ov = imag(edge(e));
        
        for s = 1:S
            for t = 1:T
                if f(iv, ov, s, t) == 1
                    z(iv, ov) = 1;
                end
            end
        end
    end
end

function res = check_feasibility(f, beta, V, S, T, SS, TT, OV, IV, E, EE, edge, ST, sigma, ROUTING, ATOMS)
    global m
    
    function res = sum_cell(cell_m)
        res = cellfun(@sum, cell_m);
    end
            
    function atom_s = check_atom_feasibility(x, f)
        global C D TP

        atom_s = 1;
        
        for e = 1:E
            i = real(edge(e));
            j = imag(edge(e));
            C_j = cell(1, 0);
            

            for c = 1:length(C{j})
                found = 0;
                for cji = 1:length(C_j)
                    if isequal(C_j{cji}, C{j}{c})
                        found = 1;
                        break;
                    end
                end

                if ~found
                    C_j{c} = C{j}{c};
                end
            end

            Value_sets = cell(1, length(C_j));

            for c = 1:length(C_j)
                path_set = C_j{c};

                if ~isempty(path_set)
                    value_set = zeros(1, length(path_set));
                    for k = 1:length(path_set)
                        if x(i, j, path_set(k)) ~= -1
                            value_set(k) = x(i, j, path_set(k));
                        end
                    end

                    Value_sets{c} = {value_set};
                end
            end

            
            value_sums = [];
            for vs = 1:length(Value_sets)
                value_set = Value_sets{vs};
                
                if iscell(value_set)
                    value_set = value_set{1};
                end
                
                val_sum = 0;
                
                for m = 1:length(value_set)
                    val_sum = val_sum + value_set(m);
                end
                
                value_sums = [value_sums, val_sum];
            end
            
            
            TOTAL_VALUE = length(value_sums(value_sums > 0));
            if TOTAL_VALUE > 1
                atom_s = 0;
            end
                 
            % if (j == 8) && (atom_s == 1)
            %     C{8}

            %     Value_sets

            %     for k = 1:length(Value_sets)
            %         value_sets_k = Value_sets{k}
            %     end
            %     value_sums
            %     TOTAL_VALUE
            % end
            
        end
    end
    
    
    %check conditions...
    %f = zeros(V, V, S, T);
    x = -ones(V, V, S);
    m = m + 1;
    
    x = set_x_on_initial_edges(x, f, OV, IV, V, S, T, SS, TT);
    x = compute_x(x, beta, EE, IV, OV, V, S, T);

    %% checking constraints 1:satisfied
    checkx=1;
    checkfx=1;
    checkf=1;
    checkfv=1; %Fow Conservation
    checkrouting = 1;
    checkatoms = 1;
    

    if ROUTING == 1
        e = 1;
        while (e <= E)
            iv = real(edge(e));
            ov = imag(edge(e));
            
            if ~isempty(find(TT == ov))
               e = e+1;
                continue;   
            end
            
            x_sum = 0;
            for s = 1:S
                val = x(iv, ov, s);
                if val ~= -1
                    x_sum = x_sum + val;
                end
            end
            
            if x_sum > 1
                %fprintf('Failed routing; iv: %d, ov: %d\n', iv, ov);
                checkrouting = 0;
            end
            
            e = e+1;
        end
    end
    
    %if we are doing ATOMS, then we need to check:
    %(1) atoms constraints
    %(2) checkf constraint (outside of if)
    %(3) flow conservation (outside of if)
    if ATOMS == 1
        checkatoms = check_atom_feasibility(x, f);
    else
        %check x at terminals
        t=1;
        while (t<=T)&&(checkx==1)
            for s=1:S
                if (ST(s,t)~=1)
                    for iv=1:length(IV{TT(t)})
                        if x(IV{TT(t)}(iv),TT(t),s)~=0
                            %fprintf('Checkx failed at iv: %d, t: %d, s: %d\n, ATOMS: %d\n, ROUTING: %d\n', IV{TT(t)}(iv), TT(t), s, ATOMS, ROUTING);
                   
                            checkx=0;
                        end
                    end
                end
            end
            t=t+1;
        end
        
        %check f<=x %Constraint (29)
        e=1;
        while (checkx==1)&&(e<=E)&&(checkfx==1)
            iv=real(edge(e)); %input node (represented by real number)
            ov=imag(edge(e)); %output node (represented by imaginary number)
            for s=1:S
                for t=1:T
                    if ST(s,t)==1
                        if f(iv,ov,s,t)>x(iv,ov,s)
                            %fprintf('Failed checkfx for iv: %d, ov:  %d, s: %d, t: %d\n', iv, ov, s, t);
                            %fprintf('fval: %d, xval: %d\n', f(iv, ov, s, t), x(iv, ov, s));
                            %f(:, :, s, t)
                            checkfx=0;
                        end
                    end
                end
            end
            e=e+1;
        end
        
    end
    
    %check f %Consraint (17)
    e=1;
    zt=zeros(E,T);
    z = zeros(V, V);
    while (checkx==1)&&(checkfx==1)&&(e<=E)&&(checkf==1)
        iv=real(edge(e));
        ov=imag(edge(e));

        for t=1:T
            if ov == TT(t)
                continue;
            end
            
            sumf=0;
            for s=1:S
                if ST(s,t)==1
                    val = f(iv,ov,s,t);
                    if val ~= -1
                        sumf=sumf+val;
                    end
                end
            end
            zt(e,t)=sumf;
            if sumf>1

                %fprintf('Failed on: iv: %d, ov: %d, t: %d, sumf: %d\n', iv, ov, t, sumf);
                checkf=0;
            end
        end
        z(iv,ov)=max(zt(e,:));
        e=e+1;
    end
    
    % check f flow conservation
    %Constraint (28)
    v=1;
    while (checkx==1)&&(checkfx==1)&&(checkf==1)&&(v<=V)&&(checkfv==1)
        for t=1:T
            for s=1:S
                if ST(s,t)==1
                    %sum of outgoing flows from v
                    sumoutf=0;
                    if sum(OV{v}~=0)>=1
                        for ovidx=1:length(OV{v})
                            ov=OV{v}(ovidx);
                            val = f(v,ov,s,t);
                            %fprintf('v: %d, ov: %d, s: %d, t: %d, val: %d\n', v, ov, s, t, val);
                            if val ~= -1
                                sumoutf=sumoutf+val;
                            end
                        end
                    end
                    
                    suminf=0;
                    if sum(IV{v}~=0)>=1
                        for ividx=1:length(IV{v})
                            iv=IV{v}(ividx);
                            val = f(iv,v,s,t);
                            if val ~= -1
                                suminf=suminf+val;
                            end
                        end
                    end
                    
                    if (sumoutf-suminf~=sigma(v,s,t))
                        %fprintf('Failed on node: %d, t: %d, s: %d\n', v, t, s);
                        %fprintf('Sumoutf: %d, suminf: %d, sigma: %d\n', sumoutf, suminf, sigma(v, s, t));
                        checkfv=0;
                    end
                end
            end
        end
        v=v+1;
    end
    
    
    %  if ROUTING == 1
    %     checkx
    %     checkfx
    %     checkf
    %     checkfv
    %     checkroutin
    %     checkatoms
    %     fprintf('----------\n');
    % end
    
    res = checkx & checkfx & checkf & checkfv & checkrouting & checkatoms;

end

function disp_z(z, V)
    for i = 1:V
        for j= 1:V
            val = z(i, j);
            if val ~= 0
                fprintf('z, i: %d, j: %d, value: %d\n', i, j, val)
            end
        end
    end
end

function z = find_optimal(V)
    global fea_usum fea_z fea_beta fea_x fea_f
    
    r = sum(sum(sum(fea_z)));
    
    for idx = 1:length(fea_z(1, 1, :))
        if GetCost(fea_z(:, :, idx), cost_mat) ~= 0
            %fea_z(:, :, idx)
        else
            fprintf('No feasible solution found\n');
        end
    end
    
    min_cost = Inf;
    min_cost_z = zeros(V, V);
    
    for idx = 1:length(fea_z(1, 1, :))
        cost = GetCost(fea_z(:, :, idx), cost_mat);
        if (cost < min_cost) && (cost > 0)
            min_cost = cost;
            min_cost_z = fea_z(:, :, idx);
        end
    end
    
    z = min_cost_z;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path cutting                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%i: source index

function disp_f(f, S, T, V, beta)
    str = '';
    for s = 1:S
        for t=1:T
            if sum(sum(f(:, :, s, t))) ~= -V*V
                for i = 1:V
                    for j = 1:V
                        val = f(i, j, s, t);
                        if val ~= -1
                            fprintf('i: %d, j: %d, s: %d, t: %d, val: %d\n', i, j, s, t, val);
                            str = strcat(str, int2str(val));
                        end
                    end
                end
            end
        end
    end
    
    str
    beta_expected = -ones(V, V, V);
    
    if strcmp(str, '111110100000101000110010111010') == 1
        beta_expected(1,3,4)=1;
        beta_expected(3,4,6)=1;
        beta_expected(3,9,10)=1;
        beta_expected(4,6,7)=1;
        beta_expected(5,4,6)=1;
        beta_expected(4,6,10)=1;
        beta_expected(2,5,4)=1;
        beta_expected(2,5,7)=1;
        beta_expected(1,3,8)=1;
        beta_expected(1,3,9)=1;
        beta_expected(4,6,9)=0;
        beta_expected(6,9,10)=0;
        beta_expected(6,9,11)=0;
        beta_expected(9,11,8)=0;
        beta_expected(3, 9, 11) = 0;
        fprintf('answer\n');
    end
end

%sets the x values for source and terminal edges (1 and 0
%respectively)
function x=set_x_on_initial_edges(x, f, OV, IV, V, P, T, SS, TT)
    for s = 1:P
        ov = OV{SS(s)};
        for i = 1:length(ov)
            for t = 1:T
                %fprintf('SS(s): %d, ov(i): %d, s: %d, t: %d, f_val: %d\n', SS(s), ov(i), s, t, f(s, ov(i), s, t));
                if f(SS(s), ov(i), s, t) ~= -1
                    x(SS(s), ov(i), :) = 0;
                    x(SS(s), ov(i), s) = 1;
                    break;
                else
                    x(SS(s), ov(i), s) = 0;
                end
            end
        end
    end
    
    %for t=1:T
    %    iv = IV{TT(t)};
    %    for i = 1:length(iv)
    %        for s = 1:P
    %            fprintf('iv(i): %d, TT(t): %d, s: %d, t: %d, f_val: %d\n', iv(i), TT(t), s, t, f(iv(i), TT(t), s, t));
    
    %            if f(iv(i), TT(t), s, t) ~= -1
    %                x(iv(i), TT(t), s) = 0;
    %            else
    %                x(iv(i), TT(t), s) = -1;
    %            end
    %        end
    %    end
    %end
    
end

%computes x from beta values
%x, beta, EE, IV, OV, V, S, T
function x=compute_x(x, beta, EE, IV, OV, V, S, T)
    global sorted_nodes

    for n = 1:length(sorted_nodes)
        v = sorted_nodes(n);
        
        iv = IV{v};
        ov = OV{v};
        if sum(iv) ~= 0 && sum(ov) ~= 0
            for j = 1:length(ov)
                for s = 1:S
                    max_x = -Inf;
                    for k = 1:length(iv)
                        
                        %fprintf('beta(%d, %d, %d) * x(%d, %d, %d)\n', iv(k), v, ov(j), iv(k), v, s);
                        x_val = beta(iv(k), v, ov(j)) * x(iv(k), v, s);
                        if x_val > max_x
                            max_x = x_val;
                        end
                    end
                    
                    %fprintf('value set, s: %d\n', s);
                    x(v, ov(j), s) = max_x;
                end
            end
        end
    end
end

function path_comb(st_imag, st_index, f, beta, SS, TT, EE, ST, IV, OV, S, T, V, E, edge, sigma, z)
    %fprintf('\n--------------------------------\n')
    
    global m;
    global fea_x
    global fea_beta
    global fea_z
    global fea_f
    global fea_idx
    global fea_usum
    
    %if the source is the last one in consideration,
    %then we can check conditions
    if st_index > length(st_imag)
        %check conditions...
        x = -ones(V, V, S);
        m = m + 1;
        
        x = set_x_on_initial_edges(x, f, OV, IV, V, S, T, SS, TT);
        x = compute_x(x, beta, EE, IV, OV, V, S, T);
        
        %% checking constraints 1:satisfied
        checkx=1;
        checkfx=1;
        checkf=1;
        checkfv=1; %Flow Conservation
        
        %check x at terminals
        t=1;
        while (t<=T)&&(checkx==1)
            for s=1:S
                if (ST(s,t)~=1)
                    for iv=1:length(IV{TT(t)})
                        if x(IV{TT(t)}(iv),TT(t),s)~=0
                            checkx=0;
                        end
                    end
                end
            end
            t=t+1;
        end
        
        %check f<=x %Constraint (29)
        e=1;
        while (checkx==1)&&(e<=E)&&(checkfx==1)
            iv=real(edge(e)); %input node (represented by real number)
            ov=imag(edge(e)); %output node (represented by imaginary number)
            for s=1:S
                for t=1:T
                    if ST(s,t)==1
                        if f(iv,ov,s,t)>x(iv,ov,s)
                            checkfx=0;
                        end
                    end
                end
            end
            e=e+1;
        end
        
        
        %check f %Consraint (17)
        e=1;
        zt=zeros(E,T);
        z = zeros(V, V);
        while (checkx==1)&&(checkfx==1)&&(e<=E)&&(checkf==1)
            iv=real(edge(e));
            ov=imag(edge(e));

            %fprintf('Failed on: %d, %d\n', iv, ov)
            for t=1:T
                sumf=0;
                for s=1:S
                    if ST(s,t)==1
                        val = f(iv,ov,s,t);
                        if val ~= -1
                            sumf=sumf+val;
                        end
                    end
                end

                zt(e,t)=sumf;
                if sumf>1 && ov ~= TT(t)
                    checkf=0;
                end
            end
            z(iv,ov)=max(zt(e,:));
            e=e+1;
        end
        
        % check f flow conservation
        %Constraint (28)
        v=1;
        while (checkx==1)&&(checkfx==1)&&(checkf==1)&&(v<=V)&&(checkfv==1)
            for t=1:T
                for s=1:S
                    if ST(s,t)==1
                        %sum of outgoing flows from v
                        sumoutf=0;
                        if sum(OV{v}~=0)>=1
                            for ovidx=1:length(OV{v})
                                ov=OV{v}(ovidx);
                                val = f(v,ov,s,t);
                                if val ~= -1
                                    sumoutf=sumoutf+val;
                                end
                            end
                        end
                        
                        suminf=0;
                        if sum(IV{v}~=0)>=1
                            for ividx=1:length(IV{v})
                                iv=IV{v}(ividx);
                                val = f(iv,v,s,t);
                                if val ~= -1
                                    suminf=suminf+val;
                                end
                            end
                        end
                        
                        if (sumoutf-suminf~=sigma(v,s,t))
                            %fprintf('Failed on node: %d, t: %d, s: %d\n', v, t, s);
                            %fprintf('Sumoutf: %d, suminf: %d, sigma: %d\n', sumoutf, suminf, sigma(v, s, t));
                            checkfv=0;
                        end
                    end
                end
            end
            v=v+1;
        end
        
        %checkx
        %checkfx
        %checkf
        %checkfv
        
        %disp_f(f, S, T, V, beta);
        % recording feasible solutions
        if (checkx==1)&&(checkfx==1)&&(checkf==1)&&(checkfv==1)
            Z_SIZE = size(z);
            FEA_Z_SIZE = size(fea_z);
            
            fea_z(:,:,fea_idx)=z;
            %fea_x(:,:,:,fea_idx)=x;
            %fea_beta(:,:,:,fea_idx)=beta;
            %fea_f(:,:,:,:,fea_idx)=f;
            
            % cal sum cost for a feasible solution
            %Usum=0;
            %for e=1:E
            %    iv=real(edge(e));
            %    ov=imag(edge(e));
            %    Usum=Usum+z(iv,ov);
            %end
            Usum = GetCost(z, cost_mat);
            
            fea_usum(fea_idx)=Usum;
            fea_idx = fea_idx + 1;
            
            %Usum
            diary off;
        end
        return
    end
    
    si = real(st_imag(st_index));
    ti = imag(st_imag(st_index));
    
    final_mat = get_paths(TT(ti), S, SS, EE, V);
    paths = final_mat{si};
    width = max(find(~cellfun(@isempty, paths)));
    
    if length(width) == 0
        return
    end
    
    %permutation vector
    perm_mat = eye(width);
    
    %loop over the rows of the matrix; should create
    %the recursive stack
    for r = 1:width
        betas = beta;
        fs = f;
        
        row = perm_mat(r, :);
        
        one_index = find(row);
        path_vec = paths{one_index};
        
        for j = 1:length(path_vec)-1
            fs(path_vec(j), path_vec(j+1), si, ti) = 1;
        end
        
        for j = 1:length(path_vec)-2
            betas(path_vec(j), path_vec(j+1), path_vec(j+2)) = 1;
        end
        
        path_comb(st_imag, st_index+1, fs, betas, SS, TT, EE, ST, IV, OV, S, T, V, E, edge, sigma);
    end
end


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensionality reduction mechanism %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, cont]=explore_vars_x_iter(x, v, s, EE, SS, V, P, T, E, edge_imag)
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
            [x,cont_x] = explore_vars_x_iter(x, e, s, EE, SS, V, P, T, E, edge_imag);
            %this means we reached the terminal node when explored
            if cont_x
                %fprintf('Travelling down (x). i: %d, j: %d, s: %d\n', e, v, s)
                x(e, v, s) = 0;
                cont = cont_x;
            end
        end
    end
end

function x=explore_vars_x(x, EE, SS, V, P, T, E, TT, ST, edge_imag)
    for ti = 1:length(TT)
        for si = 1:length(SS)
            [x, cont] = explore_vars_x_iter(x, TT(ti), si, transpose(EE), SS, V, P, T, E, edge_imag);
        end
    end
end

function [f, cont]=explore_vars_f_iter(f, v, s, t, EE, SS, V, P, T, E, edge_imag)
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
            [f,cont_x] = explore_vars_f_iter(f, e, s, t, EE, SS, V, P, T, E, edge_imag);
            %this means we reached the terminal node when explored
            if cont_x
                %fprintf('Travelling down. i: %d, j: %d, s: %d, t: %d\n', v, e, s, t)
                f(e, v, s, t) = 0;
                %should set f to 0
                cont = cont_x;
            end
        end
    end
end

function f=explore_vars_f(f, EE, SS, V, P, T, E, TT, ST, edge_imag)
    
    tt_length = length(TT);
    %fprintf('Exploring f variables...\n')
    
    for ti = 1:length(TT)
        s_vec = ST(:, ti);
        for si = 1:length(s_vec);
            if s_vec(si) == 1
                t = TT(ti);
                [f, cont] = explore_vars_f_iter(f, t, si, ti, transpose(EE),...
                    SS, V, P, T, E, edge_imag);
            end
        end
    end
    
end

