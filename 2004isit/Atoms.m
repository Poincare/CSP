function min_cost_z = Atoms(v, ss, s, tt, t, ee, e, st)
    %function Atoms()
    %graph parameters or basic info derived from parameters
    %V: number of vertices
    %SS: sources
    %S: # of sources
    %P: # of sources
    %TT: terminals
    %EE: adj matrix
    %E: number of edges
    %ST: flow-terminal matching
    %IV: input cell array
    %OV: output cell array
    %sigma: see paper
    %edge: vector of edges (real + imag representation)
    %st_imag: vector of ST (real + imag representation)
    %TP: the terminals each node can access
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
    
    V = v;
    SS = ss;
    S = s;
    P = S;
    TT = tt;
    T = t;
    EE = ee;
    E = e;
    ST = st;
    
    addpath('../exhaustive');
    
    
    %result matrices
    global z x beta f fea_f fea_idx fea_z
    
    function initialize_params()
        V = 3;
        SS = [1,2];
        S = length(SS);
        P = S;
        TT = [4];
        T = length(TT);
        EE = zeros(V, V);
        EE(1,3) = 1;
        EE(2,3) = 1;
        EE(3,4) = 1;
        
        ST(1,1) = 1;
        ST(2,1) = 1;
        
    end
    
    function initialize_info()
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
        
        %edge vector
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
        
    end
    
    %initialize_params();
    initialize_info();
    
    
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
                end
                for n=1:V
                    if EE(m,n)==1
                        beta(l,m,n)=0;
                    end
                end
            end
        end
    end
    
    
    %use DFS to set the correct variable spots for f.
    f = explore_vars_f(f, EE, SS, V, P, T, E, TT, ST, edge);
    
    %set up T(i) from 2004isit paper
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
        atom_mat = d_atoms(i, D, 1:S, 1, cell(1, 1));
        
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
    
    %x = zeros(4, 4, 2);
    %x(1,3,1) = 1;
    %x(2,3,2) = 1;
    
    fea_f = -ones(V, V, P, T, 1);
    fea_z = zeros(V, V, P, 1);

    fea_idx = 1;
    
    path_comb(1, f, beta, C);
    
    min_cost = Inf;
    min_cost_z = fea_z(:, :, 1);
   
     
    for idx = 1:length(fea_z(1, 1, :))
        cost = sum(sum(fea_z(:, :, idx)));
        if (cost < min_cost) && (cost > 0)
            min_cost = cost;
            min_cost_z = fea_z(:, :, idx);
        end
    end
end

%looks at the generated f's and determines the z values
function fea_z=generate_fea_z()
    global fea_f fea_idx
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
   
    fea_length = length(fea_f(1, 1, 1, 1, :));
    fea_z = -ones(V, V, fea_length);
    
    for idx = 1:fea_length
        z = zeros(V, V);
        f = fea_f(:, :, :, :, idx);
        if sum(sum(sum(sum(f)))) == -V*V*P*T
            continue; 
        end

        for e = 1:E
            iv = real(edge(e));
            ov = imag(edge(e));
            
            for t = 1:T
                for s = 1:P
                    f_val = f(iv, ov, s, t);
                    if f_val == 1
                        z(iv, ov) = 1;
                    end
                end
            end
        end
       
        fea_z(:, :, idx) = z
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

function res = check_path_dfs(v, si, ti, f)
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
    
    res = 0;
    if v > V || sum(OV{v}) == 0
        if v == TT(ti)
            res=1;
            return;
        else
            res = 0;
            return;
        end
    end
    
    ov = OV{v};
    for j = 1:length(ov)
        if f(v, ov(j), si, ti) == 1
            res_p = check_path_dfs(ov(j), si, ti, f);
            if res_p == 1
                res = 1;
                return
            end
        end
    end
end

function res = check_path(f)
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
    res = 1;
    
    for si = 1:S
        for ti = 1:T
            if ST(si, ti) == 1
                res_p = check_path_dfs(SS(si), si, ti, f);
                if res_p ~= 1
                    res = 0;
                    return;
                end
            end
        end
    end
end

function disp_f(f, S, T, V)
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
end

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
function x=compute_x(x, beta, EE, IV, OV, V, S, T)
    sorted_nodes = TopologicalSort(EE, IV);
    
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


%st_index: index of the st_imag vector (should be 1 for the first call)
function path_comb(st_index, f, beta, C)
    %fprintf('\n--------------------------------\n')
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
    
    %result matrices
    global z x fea_f fea_idx fea_z
    
    %if the source is the last one in consideration,
    %then we can check conditions
    if st_index > length(st_imag)
        %satisfied = 1;
        atom_s = 1;
        checkf = 1;
        checkfv = 1;
        
        x = set_x_on_initial_edges(x, f, OV, IV, V, S, T, SS, TT);
        x = compute_x(x, beta, EE, IV, OV, V, S, T);
        
        v=1;
        while (v <= V) && (checkfv==1)
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
                            fprintf('Failed on node: %d, t: %d, s: %d\n', v, t, s);
                            fprintf('Sumoutf: %d, suminf: %d, sigma: %d\n', sumoutf, suminf, sigma(v, s, t));
                            checkfv=0;
                        end
                    end
                end
            end
            v=v+1;
        end
        
        %check f %Consraint (17)
        e=1;
        zt=zeros(E,T);
        if (checkf == 1) && (checkfv == 1)
            z_c = zeros(V, V);
            while (e<=E)&&(checkf==1)&&(checkfv==1)
                iv=real(edge(e));
                ov=imag(edge(e));
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
                    if sumf>1
                        checkf=0;
                    end
                end
                z_c(iv,ov)=max(zt(e,:));
                e=e+1;
            end
            
           
            atom_s = 1;
            if (checkf == 1) && (checkfv == 1)
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
                            
                            Value_sets{c} = value_set;
                        end
                    end
                  
                    %C_j 
                    %Value_sets 
                    value_sums = cellfun(@sum, Value_sets);
                    TOTAL_VALUE = length(value_sums(value_sums > 0));
                    if TOTAL_VALUE > 1
                        atom_s = 0;
                    end
                end
            end
            
            %fprintf('-----------------\n');
            %checkf
            %atom_s
            %checkfv
            %fprintf('---------------------\n');
            
            if (checkf == 1) && (checkfv == 1) && (atom_s == 1)
                z = zeros(V, V);
                for i = 1:V
                    for j = 1:V
                        for s = 1:S
                            for t=1:T
                                f_val = f(i, j, s, t);
                                if f_val == 1
                                    z(i, j) = 1;
                                end
                            end
                        end
                    end
                end

                z;
                fea_f(:, :, :, :, fea_idx) = f;
                fea_z(:, :, fea_idx) = z;

                fea_idx = fea_idx + 1;
                save('satisfied_atoms.mat', 'f', 'x', 'beta', 'ST', 'EE', 'TT', 'SS', 'S', 'T', 'V');
            end
            
            return
        end
    end
    
    si = real(st_imag(st_index));
    ti = imag(st_imag(st_index));
    
    final_mat = get_paths(TT(ti), S, SS, EE, V);
    paths = final_mat{si};
    width = max(find(~cellfun(@isempty, paths)));
    
    if length(width) == 0
        return;
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
        
        path_comb(st_index+1, fs, betas, C);
    end
end

function atom_mat=d_atoms(i, D,atoms, ti, atom_mat)
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP
    
    %result matrices
    global z x beta f
    
    terminals = TP{i};
    
    %base case
    if ti > length(terminals)
        atom_mat = cell(1, 1);
        atom_mat{1} = atoms;
        return
    end
    
    atoms_l = intersect(atoms, D{terminals(ti)});
    if length(atoms_l) ~= 0
        atom_mat_l = d_atoms(i, D, atoms_l, ti + 1, atom_mat);
        if ~isempty(atom_mat_l)
            atom_mat = [atom_mat; atom_mat_l];
        end
    end
    
    whole = 1:S;
    d_comp = setdiff(whole, D{terminals(ti)});
    atoms_r = intersect(atoms, d_comp);
    if length(atoms_r) ~= 0
        atom_mat_r = d_atoms(i, D, atoms_r, ti + 1, atom_mat);
        if ~isempty(atom_mat_r)
            atom_mat = [atom_mat; atom_mat_r];
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

