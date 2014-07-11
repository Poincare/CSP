function Atoms()
    clear
    clc

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

    %result matrices
    global z x beta f

    function initialize_params()
        V = 5;
        SS = [1, 2];
        S = length(SS);
        P = S;
        TT = [4, 5];
        T = length(TT);
        EE = zeros(V, V);
        EE(1, 3) = 1;
        EE(1, 5) = 1;
        EE(2, 3) = 1;
        EE(3, 4) = 1;
        EE(3, 5) = 1;
        E = sum(sum(EE));

        ST(1, 1) = 1;
        ST(2, 2) = 1; 
        %V=11; %Number of Nodes
        %SS=[1,2]; %Sources
        %S=length(SS); %Number of Sources
        %P = S;
        %TT=[7,8,10]; %Terminals
        %T=length(TT); %Number of Terminals
        %EE=zeros(V,V); %Edges
        %EE(1,3)=1; %Edges
        %EE(3,4)=1;
        %EE(3,8)=1;
        %EE(11,8)=1;
        %EE(5,4)=1;
        %EE(4,6)=1;
        %EE(6,7)=1;
        %EE(6,10)=1;
        %EE(6,9)=1;
        %EE(9,11)=1;
        %EE(2,5)=1;
        %EE(5,7)=1;
        %EE(9,10)=1;
        %EE(3,9)=1;
        %E=sum(sum(EE)); %Number of Edges
%
%        ST=zeros(S,T); %Flows Each Terminal Wants
%        ST(1,1)=1;
%        ST(1,2)=1;
%        ST(1,3)=1;
%        ST(2,1)=1;
%        ST(2,3)=1;

        %Graph the Network
        %VE=sparse([1 2 3 3 4 5 5 6 6 6 3 9 9 11],[3 5 4 8 6 4 7 7 10 9 9 11 10 8],true,11,11);
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

    initialize_params();
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

        t
        D{t} = sources;
    end

    ti_size = length(TP{3})
    atom_mat = d_atoms(3, D, 1:S, 1, cell(1, 1));
    atoms = cell(1, S);
    [rows, cols] = size(atom_mat);
    atoms = cell(1, 2^ti_size);
    i = 1
    for r = 1:rows
        if length(atom_mat{r}) ~= 0
            atoms{i} = atom_mat{r}  
            i = i + 1
        end
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


%st_index: index of the st_imag vector (should be 1 for the first call)
function path_comb(st_index, f, atoms)
    fprintf('\n--------------------------------\n')
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP

    %result matrices
    global z x beta
    
    %if the source is the last one in consideration,
    %then we can check conditions
    if st_index > length(st_imag)
        %CHECK ATOM FEASIBILITY HERE
        for i = 1:V
            iv = IV{i}
            for j = 1:length(iv);
                for m = 1:length(atoms)
                    atom_set = atoms(m);
                    for t = 1:T
                        equi = zeros(1, length(atom_set));
                        for si = 1:length(atom_set)
                            s = atom_set(si);
                            equi(s) = f(j, i, s, t);
                        end
                        
                        equi_r = equi(equi ~= -1)

                        %not all 1's or all 0's
                        if sum(equi_r) ~= 0 && sum(equi_r) ~= length(equi_r)
                            %don't do anything? do something if it all works out? 
                        end
                    end
                end
            end
        end
        return
    end 

    si = real(st_imag(st_index));
    ti = imag(st_imag(st_index));

    final_mat = get_paths(TT(ti), S, SS, EE, V);
    paths = final_mat{SS(si)};
    width = max(find(~cellfun(@isempty, paths)));

    %permutation vector
    perm_mat = eye(width);

    %loop over the rows of the matrix; should create
    %the recursive stack
    for r = 1:width
        %betas = b;
        fs = f;
        row = perm_mat(r, :);

        one_index = find(row);
        path_vec = paths{one_index};

        for j = 1:length(path_vec)-1
            fs(path_vec(j), path_vec(j+1), si, ti) = 1; 
        end

        %for j = 1:length(path_vec)-2
        %    betas(path_vec(j), path_vec(j+1), path_vec(j+2)) = 1;
        %end

        path_comb(st_index+1, fs);
    end
end

function atom_mat=d_atoms(i, D,atoms, ti, atom_mat) 
    global V SS S P TT T EE E ST IV OV sigma edge st_imag TP

    %result matrices
    global z x beta f

    terminals = TP{i};
    terminals

    %base case
    if ti > length(terminals)
        atom_mat = cell(1, 1);
        atom_mat{1} = atoms; 
        return
    end   

    atoms_l = intersect(atoms, D{terminals(ti)});  
    if length(atoms_l) ~= 0
        atom_mat_l = d_atoms(i, D, atoms_l, ti + 1, atom_mat);
        atom_mat_l
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
            fprintf('Found source. v: %d, s: %d, t:%d\n', v, s, t)
            cont = 1;
            return
        end
    end

    for e = 1:length(edges)
        if edges(e) == 1
            [f,cont_x] = explore_vars_f_iter(f, e, s, t, EE, SS, V, P, T, E, edge_imag);
            %this means we reached the terminal node when explored
            if cont_x
                fprintf('Travelling down. i: %d, j: %d, s: %d, t: %d\n', v, e, s, t)
                f(e, v, s, t) = 0;
                %should set f to 0
                cont = cont_x;
            end
        end
    end
end

function f=explore_vars_f(f, EE, SS, V, P, T, E, TT, ST, edge_imag)
    
    tt_length = length(TT);
    fprintf('Exploring f variables...\n')
    
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

