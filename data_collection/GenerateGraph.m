%iteration: the current realization # (used to name variables file)j
%pairing_avg: average pairing between source and terminal
%S: number of sources
%T: number of terminals
function [cost_exhaustive, shortest_path_count, path_count, cost_exhaustive_no_expansion, cost_routing, cost_atoms]...
=GenerateGraph(iteration, pairing_avg, s, t)

%this is to get the exhaustive search code
%on the path so that we can use it from this file
addpath('../exhaustive/');
addpath('../2004isit/');

global V RS S RT T EE ST SS TT OV

%changed after virtual sources and terminals are added
V = 14;

%number of sources
S=s; 
P = S;

%number of terminals
T=t; 

%THIS RANDOMIZES THE SOURCES AND TERMINALS - NEED IN SIMULATION
RS = generateRS(1, int64(V/2))
RT = generateRT(int64(V/2)+int64(1), V)

virtuals()

ST = zeros(S, T);
ST = generateST(pairing_avg);

fea_idx = 1;
fea_z = zeros(V, V, V);

EE = generateNSFNET();
EE = setVirtualLinks(EE);
E = sum(sum(EE));

OV = computeOV(V, EE);

%we should finally have a acyclic, directed graph
for si = 1:S
    makeAcyclic(SS(si), zeros(1, V), zeros(1, 1));
end

%save the realization to a file so that it can be
%used for other trials
foldername = strcat('realizations', '-', num2str(S), '-', num2str(T), '-', num2str(pairing_avg));
mkdir(foldername);
filename = strcat(foldername, '/realizations', num2str(iteration), '.vars');
filename
save(filename, 'RS', 'RT', 'ST', 'EE');

disp_EE(EE, V);

[z_exhaustive, shortest_path_count, path_count] = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST, 1);
%fprintf('Mixing\n');
if sum(sum(z_exhaustive)) ~= 0
    cost_exhaustive = getCost(z_exhaustive);
else
    cost_exhaustive = -1;
end

%do not expand the demand set in order to enlarge feasibility range
z_exhaustive_no_expand = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST, 0);
if sum(sum(z_exhaustive_no_expand)) ~= 0
    cost_exhaustive_no_expansion = getCost(z_exhaustive_no_expand);
else
    cost_exhaustive_no_expansion = -1;
    end

    z_routing = ExhaustiveSearchRouting(V, SS, S, TT, T, EE, E, ST);  
    %fprintf('Routing\n');
    if sum(sum(z_routing)) ~= 0
        cost_routing = getCost(z_routing);
    else
        cost_routing = -1;
    end

    z_atoms = Atoms(V, SS, S, TT, T, EE, E, ST);
    %disp_z(z_atoms, V);
    %fprintf('Atoms\n');
    if sum(sum(z_atoms)) ~= 0
        cost_atoms = getCost(z_atoms);
    else
        cost_atoms = -1;
    end

    %fprintf('Cost exhaustive: %d\n', cost_exhaustive);
    %fprintf('Cost routing: %d\n', cost_routing);
    %fprintf('Cost atoms: %d\n', cost_atoms);

end

function disp_EE(EE, V)
    for i = 1:V
        for j = 1:V
            if EE(i, j) == 1
                fprintf('i :%d, j: %d, val: %d\n', i, j, EE(i,j));
            end
        end
    end 
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

function cost = getCost(z)
    cost = sum(sum(z));
end

function RS = generateRS(Vs, Ve)
    global V S RT T EE ST SS TT OV

    RS = [randperm(round(Ve)-round(Vs), S)] + Vs; 
end

function RT = generateRT(Vs, Ve)
    global V S RT T EE ST SS TT OV

    RT = [randperm(round(Ve)-round(Vs), T)] + double(Vs);
end

function ST = generateST(pairing_avg)
    global V RS S RT T EE SS TT OV

    ST = zeros(S, T);

    %every terminal is first paired with one node
    %(which is the first element of the sources sequence)
    for ti = 1:T
        ST(1, ti) = 1;
    end 
    
    P = (pairing_avg - 1)/(S - 1); 
    for ti = 1:T
        for si = 1:S
            r = rand;
            if r < P
                ST(si, ti) = 1;
            end    
        end
    end
end

function EE = directionalize(path_cells, EE)
    %take one of the paths and directionalize according to it
    for j = 1:length(path_cells)
        path = path_cells{j};

        for i = 1:length(path)-1
            EE(path(i), path(i+1)) = 1;                 
        end
    end
end

function makeAcyclic(v, visited, cycle_stack)
    global V RS S RT EE T ST SS TT OV

    visited(v) = 1;
 
    if sum(OV{v}) == 0
        return
    end
 
    ov = OV{v};
    for i = 1:length(ov)
        if visited(ov(i)) == 1
            EE(v, ov(i)) = 0;
%            j = find(cycle_stack == ov(i), v, 'last');
%            cycle = cycle_stack(j:length(cycle_stack));
%            cycle = [cycle, v];
%            cycle = [cycle, ov(i)];
%
%            for k = 1:length(cycle)-1
%                m = cycle(k);
%                n = cycle(k+1);
%
%                for ti = 1:T
%                    paths_t = get_paths_directed(TT(t), S, SS, EE, V);
%
%                    for si = 1:S 
%                        paths = paths_t{si};
%                        if ST(si, ti) == 1
%                            
%                        end
%                    end
%                end
%            end
            
        else
            cycle_stack = [cycle_stack, v];
            makeAcyclic(ov(i), visited, cycle_stack);
        end 
    end 
end

%set links to/from virtual sources/terminals
function EE=setVirtualLinks(EE)
    global V RS S RT T OV ST SS TT

    for si = 1:S
        EE(SS(si), RS(si)) = 1; 
    end

    for ti = 1:T
        EE(RT(ti), TT(ti)) = 1;
    end

end

function IV=computeIV(V, EE)
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
end

function OV=computeOV(V, EE)
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
end

function EE=generateSprint()
    global V RS S RT T EE ST SS TT
    UE = zeros(V, V);
    UE(1,5) =1;
    UE(1,8) = 1;
    UE(2, 11) = 1;
    UE(2, 7) = 1;
    UE(4, 9) = 1;
    UE(4, 5) = 1;
    UE(5, 6) = 1;
    UE(5, 9) = 1;
    UE(5, 10) = 1;
    UE(5, 11) = 1;
    UE(6, 7) = 1;
    UE(7, 11) = 1;
    UE(7, 8) = 1;
    UE(8, 9) = 1;
    UE(8, 11) = 1;
    UE(9, 10) = 1;
    UE(10, 11) = 1;

    EE = generateDirected(UE);
end

function EE=generateDirected(UE)
    global V RS S RT T EE ST SS TT

    UE = transpose(UE) + UE;
    
    IV = computeIV(V, UE);
    OV = computeOV(V, UE); 

    EE = zeros(V, V);

    for si = 1:S
        for ti = 1:T
            if ST(si, ti) == 1
                path_cells = get_paths(UE, OV, RS(si), ti, zeros(1, V));
 
                %take one of the paths and directionalize according to it
                for j = 1:length(path_cells)
                    path = path_cells{j};

                    for i = 1:length(path)-1
                        if EE(path(i+1), path(i)) ~= 1
                            EE(path(i), path(i+1)) = 1;   
                        end
                    end
                end
            end
        end
    end    
end

function EE=generateNSFNET()
    global V RS S RT T EE ST SS TT
    
    %undirected version of the NSFNET topology
    UE = zeros(V, V);
    UE(1, 2) = 1;
    UE(1, 9) = 1;
    UE(2, 4) = 1;
    UE(2, 3) = 1;
    UE(3, 6) = 1;
    UE(4, 5) = 1;
    UE(4, 11) = 1;
    UE(4, 5) = 1;
    UE(5, 6) = 1;
    UE(5, 7) = 1;
    UE(6, 8) = 1;
    UE(6, 13) = 1;
    UE(7, 9) = 1;
    UE(8, 10) = 1;
    UE(9, 10) = 1;
    UE(10, 12) = 1;
    UE(10, 14) = 1;
    UE(11, 14) = 1;
    UE(12, 13) = 1;
    UE(13, 14) = 1;

    EE = generateDirected(UE);
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
function final_mat=get_paths_directed(t, S, SS, EE, V)
    final_mat = cell(1, S);
    for i = 1:S
        cell_mat = cell(i, V);
        [cont, cell_mat, path_count, paths] = get_paths_explore_iter(...
            SS(i), t, EE, cell_mat, 1, []);
        final_mat{i} = cell_mat;
    end
end

function m = get_paths(UE, OV, v, ti, visited)
    global V RS S RT T ST SS TT
    %v

    visited(v) = 1;
  
    if v == RT(ti)
        m = cell(1, 1);
        m{1} = v;
        %fprintf('HIT, v: %d, RT(ti): %d\n', v, RT(ti));
        return
    end

    ov = OV{v};
    
%    if sum(ov) == 0
%        m = cell(0);
%        return
%    end
%
    curr = 1;
    m = cell(0, 0);

    if sum(ov) ~= 0
        for j = 1:length(ov)
            visiteds = visited;

            if visited(ov(j)) ~= 1
                m_j = get_paths(UE, OV, ov(j), ti, visiteds);
                for k = 1:length(m_j)
                    m_j{k} = [v, m_j{k}];
                end
                m = [m; m_j];
            end
        end
    end
end

function virtuals()
    global V RS S RT T EE ST SS TT

    SS = zeros(1, S); 
    for s = 1:S
        SS(s) = V+1; 
        V = V+1;
    end

    TT = zeros(1, T);
    for t = 1:T
        TT(t) = V+1;
        V = V+1;
    end

    EE = zeros(V, V);
    for s = 1:S
        EE(SS(s), RS(s)) = 1; 
    end

    for t = 1:T
        EE(RT(t), TT(t)) = 1;
    end

end
