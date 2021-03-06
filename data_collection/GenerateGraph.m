%iteration: the current realization # (used to name variables file)j
%pairing_avg: average pairing between source and terminal
%S: number of sources
%T: number of terminals
%scheme_name: "NSFNET" or "SPRINT"
function [cost_exhaustive, shortest_path_count, path_count, cost_exhaustive_no_expansion, cost_routing, cost_atoms, ST_classification, RT, complexity]...
=GenerateGraph(iteration, pairing_avg, s, t, scheme_name, scheme_suffix)
%function cost_atoms=GenerateGraph(iteration, pairing_avg, s, t, scheme_name)

clearvars -except iteration pairing_avg s t scheme_name scheme_suffix

%this is to get the exhaustive search code
%on the path so that we can use it from this file
addpath('../exhaustive/');
addpath('../2004isit/');

global V RS S RT T EE ST SS TT OV

%changed after virtual sources and terminals are added
if strcmp(scheme_name, 'NSFNET') ~= 0
    V = 14;
else
    V=11;
end

%number of sources
S=s; 
P =S;

%number of terminals
T=t; 

%THIS RANDOMIZES THE SOURCES AND TERMINALS - NEED IN SIMULATION
%RS = generateRS(8,14);

%RT = generateRT(2, 4);
RT = randsample([2,3,4], T);
%RT = [2, 3];
% RT = [1,3];
%RT
%RS
%RT = [2,4];
%RT
if strcmp(scheme_name, 'NSFNET')
    RS = [14, 13];
else
    RS = [11, 8];
end

RT
RS

%RS = generateRS(11, 14);

    % V=11; %Number of Nodes
    % RS=[1,2]; %Sources
    % S=length(SS); %Number of Sources
    % T=[7,8,10]; %Terminals
    % T=length(TT); %Number of Terminals

    %SPRINT
    %RS = [11, 10];
    %RT = [1, 2, 3];
    
    %GRID
    %RS = [1,4];
    %RT = [6, 8, 9];
    
virtuals()

ST = zeros(S, T);
ST = generateST(pairing_avg);
%ST(1,1) = 1;
%ST(1,2) = 1;
%ST(2, 1) = 1;
%ST(2,2) = 1;

%ST(1,2)=1;
%ST(1,2) =1;
%ST(2,1) = 1;
%ST(2,2) = 1;

 %ST(2, 2) = 1;
% ST(1, 3) = 1;
% ST(1, 1) = 1;
% ST(2, 1) = 1;
% ST(1, 2) = 1;
% ST(1, 3) = 1;
% ST(2, 3) = 1;


    % ST(1,1)=1;
    % ST(1,2)=1;
    % ST(1,3)=1;
    % ST(2,1)=1;
    % ST(2,3)=1;

    %SPRINT
    %ST(1, 1) = 1;
    %ST(2, 1) = 1;
    %ST(1, 3) = 1;
    %ST(2, 3) = 1;
    %ST(1, 2) = 1;
    %ST(2, 2) = 1;

fea_idx = 1;
fea_z = zeros(V, V, V);

% if scheme_name == 'NSFNET'
%     EE = generateNSFNET();
% else
%     EE = generateSprint();
% end
% EE = generateGridTopology();

if strcmp(scheme_name, 'SPRINT') && strcmp(scheme_suffix, 'DIRECTED')
    %sprint directionality
    EE = zeros(V, V);
    EE(1,2) = 1;
    EE(1, 3) = 1;
    EE(4,1) = 1;
    EE(4,5) = 1;
    EE(5,1) = 1;
    EE(6,3) = 1;
    EE(7,6) = 1;
    EE(7,4) = 1;
    EE(8,6) = 1;
    EE(9, 2) = 1;
    EE(9,4) = 1;
    %EE(9,7) = 1;
    EE(7, 9) = 1;
    
    EE(10, 5) = 1;
    EE(10, 7) = 1;
    EE(10, 6) = 1;
    %EE(10, 8) = 1;
    EE(8, 10) = 1;
    EE(11, 9) = 1;
    EE(11, 10) = 1;

elseif strcmp(scheme_name, 'NSFNET') && strcmp(scheme_suffix, 'DIRECTED')
    EE = zeros(V, V);
    EE(14, 10) = 1;
    EE(14, 11) = 1;
    %EE(14, 12) = 1;
    %EE(13, 10) = 1;
    EE(13, 12) = 1;
    EE(13, 11) = 1;
    EE(12, 6) = 1;
    EE(11, 9) = 1;
    EE(11, 8) = 1;
    EE(10, 4) = 1;
    EE(9, 6) = 1;
    EE(8, 7) = 1;
    EE(8, 2) = 1;
    EE(7, 5) = 1;
    EE(6, 5) = 1;
    EE(6, 3) = 1;
    EE(5, 4) = 1;
    EE(4, 1) = 1;
    %EE(2, 3) = 1;
    EE(1, 3) = 1;
    EE(1, 2) = 1;
end

ST_classification = ClassifyST();

    % EE=zeros(V,V); %Edges
    %     EE(1,3)=1; %Edges
    %     EE(3,4)=1;
    %     EE(3,8)=1;
    %     EE(11,8)=1;
    %     EE(5,4)=1;
    %     EE(4,6)=1;
    %     EE(6,7)=1;
    %     EE(6,10)=1;
    %     EE(6,9)=1;
    %     EE(9,11)=1;
    %     EE(2,5)=1;
    %     EE(5,7)=1;
    %     EE(9,10)=1;
    %     EE(3,9)=1;



EE = setVirtualLinks(EE);
E = sum(sum(EE));

global cost_mat
cost_mat = CreateCostMat(V);

OV = computeOV(V, EE);

%we should finally have a acyclic, directed graph
for si = 1:S
    makeAcyclic(SS(si), zeros(1, V), zeros(1, 1));
end

%save the realization to a file so that it can be
%used for other trials
foldername = strcat('realizations-', scheme_name, '-', scheme_suffix, '-', num2str(S),...
    '-', num2str(T), '-', num2str(pairing_avg));

mkdir(foldername);
filename = strcat(foldername, '/realizations', num2str(iteration), '.vars');
filename
save(filename, 'RS', 'RT', 'ST', 'EE');

complexity = zeros(1, 4);


[z_exhaustive, shortest_path_count, path_count, check_constraint_counter] = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST, 1, 0, 0, cost_mat);
%fprintf('Mixing\n');
if GetCost(z_exhaustive, cost_mat) ~= 0
    %z_exhaustive
    cost_exhaustive = GetCost(z_exhaustive, cost_mat);
    
    complexity(1, 1) = check_constraint_counter;
else
    cost_exhaustive = -1;
end

%do not expand the demand set in order to enlarge feasibility range
[z_exhaustive_no_expand, ~, ~, check_constraint_counter] = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST, 0, 0, 0, cost_mat);
if GetCost(z_exhaustive_no_expand, cost_mat) ~= 0
    cost_exhaustive_no_expansion = GetCost(z_exhaustive_no_expand, cost_mat);
    complexity(1, 2) = check_constraint_counter;
else
    cost_exhaustive_no_expansion = -1;
    end

    [z_routing,~,~, check_constraint_counter] = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST, 0, 1, 0, cost_mat);  
    fprintf('Routing\n');
    if GetCost(z_routing, cost_mat) ~= 0
        %z_routing
        cost_routing = GetCost(z_routing, cost_mat);
        complexity(1, 3) = check_constraint_counter;
    else
        cost_routing = -1;
    end

    [z_atoms,~,~,check_constraint_counter]  = ExhaustiveSearch(V, SS, S, TT, T, EE, E, ST, 0, 0, 1, cost_mat);
    %disp_z(z_atoms, V);
    %fprintf('Atoms\n');
    if GetCost(z_atoms, cost_mat) ~= 0
        cost_atoms = GetCost(z_atoms, cost_mat);

        complexity(1, 4) = check_constraint_counter;
    else
        cost_atoms = -1;
    end

    
    ST
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

function RS = generateRS(Vs, Ve)
    global V S RT T EE ST SS TT OV

    RS = [randperm(round(Ve)-round(Vs), S)] + Vs; 
end

function RT = generateRT(Vs, Ve)
    global V S RT T EE ST SS TT OV

    RT = [randperm(round(Ve)-round(Vs) + 1, T)] + round(Vs) - 1;
end

function ST = generateST(pairing_avg)
    global V RS S RT T EE SS TT OV

    ST = zeros(S, T);

    %every terminal is first paired with one node
    %(which is the first element of the sources sequence)
    for ti = 1:T
        ST(randi(S), ti) = 1;
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

function EE=generateGridTopology()
    global V RS S RT T EE ST SS TT
    
    V = 14;
    EE = zeros(V, V);
    EE(1,2) = 1;
    EE(1,4) = 1;
    EE(1,5) = 1;
    EE(2,3) = 1;
    EE(2,5) = 1;
    EE(2,6) = 1;
    EE(3,6) = 1;
    EE(4,5) = 1;
    EE(4,7) = 1;
    EE(4, 8) = 1;
    EE(5,6) = 1;
    EE(5,8) = 1;
    EE(5,9) = 1;
    EE(6,9) = 1;
    EE(7,8) = 1;
    EE(8, 9) = 1;
    
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
