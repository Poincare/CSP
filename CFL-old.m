%E: number of edges
%P: number of flows
%T: number of terminal nodes
%EE: adjacency matrix
function vars=CFL(V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
    vars = zeros(1, V*V*P*T + V*V*P);
    
    %number of variables
    N = length(vars);
    
    %possible values of variables
    D = 2;
    
    %number of clauses
    %(17): flow in {0, 1} - every f participates
    C = 5;
    
    p = zeros(N, 2);
    p(1:N, 1:D) = 1/D;
  
    %vector of variable values
    %f: first V* V * P * T entries
    %x: next V * V * P entries
    %vars(:) = -1;
    vars(1:V*V*P*T) = -1;
    vars = remove_nonexistent_edges(vars, V, P, T, EE);
    vars = explore_vars_x(vars, EE, SS, V, P, T, E, TT, ST, edge);
    vars = explore_vars_f(vars, EE, SS, V, P, T, E, TT, ST, edge);
    disp_vars(vars, V, P, T, E, edge)
    vars_after = length(vars(vars >= 0))

    %clause matrix
    %i,j is marked as 1 if variable i participates in clause j
    clause_mat = zeros(N, C);
    
    %all f's participate in clause 1
    clause_mat(1:V*V*P*T, 1) = 1;
    
    %all f's participate in clause 2
    clause_mat(1:V*V*P*T, 2) = 1;
    
    %all variables participate in clause 3
    %TODO: Not sure if this is really the case
    clause_mat(:, 3) = 1;
    
    %all x variables participate in clause 4
    clause_mat(V*V*P*T+1:N, 4) = 1;
    
    %all x variables participate in clause 5
    clause_mat(V*V*P*T+1:N, 5) = 1;
    
    [vars, p, clause_mat] = set_source_links(vars, p, clause_mat, V, P, T, EE, SS);
    
    last_time = cputime;
    iter_counter = 0;
    
    while 1
        done = 1;
        for i = 1:N 
            %check if the variable is supposed to be considered
            if vars(i) == -1
                continue;
            end
            
            %design variable
            b = 0.2;

            r = rand;
            %realize random bernoulli variable
            %remember that p(i, 1) is actually the probability of j=0
            if r <= p(i, 1)
                vars(i) = 0;
            else
                vars(i) = 1;
            end

            %evaluate the clauses and see if they are
            %satisfied
            satisfied = 1;
            if clause_mat(i, 1) == 1
                %checking clause1
                if flow_limit(vars, V, P, T, E, edge, ST) ~= 1
                    %fprintf('Flow limit failed.\n')
                    satisfied = 0;
                end
            end
            
            if clause_mat(i, 2) == 1 && satisfied
                %checking clause2
                if flow_conservation(vars,V,P,T,OV,IV,sigma,ST) ~= 1
                    %fprintf('Flow conservation failed.\n')
                    satisfied = 0;
                end
            end

            if clause_mat(i, 3) == 1 && satisfied
                %checking clause 3
                if flow_x(vars, V, P, T, E, edge, ST) ~= 1
                    %fprintf('Flow x failed.\n')
                    satisfied = 0;
                end
            end
            
            if clause_mat(i, 4) == 1 && satisfied
                %checking clause 4
                if checkx(vars, V, P, T, ST, IV, TT) ~= 1
                    %fprintf('Checkx failed.\n')
                    satisfied = 0;
                end
            end
            
            if clause_mat(i, 5) == 1 && satisfied
                %checking clause 5
                if checkbeta(vars, V, P, T, IV, E, edge) ~= 1
                    %fprintf('Checkbeta failed.\n')
                    satisfied = 0;
                end
            end
            
            %if it is satisfied, we can clear the probabilities
            if satisfied
                p(i, :) = 0;
                p(i, vars(i) + 1) = 1;
                continue;
                
            %otherwise, we have to interpolate the distribution
            else
                t = p(i,vars(i)+1);
                for j = 1:D
                    p(i, j) = (1-b)*(p(i,j)) + (b/(D-1));
                end
                
                p(i, vars(i) + 1) = (1-b)*t;
            end
            %fprintf('\n------------\n')
                    
            if rem(iter_counter, 10000) == 0
                fprintf('Iter: %d\n', iter_counter);
                disp_vars_p(vars, V, P, T, E, edge);
                fprintf('%.10f\n', (cputime - last_time)*1000)
                last_time = cputime;
            end
        end

        %check if we are completely done
        for r = 1:N
            if vars(r) == -1
                continue
            end
            
            for c = 1:D
                if p(r, c) ~= 0 && p(r, c) ~= 1
                    done = 0;
                end
            end
        end
        
        if done
           break 
        end

    end
    
    fprintf('ANSWER:: \n')
    disp_vars(vars, V, P, T, E, edge)
end

function vars=explore_vars_x(vars, EE, SS, V, P, T, E, TT, ST, edge_imag)
    function [vars, cont]=explore_vars_x_iter(vars, v, s, EE, SS, V, P, T, E, edge_imag)
        cont = 0;
        edges = EE(v, :);
        
        if sum(edges) == 0
            if v == SS(s)
                fprintf('Found source (x). v: %d, s: %d\n', v, s)
                cont = 1;
                return
            end
        end
        
        for e = 1:length(edges)
            if edges(e) == 1
                [vars,cont_x] = explore_vars_x_iter(vars, e, s, EE, SS, V, P, T, E, edge_imag);
                %this means we reached the terminal node when explored
                if cont_x
                    fprintf('Travelling down (x). i: %d, j: %d, s: %d\n', e, v, s)
                    vars = set_x_in_vars(vars, 0, e, v, s, V, P, T);
                    disp_vars(vars, V, P, T, E, edge_imag)
                    %should set f to 0
                    cont = cont_x;
                end
            end
        end
    end

    for ti = 1:length(TT)
        for si = 1:length(SS)
            [vars, cont] = explore_vars_x_iter(vars, TT(ti), si, transpose(EE), SS, V, P, T, E, edge_imag)
            vars;
        end
    end
end

function vars=explore_vars_f(vars, EE, SS, V, P, T, E, TT, ST, edge_imag)
    function [vars, cont]=explore_vars_f_iter(vars, v, s, t, EE, SS, V, P, T, E, edge_imag)
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
                [vars,cont_x] = explore_vars_f_iter(vars, e, s, t, EE, SS, V, P, T, E, edge_imag);
                %this means we reached the terminal node when explored
                if cont_x
                    fprintf('Travelling down. i: %d, j: %d, s: %d, t: %d\n', v, e, s, t)
                    vars = set_f_in_vars(vars, 0, e, v, s, t, V, P, T);
                    disp_vars(vars, V, P, T, E, edge_imag)
                    %should set f to 0
                    cont = cont_x;
                end
            end
        end
    end
    
    tt_length = length(TT);
    fprintf('Exploring f variables...\n')
    
    for ti = 1:length(TT)
        s_vec = ST(:, ti);
        for si = 1:length(s_vec);
            if s_vec(si) == 1
                t = TT(ti);
                [vars, cont] = explore_vars_f_iter(vars, t, si, ti, transpose(EE), SS, V, P, T, E, edge_imag)
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

%clause 1: Constraint (17)
function checkf=flow_limit(vars, V, P, T, E, edge, ST)
    checkf = 1;
    for e = 1:E
        iv=real(edge(e));
        ov=imag(edge(e));
        for t=1:T
            sumf=0;
            for s=1:P
                fval = get_f_from_vars(vars, iv,ov,s,t, V,P,T);
                if ST(s,t)==1 && fval ~= -1
                    sumf=sumf + fval;
                end
            end
            if sumf>1
                checkf=0;
            end
        end
    end
end

%clause2: Constraint (28)
function checkfv=flow_conservation(vars, V, P, T, OV, IV, sigma, ST)
    checkfv = 1;
    for v=1:V
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
                        checkfv=0;
                    end
                end
            end
        end
    end
end

%clause 3: Constraint (29)
function checkfx=flow_x(vars, V, P, T, E, edge, ST)
    checkfx=1;
    for e=1:E
        iv=real(edge(e)); %input node (represented by real number)
        ov=imag(edge(e)); %output node (represented by imaginary number)
        for s=1:P
            for t=1:T
                if ST(s,t) == 1
                    fval = get_f_from_vars(vars,iv,ov,s,t,V,P,T);
                    xval = get_x_from_vars(vars,iv,ov,s,V,P,T);
                    if fval > xval && fval ~= -1 && xval ~= -1
                        checkfx=0;
                    end
                end
            end
        end
    end
end


%clause 4: Constraint (32)
function res=checkx(vars, V, P, T, ST, IV, TT)
    res = 1;
    for t = 1:T
         for s=1:P
             if (ST(s,t)~=1)
                 for iv=1:length(IV{TT(t)})
                     xval = get_x_from_vars(vars, IV{TT(t)}(iv),TT(t),s,V,P,T);
                     if xval ~=0 && xval ~= -1
                         res=0;
                     end
                 end
             end
         end 
     end
end

%clause 5: Constraint (33)
%this involves trying multiple beta values
function res=checkbeta(vars, V, P, T, IV, E, edges)
    res = 1;
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
    for e = 1:E
        i = real(edge(e));
        j = imag(edge(e));
        
        for p = 1:P
            for t = 1:T
                val = get_f_from_vars(vars, i, j, p, t, V, P, T);
                fprintf('f; i: %d, j: %d, p: %d, t: %d, val: %d\n', i, j, p, t, val)
            end
        end
    end
    
    for e = 1:E
        i = real(edge(e));
        j = imag(edge(e));
        for p = 1:P
            val = get_x_from_vars(vars, i, j, p, V, P, T);
            fprintf('x; i: %d, j: %d, p: %d, val: %d\n', i, j, p, val);
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