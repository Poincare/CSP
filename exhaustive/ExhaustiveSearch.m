function ExhaustiveSearch()
    global m 
    global fea_idx
    global fea_x
    global fea_beta
    global fea_z
    global fea_f
    global fea_usum;
    
    fea_idx = 1;
    m=1;

%    V = 3;
%    SS=[1,2];
%    S = length(SS);
%    P = S;
%    TT = [3];
%    T = length(TT);
%    EE = zeros(V, V);
%    EE(1, 3) = 1;
%    EE(2, 3) = 1;
%    E = sum(sum(EE));

%    ST = zeros(S, T);
%    ST(1, 1) = 1;
%    ST(2, 1) = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up parameters                     %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V=11; %Number of Nodes
    SS=[1,2]; %Sources
    S=length(SS); %Number of Sources
    P = S;
    TT=[7,8,10]; %Terminals
    T=length(TT); %Number of Terminals
    EE=zeros(V,V); %Edges
    EE(1,3)=1; %Edges
    EE(3,4)=1;
    EE(3,8)=1;
    EE(11,8)=1;
    EE(5,4)=1;
    EE(4,6)=1;
    EE(6,7)=1;
    EE(6,10)=1;
    EE(6,9)=1;
    EE(9,11)=1;
    EE(2,5)=1;
    EE(5,7)=1;
    EE(9,10)=1;
    EE(3,9)=1;
    E=sum(sum(EE)); %Number of Edges

    ST=zeros(S,T); %Flows Each Terminal Wants
    ST(1,1)=1;
    ST(1,2)=1;
    ST(1,3)=1;
    ST(2,1)=1;
    ST(2,3)=1;

    %Graph the Network
    VE=sparse([1 2 3 3 4 5 5 6 6 6 3 9 9 11],[3 5 4 8 6 4 7 7 10 9 9 11 10 8],true,11,11);

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
                edge(e)=l+m*i;
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

    %use DFS to set the correct variable spots for f.    
    f = explore_vars_f(f, EE, SS, V, P, T, E, TT, ST, edge);
    path_comb(st_imag, 1, f, beta, SS, TT, EE, ST, IV, OV, S, T, V, E, edge, sigma, z);

    diary 'optimal.txt'
    find_optimal()
    diary off
end


function find_optimal()
    global fea_usum fea_z fea_beta fea_x fea_f

    min_cost_idx = find(fea_usum, min(fea_usum))
    z = fea_z(:, :, min_cost_idx)
    beta = fea_beta(:, :, min_cost_idx)
    x = fea_x(:, :, :, min_cost_idx)
    f = fea_f(:, :, :, min_cost_idx)
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
                fprintf('SS(s): %d, ov(i): %d, s: %d, t: %d, f_val: %d\n', SS(s), ov(i), s, t, f(s, ov(i), s, t));
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

                        fprintf('beta(%d, %d, %d) * x(%d, %d, %d)\n', iv(k), v, ov(j), iv(k), v, s);
                        x_val = beta(iv(k), v, ov(j)) * x(iv(k), v, s);
                        if x_val > max_x
                            max_x = x_val;
                        end 
                    end
                    fprintf('value set, s: %d\n', s);
                    x(v, ov(j), s) = max_x 
                end
            end
        end
    end

end

function path_comb(st_imag, st_index, f, beta, SS, TT, EE, ST, IV, OV, S, T, V, E, edge, sigma, z)
    fprintf('\n--------------------------------\n')

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

        %TODO generalize this
        % x for source edges *************brute force way*****
        x = set_x_on_initial_edges(x, f, OV, IV, V, S, T, SS, TT)

        %TODO generalize this
        %generate x based on beta ********brute force way**********
        %Constraint (31)
        %step 1

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
         while (checkx==1)&&(checkfx==1)&&(e<=E)&&(checkf==1)
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
                            fprintf('Failed on node: %d, t: %d, s: %d\n', v, t, s);
                            fprintf('Sumoutf: %d, suminf: %d, sigma: %d\n', sumoutf, suminf, sigma(v, s, t));
                             checkfv=0;
                         end
                     end
                 end
             end
             v=v+1;
         end

        checkx
        checkfx
        checkf
        checkfv

        disp_f(f, S, T, V, beta);
         % recording feasible solutions
         if (checkx==1)&&(checkfx==1)&&(checkf==1)&&(checkfv==1)
             fea_z(:,:,fea_idx)=z;
             fea_x(:,:,:,fea_idx)=x;
             fea_beta(:,:,:,fea_idx)=beta;
             fea_f(:,:,:,:,fea_idx)=f;

             % cal sum cost for a feasible solution
             Usum=0;
             for e=1:E
                 iv=real(edge(e));
                 ov=imag(edge(e));
                 Usum=Usum+z(iv,ov);
             end

             fea_usum(fea_idx)=Usum;
             fea_idx = fea_idx + 1;

            Usum
            diary off;
         end     
        return
    end 

    si = real(st_imag(st_index));
    ti = imag(st_imag(st_index));

    final_mat = get_paths(TT(ti), S, SS, EE, V);
    paths = final_mat{SS(si)};
    paths
    width = max(find(~cellfun(@isempty, paths)));
    width

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
            fprintf('Found source (x). v: %d, s: %d\n', v, s)
            cont = 1;
            return
        end
    end
        
    for e = 1:length(edges)
        if edges(e) == 1
            [x,cont_x] = explore_vars_x_iter(x, e, s, EE, SS, V, P, T, E, edge_imag);
            %this means we reached the terminal node when explored
            if cont_x
                fprintf('Travelling down (x). i: %d, j: %d, s: %d\n', e, v, s)
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
    
     
