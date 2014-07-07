function ExhaustiveSearch()
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
    view(biograph(VE)); 

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
    f = explore_vars_f(f, EE, SS, V, P, T, E, TT, ST, edge)

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
            [x, cont] = explore_vars_x_iter(x, TT(ti), si, transpose(EE), SS, V, P, T, E, edge_imag)
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
                [f, cont] = explore_vars_f_iter(f, t, si, ti, transpose(EE), SS, V, P, T, E, edge_imag);
            end
        end
    end
    
end
    
     
