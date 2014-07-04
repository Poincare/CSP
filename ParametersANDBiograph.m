function ParametersANDBiograph()
    clear
    clc

    %%Parameters
    V=11; %Number of Nodes
    SS=[1,2]; %Sources
    S=length(SS); %Number of Sources
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
    %view(biograph(VE)); 

    %input node
    IV=cell(V);
    for v=1:V
        %make first column 0
        IV{v}=0;

        %loop over cols
        for iv=1:V

            %if an edge exists from iv to v
            if EE(iv,v)==1

                %if there are no nodes yet, make a list
                if sum(IV{v})==0
                    IV{v}=iv;

                %otherwise, put into existing list
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



    z=-ones(V,V);
    x=-ones(V,V,S);
    beta=-ones(V,V,V);
    f=-ones(V,V,S,T);

    for l=1:V %column of EE matrix
        for m=1:V %row of EE matrix
            if EE(l,m)==1 %if l,m work in the EE matrix(when there is a 1 instead of 0)
                z(l,m)=0; %show 0 instead of -1
                %for s=1:S
                    %x(l,m,s)=0;
                    %for t=1:T
                    %    if ST(s,t)==1
                    %        f(l,m,s,t)=0;
                    %    end
                    %end
                %end
                for n=1:V
		  if EE(m,n)==1
                      beta(l,m,n)=0;
                    end
                end
            end
        end
    end
    
    f = explore_vars_f(f, EE, SS, V, S, T, E, TT, ST, edge);

    %reduces x dimensionality
    x = explore_vars_x(x, EE, SS, V, S, T, E, TT, ST, edge);
    
    return

    % x for source edges *************brute force way*****
    x(1,3,1)=1;
    x(2,5,2)=1;

    Z=E; %Number of Edges
    X=E*S;
    Beta=sum(sum(sum(beta==zeros(V,V,V))));%Use to find x
    F=E*sum(sum(ST)); %Flows for each terminal
    Tcomb=2^Beta*2^F; %Possible combinations: consider beta and f, othersdetermined by them

    fea_idx=1;
    for tcomb=1:2 %Tcomb
        %% assign values of a combination
        %for tcomb=1:2%Tcomb %Tcomb=possible combinations
        display('tcomb');
        tcomb
        display('fea_idx');
        fea_idx

        tcombvec=change_10_K_N(Tcomb-tcomb+1,Beta+F,2); %change to binary
        tcombidx=1; %tcombidx=Tcomb index
        for l=1:V
            for m=1:V
                for n=1:V
                    if beta(l,m,n)~=-1
                        beta(l,m,n)=tcombvec(tcombidx);
                        tcombidx=tcombidx+1;
                    end
                end
            end
        end
        %tcombidx==Beta+1

        for l=1:V
            for m=1:V
                for s=1:S
                    for t=1:T
                        if f(l,m,s,t)~=-1
                            f(l,m,s,t)=tcombvec(tcombidx);
                            tcombidx=tcombidx+1; 
                        end
                    end
                end
            end
        end

        %tcombidx==Beta+F+1

        %%test whether it is a feasible combination

        %generate x based on beta ********brute force way**********
        %Constraint (31)
        %step 1

        for s=1:S %x on the left is the edge and the x on the right is the preceeding edge
            x(3,4,s)=beta(1,3,4)*x(1,3,s);
        end

        for s=1:S
            x(3,8,s)=beta(1,3,8)*x(1,3,s);
        end

        for s=1:S
            x(3,9,s)=beta(1,3,9)*x(1,3,s);
        end

        for s=1:S
            x(5,4,s)=beta(2,5,4)*x(2,5,s);
        end

        for s=1:S
            x(5,7,s)=beta(2,5,7)*x(2,5,s);
        end
        %step 2

        for s=1:S
            x(4,6,s)=max(beta(3,4,6)*x(3,4,s),beta(5,4,6)*x(5,4,s));
        end

        %step 3
        for s=1:S
            x(6,7,s)=beta(4,6,7)*x(4,6,s);
        end

        for s=1:S
            x(6,10,s)=beta(4,6,10)*x(4,6,s);
         end

         for s=1:S
            x(6,9,s)=beta(4,6,9)*x(4,6,s);
         end

         %step 4
        for s=1:S
            x(9,11,s)=max(beta(6,9,11)*x(6,9,s),beta(3,9,11)*x(3,9,s));
        end

         for s=1:S
            x(9,10,s)=max(beta(6,9,10)*x(6,9,s),beta(3,9,10)*x(3,9,s));
         end

         %step 5
         for s=1:S
            x(11,8,s)=beta(9,11,8)*x(9,11,s);
         end

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
                         sumf=sumf+f(iv,ov,s,t);
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
         %Constraint (27)
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
                                 sumoutf=sumoutf+f(v,ov,s,t);
                             end
                         end

                         suminf=0;
                         if sum(IV{v}~=0)>=1
                             for ividx=1:length(IV{v})
                                 iv=IV{v}(ividx);
                                 suminf=suminf+f(iv,v,s,t);
                             end
                         end

                         if (sumoutf-suminf~=sigma(v,s,t))
                             checkfv=0;
                         end
                     end
                 end
             end
             v=v+1;
         end

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

             fea_Usum(fea_idx)=Usum;
         end     
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    