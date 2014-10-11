%%This file tries out a few simple topologies
%%to test CFL.m
clc
clear

%node 1 -------> node 2
V = 2;
P = 1;
S = P; %change of notation

%number of edges
E = 1;

%adj matrix
EE = zeros(V, V);
EE(1, 2) = 1;

%source(s)
SS=[1];

%terminal(s)
TT = [2];
T = length(TT);

%Flows Each Terminal Wants
ST=zeros(S,T);
ST(1, 1) = 1;

%set up OV and IV
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

%sigma
sigma=zeros(V,P,T);
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

%e
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
            for s=1:S
                x(l,m,s)=0;
                for t=1:T
                    if ST(s,t)==1
                        f(l,m,s,t)=0;
                    end
                end
            end
            for n=1:V
                if EE(m,n)==1
                    beta(l,m,n)=0;
                end
            end
        end
    end
end

f


CFL(V, P, T, EE, E, SS, OV, IV, sigma,ST, edge, TT)
