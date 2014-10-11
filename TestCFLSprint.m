clc
clear

%%Parameters
V=11; %Number of Nodes
SS=[12, 13]; %Sources
S=length(SS); %Number of Sources
P=S
TT=[14,15]; %Terminals
T=length(TT); %Number of Terminals


%VIRTUAL TERMINALS + SOURCES
V = 11+S+T;

EE=zeros(V,V); %Edges
EE(12, 11) = 1;
EE(13, 8) = 1;
EE(6, 14) = 1;
EE(9, 15) = 1;
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

E=sum(sum(EE)); %Number of Edges

ST=zeros(S,T); %Flows Each Terminal Wants
ST(1,1) = 1;
ST(1,2) = 1;
ST(2,1) = 1;
ST(2,2) = 1;

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


OptimalCFL(V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
