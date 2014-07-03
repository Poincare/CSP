clc

V = 3;
P = 1;
S = P ; %change of notation

%number of edges
E = 2;

%adj matrix
EE = zeros(V, V);
EE(1, 3) = 1;
EE(1, 2) = 1;

%source(s)
SS=[1];

%terminal(s)
TT = [2, 3];
T = length(TT);

%Flows Each Terminal Wants
ST=zeros(S,T);
ST(1, 1) = 1;
ST(1, 2) = 1;

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

IV

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
            edge(e)=l+m*1i;
            e=e+1;
        end
    end
end

edge

CFL(V, P, T, EE, E, SS, OV, IV, sigma, ST, edge, TT)
