
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
z=-ones(V,V);
x=-ones(V,V,S);
beta=-ones(V,V,V);
f=-ones(V,V,S,T);
x(1,3,1)=1;
x(1,3,2)=0;
x(2,5,2)=1;
x(2,5,1)=0;

Z=E; %Number of Edges
X=E*S;
Beta=[1,3,8];%Use to find x
F=E*sum(sum(ST)); %Flows for each terminal

beta(1,3,4)=1;
beta(3,4,6)=1;
beta(3,9,10)=1;
beta(4,6,7)=1;
beta(5,4,6)=1;
beta(4,6,10)=1;
beta(2,5,4)=1;
beta(2,5,7)=1;
beta(1,3,8)=1;
beta(1,3,9)=1;
beta(4,6,9)=0;
beta(6,9,10)=0;
beta(6,9,11)=0;
beta(9,11,8)=0;



for s=1:S %x on the left is the edge and the x on the right is the preceeding edge
        x(3,4,s)=beta(1,3,4)*x(1,3,s);
    end
    
    for s=1
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
                     