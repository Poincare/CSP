V=11; %Number of Nodes
SS=[1,2]; %Sources
S=length(SS); %Number of Sources
P=S;
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
ST(2,1) =1;
ST(1,2)=1;
ST(1,3)=1;
ST(2, 3) = 1;


Atoms(V, SS, S, TT, T, EE, E, ST)
