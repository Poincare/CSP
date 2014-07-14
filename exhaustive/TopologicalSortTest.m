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

ST=zeros(S,T); %Flows Each Terminal Wants
ST(1,1)=1;
ST(1,2)=1;
ST(1,3)=1;

TopologicalSort(EE, IV)
