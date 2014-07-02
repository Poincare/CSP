function DFS()
    V = 4;
    EE = zeros(5, 5);
    EE(1, 3) = 1;
    EE(2, 3) = 1;
    EE(1, 4) = 1;
    EE(4, 5) = 1;
    EE(3, 5) = 1;
    
    SS = [1, 2];
    s = 1;
    
    explore_wrapper(5, s, transpose(EE), SS, V)
end

function explore_wrapper(v, s, EE, SS, V)
    mat = zeros(V, V);
    
    function stack=explore(v, s, EE, SS, V)
        edges = EE(v, :);
        stack = [v];
        
        if sum(edges) == 0
            if v == SS(s)
                stack
            else
                %set all the right things to -1
                stack = [];
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

    explore(v, s, EE, SS, V)
end