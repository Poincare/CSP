%EE: adjacency matrix
%IV: input nodes cell array
function sorted_res=TopologicalSort(EE, IV)
    V = length(EE);
    DM = -ones(1, V);
    for i = 1:V
        iv = IV{i};
        if sum(iv) == 0
            DM(i) = 0;
        else
            DM(i) = length(iv);
        end
    end

    sorted_res = TopologicalSortIter(EE, DM, []);
end

%EE: adjacency matrix
%DM: degree matrix (initially 0's for source nodes, 1's for nonsource nodes)
%sorted: Sorted list of node indices (initially [])
function sorted_res =TopologicalSortIter(EE, DM, sorted_res)
    if length(DM(DM >= 0)) <= 0
        return
    end 

    min_val = min(DM(DM >= 0));
    source_node = find(DM == min_val);
    source_node = source_node(1);
    DM(source_node) = DM(source_node) - 1;
    for i = 1:length(EE(source_node, :))
        if EE(source_node, i) == 1
            DM(i) = DM(i) - 1;
        end
    end 

    sorted_res = [sorted_res, source_node];
    sorted_res = TopologicalSortIter(EE, DM, sorted_res);
end
