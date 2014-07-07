function DFSPathCounting()
    EE = zeros(5, 5);
    EE(1, 2) = 1;
    EE(2, 4) = 1;
    EE(2, 3) = 1;
    EE(4, 5) = 1;
    EE(5, 3) = 1;

    view(biograph(EE));
    paths = explore(3, EE)
end

function [paths, cont, cell_mat, path_count]=explore_iter(s, t, EE, paths, cell_mat, path_count)
    edges = EE(s, :);
    paths = [paths, s];

    %edge exists between the starting point and terminal;
    %that means we are done with this iteration of DFS
    if s == t
        cont = 1;
        cell_mat(path_count + 1, :) = paths
        path_count = path_count + 1 
        return
    end

    cont = 0;
    edges
    for e = 1:length(edges)
        if edges(e) == 1 
            [paths, cont_x, cell_mat, path_count] = explore_iter(e, t, EE, paths, cell_mat, path_count);
            if cont_x == 1
                cont =  1;
            end 
        end
    end
end


function paths=explore(t, EE)
    cell_mat = cell(5, 5) 
    paths= explore_iter(1, t, EE, [], cell_mat, 0);
end

