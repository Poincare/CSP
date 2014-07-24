function cost_mat=CreateCostMat(V)
    cost_mat = zeros(V, V);

    function randomCostMat()
        %randomize the cost matrix so that
        %we get some bottlenecks sometimes
        cost_mat = randi(20, V, V);       
    end

    function unitCostMat()
        cost_mat = ones(V, V);
    end

    randomCostMat();
end