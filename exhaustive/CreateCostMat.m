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

    function oneTenMat()
        cost_mat = (randi(2, V, V) - ones(V,V))*9 + ones(V,V);
    end
    
    unitCostMat();
    cost_mat(10,6) = 20;
    cost_mat(10,5) = 20;
    cost_mat(9,4) = 20;
    
end