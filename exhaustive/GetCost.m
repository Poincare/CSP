
function cost = GetCost(z, cost_mat)
    % cost = sum(sum(z))

    %check if the cost matrix has been filled in
    if sum(sum(cost_mat)) == 0
        cost = sum(sum(z));
        return
    end

    cost_mat_size = size(cost_mat);
    [z_x, z_y] = size(z);

    cost_mat = cost_mat(1:z_x, 1:z_y);
    cost = sum(sum(cost_mat.*z));
end
