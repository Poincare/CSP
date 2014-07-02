function TestParallel
    x = zeros(2, 2);
    spmd
        x(labindex) = labindex;
    end
    
    x
end

function res=sq_mat(mat)
    res = mat.^2;
end