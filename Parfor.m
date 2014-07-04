function Parfor()
    vars = zeros(1, 10);
    satisfied = 0;
    
    while ~satisfied
        for i = 1:10
            spmd
                vars(i) = rand
                if clause(vars) ~= 1
                    satisfied = 1;
                end
            end
        end
    end
    vars
end

function res = clause(vars)
    if vars(1) < 0.5
        res = 1;
    else
        res = 0;
    end
end