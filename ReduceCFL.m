function [vars, p]=ReduceSearch(vars, p, v, V, P, T, OV)
    i = v;
    ov = OV{v};

    if sum(ov) == 0
        return
    end

    for ovi = 1:length(ov)  
        j = ov(ovi);

        found_one = 0;
        for p = 1:P
            for t = 1:T
                if get_f_from_vars(vars, i, j, p, t, V, P, T) == 1
                    found_one = 1;
                    break;
                end
            end

            if found_one
                vars = set_x_in_vars(vars, 1, i, j, p, V, P, T);
                p(linear_index_x(i, j, p, V, P, T), 1) = 0;
                p(linear_index_x(i, j, p, V, P, T), 2) = 1;
            end
        end
    end

    
end
