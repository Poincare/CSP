function Benchmarking()
    V = 6;
    P = 2;
    T = 1;
    
    %vars = zeros(1, V*V*P*T + V*V*P);
    %vars(:) = rand;
    
    %sum = 0;
    %for i = 1:10000000
    %    x = get_f_from_vars(vars, 1, 1, 1, 1, V, P, T);
    %    sum = sum+ x;
    %end
    
    vars = zeros(2, V, V, P, T);
    vars(1, :, :, :, :) = rand;
    
    sum = 0;
    for i = 1:10000000
        x = vars(1, 1, 1, 1, 1);
        sum = sum + x;
    end
end

%get value of f from the long "vars" vector
function res=get_f_from_vars(vars, i, j, p, t, V, P, T)
    res = vars((i-1)*V*P*T + (j-1)*P*T + (p-1)*T + t);
end