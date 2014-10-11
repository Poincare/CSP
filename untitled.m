load('p_cell_mat_fun.mat', '-mat');

V = 8;
i = 3;
j = 4;
S = 2;
P = S;
T = 2;

i1 = 77;
i2 = 78;
i3 = 79;
i4 = 80;

i1_values = [];
i2_values = [];
i3_values = [];
i4_values = [];

for pi = 1:length(p_cell_mat)
    i1_values = [i1_values, p_cell_mat{pi}(i1,2)];
    i2_values = [i2_values, p_cell_mat{pi}(i2, 2)];
    i3_values = [i3_values, p_cell_mat{pi}(i3, 2)];
    i4_values = [i4_values, p_cell_mat{pi}(i4, 2)];
end

i1_len = length(i1_values);
i2_len = length(i2_values);
i3_len = length(i3_values);
i4_len = length(i4_values);

len = 15000

plot(1:i1_len, i1_values, 1:i2_len, i2_values, 1:i3_len, i3_values, 1:i4_len, i4_values)

i1_values_m = [i1_values, ones(1, len - i1_len)* i1_values(i1_len)] ;
i2_values_m = [i2_values, ones(1, len - i2_len)* i2_values(i2_len)] ;
i3_values_m = [i3_values, ones(1, len - i3_len)* i3_values(i3_len)] ;
i4_values_m = [i4_values, ones(1, len - i4_len)* i4_values(i4_len)] ;

plot(1:len, i1_values_m, 1:len, i2_values_m, 1:len, i3_values_m, 1:len, i4_values_m)


