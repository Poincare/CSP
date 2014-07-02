V = 3;
EE = zeros(V, V);
EE(1,3) = 1;
EE(2, 3) = 1;
view(biograph(EE));
view(biograph(transpose(EE)));