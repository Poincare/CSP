EE = zeros(7, 7);
EE(1, 2) = 1;
EE(1, 3) = 1;
EE(2, 4) = 1;
EE(3, 4) = 1;
EE(4, 5) = 1;
EE(5, 6) = 1;
EE(5, 7) = 1;

view(biograph(EE));