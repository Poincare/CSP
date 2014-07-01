m = zeros(100, 100);
for i = 1:size(m, 1)
    for j = 1:size(m, 2)
        m(i, j) = 2*i;
    end
end
imagesc(m)

        