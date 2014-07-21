function ReadCCH(filename)
    file_h = fopen(filename, 'r');
    tline = fgetl(file_h)
    while tline ~= -1
        sscanf(tline, '%s
        tline = fgetl(file_h);
    end
end

