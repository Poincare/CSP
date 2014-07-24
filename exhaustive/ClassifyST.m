%decides whether ST represents a unicast, multicast or general
%set of connections
function classification=ClassifyST()
    global S T ST

    D = cell(1, T);
    for ti = 1:T
        sources = [];

        for si = 1:S
            if ST(si, ti) == 1
                sources = [sources, si];
            end
        end

        D{ti} = sources;
    end

    %compute intersection
    whole = D{1};
    intersected = whole;
    for ti = 1:T
        intersected = intersect(intersected, D{ti});
    end

    if isequal(whole, intersected)
        classification = 'multicast';
    else if isempty(intersected)
        classification = 'unicast';
    else
        classification = 'general';
    end
end