%number of variables
D = 2;

%this would have to be done slightly differently
%if this algorithm was distributed
x = zeros(1, D);

%initialize p
p = zeros(D, 2);
p(1:D, 1:2) = 1/D;

%should be run in parallel
while 1
    done = 1;
    for i = 1:D
        %design variable
        b = 0.2;

        %number of variables
        D = 2;

        %set up "M" (in CSP paper)
        %a "1" in the ith entry means that x(i) participates in
        %the given clause
        clauses = zeros(1, D);
        clauses(1, 1) = 1;
        clauses(1, 2) = 1;
        clauses(2, 2) = 1;

        r = rand;

        %realize random bernoulli variable
        %remember that p(i, 1) is actually the probability of j=0
        if r <= p(i, 1)
            x(i) = 0;
        else
            x(i) = 1;
        end

        %evaluate the clauses and see if they are
        %satisfied
        satisfied = 1;
        if clauses(i, 1) == 1
            %checking clause1
            if x(1) ~= 1
                satisfied = 0;
            end
        end

        if clauses(i, 2) == 1
           if x(2) ~= 1
               satisfied = 0;
           end
        end

        disp(p);
        disp('Satisfied: ');
        disp(satisfied);
        fprintf('X: %d, I: %d\n', x(i), i);

        %if it is satisfied, we can clear the probabilities
        if satisfied
            p(i, :) = 0;
            p(i, x(i) + 1) = 1;
                
            break;
        %otherwise, we have to interpolate the distribution
        else
            p(i, :) = (1-b)*(p(i,x(i) + 1)) + b*(1/(D-1));
            p(i, x(i) + 1) = (1-b)*(p(i,x(i) + 1));
        end
    end
    
    %check if we are completely done
    for r = 1:D
        for c = 1:2
            if p(r, c) ~= 0 && p(r, c) ~= 1
                done = 0;
            end
        end
    end
    
    if done
       break 
    end
end

fprintf('FINAL X: ');
disp(x);