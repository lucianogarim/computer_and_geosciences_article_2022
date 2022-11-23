function distance = distance_calculation (X,Y)
    
    distance = zeros(length(X), length(Y));
    
    for i = 1 : length(X)
        for j = 1 : length(Y)
            distance(i,j) = sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
        end
    end
    
end