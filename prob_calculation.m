function [prob,new_facies] = prob_calculation(lambdas, indicators)
    
    prob = zeros(1,size(indicators,2));
    for i = 1 : size(indicators,2)
        prob(i) = sum(indicators(:,i).*lambdas) + (1-sum(lambdas))*0.5;
    end
    
    if unique(prob)==0.5000
        new_facies = 0;
        prob = 0;
    else
        [~,posicao] = max(prob);
        new_facies = posicao;
        prob = max(prob);
    end
    
end
    
   