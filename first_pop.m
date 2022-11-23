function [new_population] = first_pop(dir, h0, x ,y, nPopulationSize)

% valores de busca para a direção de ondas
int_dir = [0,45,90,135,180,225,270,315];

% adiciona a população os parâmetros informados pelo usuário
new_population(1,:) = [dir h0 x y];

% taxa de variação para altura de ondas, parâmetros de profundidade e
% energia.
dx = 0.01:0.01:0.5;
k = 2;
while( size(new_population,1) < nPopulationSize )
    
        if rand>0.5
            new_alt = h0+dx(randi(length(dx), 1))*h0;
            new_x = x+dx(randi(length(dx), 1))*x;
            new_y = y+dx(randi(length(dx), 1))*y;
        else
            new_alt = h0-dx(randi(length(dx), 1))*h0;
            new_x = x-dx(randi(length(dx), 1))*x;
            new_y = y-dx(randi(length(dx), 1))*y;
        end
        
         if rand>0.5
            new_alt = h0+dx(randi(length(dx), 1))*h0;
            new_x = x+dx(randi(length(dx), 1))*x;
            new_y = y+dx(randi(length(dx), 1))*y;
        else
            new_alt = h0-dx(randi(length(dx), 1))*h0;
            new_x = x-dx(randi(length(dx), 1))*x;
            new_y = y-dx(randi(length(dx), 1))*y;
        end
        

    new_population(k,:) = [int_dir(randi(length(int_dir), 1)) new_alt new_x new_y];
    k = k+1;
end


