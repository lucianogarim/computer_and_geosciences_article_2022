function [geracao] = reproduction(selecao, probabilityOfCrossOver,p_cross)

% valores de busca para a direção de ondas
int_dir = [0,45,90,135,180,225,270,315];
geracao = [];

% taxa de variação para altura de ondas dx
%dx=0.1:0.01:0.5;
dx = 0.01:0.001:0.1;
% taxa de variação para parâmetros de profundidade e energia dxx
dxx=0.01:0.001:0.1;
%p_cross=randi([1,8]);
while size(selecao,1) ~= 0
    % sorteia dois indivíduos diferentes
    pos1 = randi(size(selecao,1));
    pai1 = selecao(pos1,1:end);
    pos2 = pos1;
    while (pos1 == pos2)
        pos2 = randi(size(selecao,1));
        pai2 = selecao(pos2,1:end);
    end
    % remove os pais do pool de seleção
    selecao([pos1, pos2],:) = [];
    
    % sorteia um valor entre 0.0 e 1.0 de uma distribuição uniforme
    p = rand;
    
    if p < probabilityOfCrossOver
        
        filho1 = pai1;
        filho2 = pai2;
        % faz o cruzamento, gerando dois descendentes
        filho1(1:p_cross) = pai2(1:p_cross);
        filho2(1:p_cross) = pai1(1:p_cross);
        
        % faz a mutação tanto nos pais quanto nos filhos
        mutacao = [pai1; pai2; filho1; filho2];
        
        % probabilidade de mutação (sem controle do usuário)
        if rand <0.9
            
            mutacao(1,1) = int_dir(randi(length(int_dir), 1));
            
            for i=1:size(mutacao,1)
                if rand>0.5
                    mutacao(i,2) = mutacao(i,2)+randi([-1,1])*dx(randi(length(dx), 1))*mutacao(i,2);
                    mutacao(i,3:end) = mutacao(i,3:end)+randi([-1,1])*dxx(randi(length(dxx), 1))*mutacao(i,3:end);
                else
                    mutacao(i,2) = mutacao(i,2)-randi([-1,1])*dx(randi(length(dx), 1))*mutacao(i,2);
                    mutacao(i,3:end) = mutacao(i,3:end)-randi([-1,1])*dxx(randi(length(dxx), 1))*mutacao(i,3:end);
                end
            end
            
            geracao = vertcat(geracao,mutacao);
        else
            geracao = vertcat(geracao,mutacao);
        end
    else
        mutacao = [pai1; pai2];
        
        % probabilidade de mutação (sem controle do usuário)
        if rand <0.9
            
            mutacao(1,1) = int_dir(randi(length(int_dir), 1));
            
            for i=1:size(mutacao,1)
                if rand>0.5
                    mutacao(i,2) = mutacao(i,2)+dx(randi(length(dx), 1))*mutacao(i,2);
                    mutacao(i,3:end) = mutacao(i,3:end)+dxx(randi(length(dxx), 1))*mutacao(i,3:end);
                else
                    mutacao(i,2) = mutacao(i,2)-dx(randi(length(dx), 1))*mutacao(i,2);
                    mutacao(i,3:end) = mutacao(i,3:end)-dxx(randi(length(dxx), 1))*mutacao(i,3:end);
                end
            end
            
            geracao = vertcat(geracao,mutacao);
        else
            geracao = vertcat(geracao,mutacao);
        end
    end
end

end