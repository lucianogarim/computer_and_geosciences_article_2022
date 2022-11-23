function [selecao] = selection(generation, n_sel)

caixa = generation;
selecao = [];

for i = 1:n_sel
    % sorteia dois indiv�duos diferentes para o torneio
    pos1 = randi(size(caixa,1));
    candidato1 = caixa(pos1,1:end);
    pos2 = pos1;
    while (pos2 == pos1)
         pos2 = randi(size(caixa,1));
         candidato2 = caixa(pos2,1:end);
    end
    % adiciona os indiv�duos selecionados ao torneio
    torneio = [[pos1 candidato1]; [pos2 candidato2]];
    % ordena os indiv�duos na lista do torneio por ordem decrescente de
    % fitness
    torneio = sortrows(torneio,size(caixa,2)+1,'descend');
    % adiciona o vencedor ao pool de sele��o
    selecao = vertcat(selecao, torneio(1,2:end-1));
    caixa(torneio(1,1),:) = [];
end