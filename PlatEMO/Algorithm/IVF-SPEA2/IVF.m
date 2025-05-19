% Função principal: IVF
%function Offspring = IVF(Problem, Population, C, M, Parameter, Fitness)
function [Population, Zmin, IVF_Gen_FE, IVF_Total_FE, Mating_N] = IVF(Problem, Population, Fitness, Forca, Distancia, ivf_rate, C, M, V, Cycles, IVF_Total_FE, N_Offspring, EARN, SPEA2_Gen, N_Obj_Limit)
    % IVF realiza a seleção e geração de descendentes através de mutação e crossover.
    %
    %   Offspring = IVF(Problem, Population, C, Parameter, Fitness)
    %
    %   Problem - Instância do problema, usada para avaliação e limites.
    %   Population - População atual (array de objetos SOLUTION).
    %   C - Percentual de indivíduos a serem selecionados na coleta.
    %   Parameter - Parâmetros para as operações de crossover e mutação (proC, disC, proM, disM).
    %   Fitness - Vetor de fitness correspondente à população.

    %Guardar a quantidade de avaliações realizadas na run para controlar o total consumido do hospedeiro (FE) e não estourar o total
    fprintf("\n\n\n\n\n ================================================================================= Entrada no IVF =================================================================================\n\n\n\n\n ") 
    
    FE_Before_IVF = Problem.FE;
    Zmin = 0;

    Salvar_Arquivo_CSV = false;

    %Gatilho de execução
    %Verificar o IVF deve ou não ser executado
    
    %Gatilho Sávio

    %if IVF_Total_FE > ivf_rate * Problem.FE
        %Atualizar total de avaliações consumidar e outras informações que são retornadas pela função IVF
        %IVF_Gen_FE = 0;    Como o IVF não foi executado não houve consumo de FE nessa gen
        %IVF_Total_FE = IVF_Total_FE + IVF_Gen_FE;  %Total de FE geradas pelo InVitro
        %Mating_N = Problem.N   %O algoritmo infitrião deverá gerar a quantidade padrão de filhos
        %Esses valores acima serão retornados ao algoritmo anfitrião, juntos ao Population e Zmin
        %return;    %Sair do IVF
    %end;

    fprintf("Avaliações utilizadas até agora pelo Problema no Geral")
    Problem.FE
    fprintf("Avaliações utilizadas até agora pelo IVF")
    IVF_Total_FE

    if IVF_Total_FE > ivf_rate * Problem.FE
        %Atualizar total de avaliações consumidar e outras informações que são retornadas pela função IVF
        IVF_Gen_FE = 0;    %Como o IVF não foi executado não houve consumo de FE nessa gen
        IVF_Total_FE = IVF_Total_FE + IVF_Gen_FE;  %Total de FE geradas pelo InVitro
        Mating_N = Problem.N   %O algoritmo infitrião deverá gerar a quantidade padrão de filhos
        %Esses valores acima serão retornados ao algoritmo anfitrião, juntos ao Population e Zmin
        fprintf("Nessa geração não será utilizado o IVF para contribuição genética. Pulando...")
        return;    %Sair do IVF
    end;

    %Caso não tenha ocorrido o return o IVF será executado.
    

    %População para armazenar todos os filhos IVF gerados
    All_IVF_Offspring = [];

    %ETAPAS INVITRO

    % ================= Inicio Coleta ===================
    
    % fprintf('-----------------------------\n');
    % fprintf('Inicio Coleta\n');
    IVF_Population = [];

    % Ordena os valores de aptidão em ordem crescente e obtém os índices ordenados
    [SortedFitness, SortedIndices] = sort(Fitness, 'ascend');

    % Determina o número de pais, garantindo que seja pelo menos Problem.M + 1
    NumParents = max(round(Problem.N * C), Problem.M + 1);
    NumParents = min(NumParents, Problem.N); % Garante que não exceda o tamanho da população

    % Define o intervalo máximo para seleção do pai
    Max_Range = min(NumParents * 2, Problem.N); % Garante que Max_Range não exceda a população

    % Seleciona aleatoriamente o índice do pai dentro do intervalo permitido
    Pos_Pai=randi([1, Max_Range]);
    Father_Index = SortedIndices(Pos_Pai);
    Father = Population(Father_Index);
    FatherFitness = Fitness(Father_Index);

    % Seleciona os NumParents melhores indivíduos
    ParentsIndices = SortedIndices(1:NumParents);

    % Garante que as mães sejam os NumParents - 1 melhores indivíduos diferentes do pai
    if ismember(Father_Index, ParentsIndices)
        % Pai está dentro dos candidatos - apenas removê-lo
        MothersIndices = ParentsIndices(ParentsIndices ~= Father_Index);
    else
        % Pai está fora dos candidatos - pega os NumParents - 1 melhores
        MothersIndices = ParentsIndices(1:NumParents - 1);
    end

    % Coleta as mães e seus fitness
    Mothers = Population(MothersIndices);
    MothersFitness = Fitness(MothersIndices);

    % length(Father)
    % length(Mothers)

    % fprintf(' ===================== Fim Coleta ======================\n\n\n')

    % ================= Fim Coleta ===================

    % ================= Início EAR ===================

    %=====Mutando as Mães=========

    if M == 0
        % 🔹 Caso seja AR (Sem mutação das mães)
        % fprintf('É AR!!!!!!!!\n');
        MutatedMothers = Mothers;
        IVF_Population = [IVF_Population, Father, MutatedMothers];
    
    else
        if EARN == 1
            % 🔹 Caso seja EARN: Gerar novas soluções aleatoriamente no lugar das mães
            % fprintf('É EARN!!!!!!!!\n');
            % fprintf('Mothers antes de EARN:\n')
            % Mothers.decs

            MutatedMothers_decs = Problem.Initialization_with_no_evaluation(NumParents-1);
            Mothers.setDec(MutatedMothers_decs);
            % fprintf('Mothers depois de EARN:\n')
            % Mothers.decs
    
        else
            % 🔹 Caso contrário, aplica mutação polinomial normalmente (EAR)
            % fprintf('É EAR!!!!!!!!\n');
            [MutatedMothers_decs, MutateMothersIdx] = PolynomialMutation(Mothers, M, V, Problem.lower, Problem.upper, 1, 20);
            Mothers.setDec(MutatedMothers_decs);
            for i = 1:length(Mothers)
                Mothers(i).ivf = true;
                if ismember(i, MutateMothersIdx)
                    Mothers(i).mae_mutada = true;
                else
                    Mothers(i).mae = true;
                end
            end
        end
    
        % 🔹 Criar novas soluções a partir dos valores gerados
        MutatedMothers = Mothers;
    
        % 🔹 Atualizar a população IVF
        IVF_Population = [IVF_Population, Father, MutatedMothers];
    end

    %=====Fim das mutações=========
    % fprintf("Decs do Pai:");
    % Father.decs
    % fprintf("Decs das Mães:");
    % MutatedMothers.decs
    
    % fprintf("==========================");
    % fprintf("IVF Population:");
    % IVF_Population
    % fprintf("==========================\n");


    %======Reprodução Assistida======
    % fprintf("==========================");
    % fprintf("Entrada na Reprodução Assistida");
    % fprintf("==========================\n");

    %Definindo o limite máximo de avaliações por execução
    Limite_Maximo_Avals = Problem.N;
    Avaliacoes_Por_Ciclo = length(MutatedMothers);
    IVF_Gen_FE = Problem.FE - FE_Before_IVF;
    Current_Father = Father;
    % fprintf("Será observado o objetivo: \n");
    % obj_alvo = randi(size(Current_Father.obj,2)); % Sorteia um índice válido do número de objetivos
    % obj_alvo

    % Cycles

    %=================================== Criando Arquivo ==================================%
    if Salvar_Arquivo_CSV
        csv_filename = 'historico_populacao.csv';
    
        % Criar cabeçalho do CSV se for a primeira execução
        num_objetivos = length(Population(1).obj);
        if ~isfile(csv_filename)
            fileID = fopen(csv_filename, 'w'); % 'w' para criar um novo arquivo
            fprintf(fileID, 'Geracao_Anfitriao,SubGeracao_IVF,Individuo,');
            
            % Criar cabeçalhos para os objetivos (ex: Objetivo_1, Objetivo_2, ...)
            for obj_idx = 1:num_objetivos
                fprintf(fileID, 'Objetivo_%d,', obj_idx);
            end
            
            % Cabeçalhos finais
            fprintf(fileID, 'Tag\n');
            fclose(fileID);
        end
    end

    PopComparacao = Population;

    for IVF_Cycle=1:1:Cycles
    %while IVF_Cycle <= IVF_Max_Cycles
        IVF_Cycle
        
        %Verifica se adicionar mais avaliações excederá o limite
        if IVF_Gen_FE + Avaliacoes_Por_Ciclo > Limite_Maximo_Avals
            fprintf("===========================\n Atingiu o máximo de Limite de FEs\n======================\n")
            break;
        end

        %Realiza os cruzamento do pai com todas as mães
        % Current_Father
        % MutatedMothers 
        
        IVF_Offspring = IVF_Recombination(Current_Father, MutatedMothers, Cycles, Problem, N_Offspring);

        IVF_Offspring = Problem.Evaluation(IVF_Offspring);

        for i = 1:length(IVF_Offspring)
            % i
            IVF_Offspring(i).ivf = true;
            if MutatedMothers(i).mae_mutada
                IVF_Offspring(i).filho_mae_mutada = true;
            else
                IVF_Offspring(i).filho = true;
            end
            % MutatedMothers(i)
            % IVF_Offspring(i)
        end

        % fprintf("Decs dos Filhos:");
        % IVF_Offspring.decs
        % IVF_Offspring

        % IVF_Gen_FE = Problem.FE - FE_Before_IVF;
        % IVF_Gen_FE

        
        % IVF_Offspring

        IVF_Gen_FE = Problem.FE - FE_Before_IVF;
        % IVF_Gen_FE

        % fprintf("==========================");
        % fprintf("IVF Offspring:");
        % IVF_Offspring
        % fprintf("==========================");

        All_IVF_Offspring = [All_IVF_Offspring, IVF_Offspring];

        % All_IVF_Offspring

        % fprintf("==========================\n");
        % fprintf("Comparando pai e filhos gerados:\n");
        
        Current_Father.ivf = true;
        Current_Father.pai = true;

        % for i = 1:length(IVF_Offspring)
        %     IVF_Offspring(i).ivf = true;
        %     IVF_Offspring(i).filho = true;
        % end

        PopComparacao = [PopComparacao, IVF_Offspring];

        if Salvar_Arquivo_CSV
            % **Salvar no CSV antes da seleção ambiental**
            fileID = fopen(csv_filename, 'a'); % Abrir para adicionar novas linhas
            for i = 1:length(PopComparacao)
                % Escrever geração e sub-geração
                fprintf(fileID, '%d,%d,', SPEA2_Gen, IVF_Cycle);
                
                % Escrever vetor de decisões (concatenado em string)
                decs_str = sprintf('%.4f ', PopComparacao(i).decs);
                fprintf(fileID, '"%s",', strtrim(decs_str)); % Usando aspas para evitar erros com CSV
                
                % Escrever os valores de objetivo
                for obj_idx = 1:num_objetivos
                    fprintf(fileID, '%.4f,', PopComparacao(i).obj(obj_idx));
                end
                
                % Determinar a tag do indivíduo
                if ~isempty(PopComparacao(i).pai) && PopComparacao(i).pai
                    tag = 'pai';
                elseif ~isempty(PopComparacao(i).mae) && PopComparacao(i).mae
                    tag = 'mae';
                elseif ~isempty(PopComparacao(i).filho) && PopComparacao(i).filho
                    tag = 'filho';
                elseif ~isempty(PopComparacao(i).filho_mae_mutada) && PopComparacao(i).filho_mae_mutada
                    tag = 'filho_mae_mutada';
                else
                    tag = 'anfitriao';
                end
                
                % Escrever a tag no CSV
                fprintf(fileID, '%s\n', tag);
            end

            fclose(fileID); % Fechar arquivo após cada ciclo]
        end

        
        [PopComparacao, PopComparacaoFitness] = EnvironmentalSelection(PopComparacao, Problem.N);

        if N_Obj_Limit == 0
            % Selecionar apenas o pai (indivíduo com pai == 1)
            indices_pai = find(cellfun(@(x) ~isempty(x) && x, {PopComparacao.pai}));

            % indices_pais

            % Exibir valores de fitness dos pais sobreviventes
            if ~isempty(indices_pai)
                % disp('Pais sobreviventes e seus valores de fitness:');
                Fitness_Father = PopComparacaoFitness(indices_pai)
                for i = 1:length(indices_pai)
                    idx = indices_pai(i);
                    % fprintf('Pai %d: Fitness = %.4f\n', i, PopComparacaoFitness(idx));
                    % PopComparacao(idx)
                end
            else
                Fitness_Father = 100000
                % disp('Nenhum pai sobreviveu.');
            end

            % Filtrar índices dos filhos sobreviventes
            indices_filhos = find(cellfun(@(x) ~isempty(x) && x, {PopComparacao.filho}));
            % indices_filhos

            % Exibir valores de fitness dos filhos sobreviventes
            if ~isempty(indices_filhos)
                % disp('Filhos sobreviventes e seus valores de fitness:');
                Fitness_Offspring = PopComparacaoFitness(indices_filhos);
                for i = 1:length(indices_filhos)
                    idx = indices_filhos(i);
                    % fprintf('Filho %d: Fitness = %.4f\n', i, PopComparacaoFitness(idx));
                    % PopComparacao(idx)
                end
            else
                % disp('Nenhum filho sobreviveu.');
                break;
            end

            % Exibir os resultados
            % disp(['Fitness do Pai: ', num2str(Fitness_Father)]);
            % disp(['Fitness dos Filhos: ', num2str(Fitness_Offspring)]);

            % Encontrar o índice do filho com menor fitness
            [min_fitness, melhor_idx_filho] = min(Fitness_Offspring);

            % Ajustar índice para referenciar corretamente `IVF_Offspring`
            idx_melhor_filho = indices_filhos(melhor_idx_filho);

            for i = 1:length(PopComparacao)
                PopComparacao(i).pai = [];
                PopComparacao(i).filho = [];
                PopComparacao(i).filho_mae_mutada = [];
            end
        

            % Verificar se o melhor filho tem fitness menor que o pai
            if min_fitness < Fitness_Father
                % disp('O melhor filho tem fitness menor que o pai:');
                % fprintf('Filho %d com Fitness = %.4f, Pai com Fitness = %.4f\n', idx_melhor_filho, min_fitness, Fitness_Father);
                % disp(IVF_Offspring(idx_melhor_filho).decs()); % Mostrar as variáveis de decisão

                % fprintf("Valores OBJ do Pai:\n");
                disp(Current_Father.obj);

                % fprintf("Valores OBJ do Filho:\n");
                disp(PopComparacao(idx_melhor_filho).obj);

                % Atualizar o pai para o melhor filho
                Current_Father = PopComparacao(idx_melhor_filho);
                Fitness_Father = min_fitness; % Atualizar fitness do pai
            else
                % disp('Nenhum filho tem fitness menor que o pai.');
                break;
            end
        else
            % disp('=================== Seleção por N_Obj_Limit ====================');
            % Ordenar população por fitness
            [~, sortedIndices] = sort(PopComparacaoFitness);
            sortedPop = PopComparacao(sortedIndices);


            indices_pai = find(cellfun(@(x) ~isempty(x) && x, {sortedPop.pai}));
            Fitness_Father = PopComparacaoFitness(indices_pai);

            indices_filhos = find(cellfun(@(x) ~isempty(x) && x, {sortedPop.filho}));

            % Determinar o intervalo de seleção
            selectionRange = 1:min(Max_Range, length(sortedPop));
            % selectionRange
        
            % Inicializar flag para novo pai
            newFatherSelected = false;
        
            % Iterar sobre os melhores indivíduos para encontrar um filho
            for i = selectionRange
                % i
                idx = sortedIndices(i); % Recuperar índice original
                % idx
                if any(PopComparacao(idx).filho) || any(PopComparacao(idx).filho_mae_mutada)
                    % Filho encontrado dentro do limite de seleção
                    fprintf('Novo Pai selecionado: Filho na posição %d com Fitness = %.4f substituira Pai Fitness = %.4f, era o %d melhor individuo\n', idx, PopComparacaoFitness(idx), Fitness_Father, i);
                    Current_Father = PopComparacao(idx);
                    Fitness_Father = PopComparacaoFitness(idx);
                    newFatherSelected = true;                
                    for j = 1:length(PopComparacao)
                        PopComparacao(j).pai = [];
                        PopComparacao(j).filho = [];
                        PopComparacao(j).filho_mae_mutada = [];
                    end
                    break;
                end
            end

            for i = 1:length(PopComparacao)
                PopComparacao(i).pai = [];
                PopComparacao(i).filho = [];;
                PopComparacao(i).filho_mae_mutada = [];
            end
            
            if ~newFatherSelected
                % disp('Nenhum filho dentro do limite de seleção.');
                % disp('=================== Fim da Seleção por N_Obj_Limit ====================');
                break;
            end
            % disp('=================== Fim da Seleção por N_Obj_Limit ====================')
        end

    end

    Population = PopComparacao;

    for i = 1:length(Population)
        Population(i).ivf = true;
        Population(i).filho = [];
        Population(i).pai = [];
        Population(i).mae = [];
    end

    Mating_N = max(Problem.N - IVF_Gen_FE, 2);
    
    % fprintf("==========================");
    % fprintf("Saida da Reprodução Assistida");
    % fprintf("==========================");
    
    %======Fim da Reprodução Assistida======

    Zmin = 123412;
    IVF_Total_FE = IVF_Total_FE + IVF_Gen_FE;


    fprintf("\n\n\n\n\n ================================================================================= Saida do IVF =================================================================================\n\n\n\n\n ") 
% ================= Fim EAR ===================
end


function [MutatedMothers, MutateMothersIdx] = PolynomialMutation(Mothers, M, V, lower, upper, proM, disM)
    % PolynomialMutation - Aplica mutação polinomial em um subconjunto de mães e variáveis.
    %
    %   M - Fração de mães que sofrerão mutação (ex: 0.5 = 50% das mães).
    %   V - Fração de variáveis mutadas dentro de cada mãe escolhida (ex: 0.3 = 30% dos genes).
    
    if isa(Mothers, 'SOLUTION')
        Mothers_decs = Mothers.decs;  % Assume que 'decs' contém os valores das variáveis de decisão
    end

    [N,D] = size(Mothers_decs); % Número de indivíduos (N) e variáveis (D)
    MutatedMothers = Mothers_decs;  % Inicializa os filhos com os pais

    % 🔹 Escolher quais mães serão mutadas com probabilidade M
    NumMothersToMutate = round(M * N);  % Número de mães a serem mutadas

    % Selecionar as últimas 'NumMothersToMutate' mães
    MutateMothersIdx = (N - NumMothersToMutate + 1):N;

    % 🔹 Escolher quais variáveis dentro dessas mães serão mutadas com probabilidade V
    Site = false(N,D);  % Matriz booleana inicial
    for i = MutateMothersIdx
        NumVarsToMutate = round(V * D);  % Número de variáveis a mutar
        VarIndices = randperm(D, NumVarsToMutate);  % Escolhe as variáveis aleatórias
        Site(i, VarIndices) = true;  % Marca essas variáveis para mutação
    end

    mu = rand(N,D);  % Números aleatórios para a mutação
    Lower = repmat(lower, N, 1);
    Upper = repmat(upper, N, 1);

    % 🔹 Garantir que os valores permaneçam dentro dos limites
    MutatedMothers = min(max(MutatedMothers, Lower), Upper);

    % 🔹 Aplicar mutação para valores menores que 0.5
    temp = Site & mu <= 0.5;
    MutatedMothers(temp) = MutatedMothers(temp) + (Upper(temp) - Lower(temp)) .* ...
        ((2 .* mu(temp) + (1 - 2 .* mu(temp)) .* (1 - (MutatedMothers(temp) - Lower(temp)) ./ ...
        (Upper(temp) - Lower(temp))) .^ (disM + 1)) .^ (1 / (disM + 1)) - 1);

    % 🔹 Aplicar mutação para valores maiores que 0.5
    temp = Site & mu > 0.5;
    MutatedMothers(temp) = MutatedMothers(temp) + (Upper(temp) - Lower(temp)) .* ...
        (1 - (2 .* (1 - mu(temp)) + 2 .* (mu(temp) - 0.5) .* (1 - (Upper(temp) - MutatedMothers(temp)) ./ ...
        (Upper(temp) - Lower(temp))) .^ (disM + 1)) .^ (1 / (disM + 1)));

    % 🔹 Mostrar mães antes e depois da mutação
    % fprintf("\nMães Selecionadas para Mutação (%d de %d):\n", NumMothersToMutate, N);
    % disp(Mothers_decs(MutateMothersIdx, :));
    % fprintf("\nMães Após a Mutação:\n");
    % disp(MutatedMothers(MutateMothersIdx, :));
end

    


% ==================================================================================================


% Função EAR - IVF_Recombination
function Offspring = IVF_Recombination(Father, Mothers, Cycles, Problem, N_Offspring)
% IVF_Recombination realiza a recombinação de um indivíduo pai com vários indivíduos mães.
%
% Entradas:
%   - Father: Objeto SOLUTION representando o pai da recombinação.
%   - Mothers: Objeto SOLUTION representando os indivíduos mães.
%   - Cycles: (Não utilizado nesta implementação, mas pode ser útil em versões futuras).
%   - Problem: Estrutura contendo informações do problema (não utilizada diretamente aqui).
%
% Saída:
%   - Offspring: Matriz ((N-1) x D) contendo os novos indivíduos gerados pela recombinação.

    % 🔹 Parâmetros de cruzamento (100% de chance de cruzamento)
    proC = 1;   % Probabilidade de cruzamento (100%)
    disC = 20;  % Índice de distribuição SBX (define a intensidade da variação)

    % 🔹 Obtém os valores numéricos das variáveis de decisão do pai e das mães
    Father_dec  = Father.decs;   % Extrai a matriz de decisão do pai (1 x D)
    Mothers_dec = Mothers.decs;  % Extrai a matriz de decisão das mães (N x D)

    % Father_dec
    % Mothers_dec

    % 🔹 Número de mães (N) e número de variáveis de decisão (D)
    [N, D] = size(Mothers_dec);

    % [N, D]

    % 🔹 Repete o vetor do pai N vezes para combinar com cada mãe
    Fathers_dec = repmat(Father_dec, N, 1);

    % Fathers_dec

    % 🔹 Inicializa a matriz beta para controlar a variação gerada pela recombinação
    beta = zeros(N, D);

    % beta

    % 🔹 Gera valores aleatórios para cada variável e indivíduo
    mu = rand(N, D);

    % mu

    % 🔹 Calcula o beta para indivíduos com mu ≤ 0.5 (primeira equação do SBX)
    beta(mu <= 0.5) = (2 * mu(mu <= 0.5)).^(1 / (disC + 1));

    % 🔹 Calcula o beta para indivíduos com mu > 0.5 (segunda equação do SBX)
    beta(mu > 0.5) = (2 - 2 * mu(mu > 0.5)).^(-1 / (disC + 1));

    % 🔹 Aplica um fator aleatório (-1 ou 1) para diversificação
    beta = beta .* (-1).^randi([0,1], N, D);

    % 🔹 50% das variáveis não sofrem modificação (mantêm beta = 1)
    beta(rand(N, D) < 0.5) = 1;

    % 🔹 Como proC = 1, essa linha não altera nada, mas em outros casos, 
    %    ela impediria o cruzamento em alguns indivíduos
    beta(repmat(rand(N, 1) > proC, 1, D)) = 1;

    % 🔹 Gera dois filhos para cada cruzamento:
    %    - Primeiro filho: média dos pais + variação beta
    %    - Segundo filho: média dos pais - variação beta
    if N_Offspring == 1
        % Gerar um único filho aleatório diretamente
        if rand() < 0.5
            Offspring = (Fathers_dec + Mothers_dec) / 2 + beta .* (Fathers_dec - Mothers_dec) / 2;
        else
            Offspring = (Fathers_dec + Mothers_dec) / 2 - beta .* (Fathers_dec - Mothers_dec) / 2;
        end
    else
        Offspring = [(Fathers_dec + Mothers_dec) / 2 + beta .* (Fathers_dec - Mothers_dec) / 2;
                        (Fathers_dec + Mothers_dec) / 2 - beta .* (Fathers_dec - Mothers_dec) / 2];
    end

end
    

% Funções de Mutação

% Mutation para variáveis reais
function Offspring = GAMutationReal(Mother, lower, upper, proM, disM)
    Mother = Mother(:);
    lower = lower(:);
    upper = upper(:);

    N = length(Mother);
    Site = rand(N, 1) < proM / N;
    mu = rand(N, 1);
    Offspring = Mother;

    temp = Site & mu <= 0.5;
    Offspring(temp) = Mother(temp) + (upper(temp) - lower(temp)) .* ...
        ((2 .* mu(temp) + (1 - 2 .* mu(temp)) .* (1 - (Mother(temp) - lower(temp)) ./ ...
        (upper(temp) - lower(temp))).^(disM + 1)).^(1 / (disM + 1)) - 1);

    temp = Site & mu > 0.5;
    Offspring(temp) = Mother(temp) + (upper(temp) - lower(temp)) .* ...
        (1 - (2 .* (1 - mu(temp)) + 2 .* (mu(temp) - 0.5)) .* (1 - (upper(temp) - Mother(temp)) ./ ...
        (upper(temp) - lower(temp))).^(disM + 1)).^(1 / (disM + 1)) - 1;
end

% Mutation para variáveis de rótulo
function Offspring = GAMutationLabel(Mother, proM)
    N = length(Mother);
    Site = rand(1, N) < proM / N;
    Rand = randi([0, 1], 1, N);
    Offspring = Mother;
    Offspring(Site) = Rand(Site);
end

% Mutation para variáveis binárias
function Offspring = GAMutationBinary(Mother, proM)
    N = length(Mother);
    Site = rand(1, N) < proM / N;
    Offspring = Mother;
    Offspring(Site) = ~Mother(Site);
end

% Mutation para variáveis de permutação
function Offspring = GAMutationPermutation(Mother, proM)
    N = length(Mother);
    Offspring = Mother;
    numSwaps = round(proM * N / 100);
    for swap = 1:numSwaps
        idx = randperm(N, 2);
        Offspring([idx(1), idx(2)]) = Offspring([idx(2), idx(1)]);
    end
end

% Funções de Crossover

% Crossover para variáveis reais
function Offspring = GACrossoverReal(Father, Mother, lower, upper, proC, disC)
    N = length(Father);
    beta = rand(1, N);

    beta(beta <= 0.5) = (2 * beta(beta <= 0.5)).^(1 / (disC + 1));
    beta(beta > 0.5) = (2 - 2 * beta(beta > 0.5)).^(-1 / (disC + 1));

    beta(rand(1, N) < 0.5) = 1;

    Offspring1 = (Father + Mother) / 2 + beta .* (Father - Mother) / 2;
    Offspring2 = (Father + Mother) / 2 - beta .* (Father - Mother) / 2;

    Offspring = [Offspring1; Offspring2];

    if size(Offspring, 1) ~= 2
        error('GACrossoverReal deve retornar duas colunas (descendentes).');
    end
end

% Crossover para variáveis de rótulo
function Offspring = GACrossoverLabel(Father, Mother, proC)
    N = length(Father);
    k = rand(1, N) < 0.5;
    k(rand(1, N) > proC) = false;
    Offspring1 = Father;
    Offspring2 = Mother;
    Offspring1(k) = Mother(k);
    Offspring2(k) = Father(k);
    Offspring = [Offspring1; Offspring2];
end

% Crossover para variáveis binárias
function Offspring = GACrossoverBinary(Father, Mother, proC)
    N = length(Father);
    k = rand(1, N) < 0.5;
    k(rand(1, N) > proC) = false;
    Offspring1 = Father;
    Offspring2 = Mother;
    Offspring1(k) = Mother(k);
    Offspring2(k) = Father(k);
    Offspring = [Offspring1; Offspring2];
end

% Crossover para variáveis de permutação
function Offspring = GACrossoverPermutation(Father, Mother)
    N = length(Father);
    point1 = randi([1, N-1]);
    point2 = randi([point1+1, N]);

    child1 = zeros(1, N);
    child2 = zeros(1, N);

    child1(point1:point2) = Father(point1:point2);
    child2(point1:point2) = Mother(point1:point2);

    child1 = fill_child(child1, Mother, point2);
    child2 = fill_child(child2, Father, point2);

    Offspring = [child1; child2];
end

% Função auxiliar para Order Crossover (OX)
function child = fill_child(child, parent, point)
    N = length(child);
    current_pos = mod(point, N) + 1;
    parent_pos = mod(point, N) + 1;
    while any(child == 0)
        gene = parent(parent_pos);
        if ~ismember(gene, child)
            child(current_pos) = gene;
            current_pos = mod(current_pos, N) + 1;
        end
        parent_pos = mod(parent_pos, N) + 1;
    end
end
