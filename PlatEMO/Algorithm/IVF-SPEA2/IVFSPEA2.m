classdef IVFSPEA2 < ALGORITHM
    % <2024> <multi> <real/integer/label/binary/permutation>
    % Strength Pareto evolutionary algorithm 2
    % C  --- 0.11 ---  Coleta
    % R  --- 0.1  ---  Ratio
    % M  --- 0.5  ---  Quantidade de mães a serem mutadas
    % V  --- 0.5  ---  Quantidade de variáveis de decições a ser mutadas
    % Cycles  --- 3 ---  Máximo de Ciclos
    % S  --- 1 ---  Steady State
    % N_Offspring  --- 1 ---  Número de Filhos gerados pelo SBX
    % EARN  --- 0 ---  0 Caso seja outro EAR, 1 para ser EARN
    % N_Obj_Limit --- 0 --- Limite de posição para assumir Current Father, 0 ativa seleção normal

    %------------------------------- Reference --------------------------------
    % E. Zitzler, M. Laumanns, and L. Thiele, SPEA2: Improving the strength
    % Pareto evolutionary algorithm, Proceedings of the Conference on
    % Evolutionary Methods for Design, Optimization and Control with
    % Applications to Industrial Problems, 2001, 95-100.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [C,R, M, V, Cycles, S, N_Offspring, EARN, N_Obj_Limit] = Algorithm.ParameterSet(13,0.3,0.5,0.5, 2, 1, 1, 0, 0);

            %% Geração da população inicial aleatória
            Population = Problem.Initialization();
            [Fitness, Forca, Distancia]    = CalFitness(Population.objs);
            
            Mating_N = Problem.N;

            IVF_Gen_FE = 0; %%Controlar quantidade da geração de evaluations do IVF (FE=Function Evaluation)
            IVF_Total_FE = 0; %%Controlar quantidade da Total de evaluations do IVF (FE=Function Evaluation)
            SPEA2_Gen = 1

            %%Definir tamanho o tamanho da coleta a partiar do tamanho da população
            Popultation_Size = length(Population); %%Tamanho do vetor de população

            %%Calcula o tamanho da coleta como um valor inteiro
            % ivf_c_size_p_result = floor(Popultation_Size * ivf_c_size_p);
            % ivf_collect_size = max(ivf_c_size, ivf_c_size_p_result)

            %% Otimização
            while Algorithm.NotTerminated(Population)
                Parameter = {1, 20, 1, 20};
                
                [Population, Zmin, IVF_Gen_FE, IVF_Total_FE, Mating_N] = IVF(Problem, Population, Fitness, Forca, Distancia, R, C, M, V, Cycles, IVF_Total_FE, N_Offspring, EARN, SPEA2_Gen, N_Obj_Limit);

                fprintf("Número de FE gastos durante IVF:\n");
                IVF_Gen_FE

                fprintf("Número de FE gastos durante IVF ao total:\n");
                IVF_Total_FE
                
                fprintf("Tamanho da População junto aos SuperBebes:\n");
                length(Population)

                [Fitness,~,~]    = CalFitness(Population.objs);
                

                % Gerar Offspring2 usando Tournament Selection e OperatorGA
                MatingPool = TournamentSelection(2, Mating_N, Fitness);
                
                fprintf("Tamanho da MatingPool:\n");
                length(MatingPool)

                
                Offspring  = OperatorGA(Problem, Population(MatingPool));
                
                [Population, Fitness] = EnvironmentalSelection([Population, Offspring], Problem.N);

                % --- BEGIN ADDED CODE: Save generation objectives --- 
                output_filename = 'generation_objectives.csv';
                num_objectives = Problem.M; % Get number of objectives

                % Check if file exists and write header if needed
                if SPEA2_Gen == 1 || ~isfile(output_filename)
                    fileID = fopen(output_filename, 'w');
                    header = 'Generation,Individual';
                    for obj_idx = 1:num_objectives
                        header = [header, sprintf(',Objective%d', obj_idx)];
                    end
                    fprintf(fileID, '%s\n', header);
                    fclose(fileID);
                end

                % Append data for the current generation
                fileID = fopen(output_filename, 'a');
                for i = 1:length(Population)
                    fprintf(fileID, '%d,%d', SPEA2_Gen, i);
                    fprintf(fileID, ',%.6f', Population(i).obj);
                    fprintf(fileID, '\n');
                end
                fclose(fileID);
                % --- END ADDED CODE --- 

                fprintf("Tamanho da População depois da sobrevivência:\n");
                length(Population)

                fprintf("Gasto total de FE:\n");
                Problem.FE

                SPEA2_Gen = SPEA2_Gen + 1;
            end
        end
    end
end