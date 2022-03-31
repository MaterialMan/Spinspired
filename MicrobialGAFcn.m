
%% Evolve substrate for a specific task
% This script can be used to evolve any reservoir directly to a task. It
% uses the steady-state Microbial Genetic Algorithm to evolve the best
% solution.

% Author: M. Dale
% Date: 03/07/19
function population = MicrobialGAFcn(ID,res_type,node_size,task,pop_size,num_gens,hetero_exchange_interaction)

%delete(gcp('nocreate'))
%parpool('local');

addpath(genpath('/mnt/lustre/users/md596/Magnets/multilayer'));

rng(ID,'twister');

%% Setup
config.test = ID;
config.parallel = 0;                        % use parallel toolbox

% type of network to evolve
config.res_type = res_type;             % state type of reservoir(s) to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network with multiple functions), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = node_size;                   % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each.
config = selectReservoirType(config);         % collect function pointers for the selected reservoir type

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = pop_size;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)
config.total_gens = num_gens;                    % number of generations to evolve
config.mut_rate = 0.01;                       % mutation rate
config.deme_percent = 0.2;                   % speciation percentage; determines interbreeding distance on a ring.
config.deme = round(config.pop_size*config.deme_percent);
config.rec_rate = 0.5;                       % recombination rate
config.error_to_check = 'train&val&test';    % printed error includes all three errors

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = task;          % Task to evolve for

config.hetero_exchange_interaction = hetero_exchange_interaction;

% get any additional params. This might include:
% details on reservoir structure, extra task variables, etc.
config = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% general params
config.gen_print = 10;                       % after 'gen_print' generations print task performance and show any plots
config.start_time = datestr(now, 'HH:MM:SS');
config.save_gen = config.gen_print;                       % save data at generation = save_gen
config.figure_array = [figure figure];

% Only necessary if wanting to parallelise the microGA algorithm
config.multi_offspring = 0;                 % multiple tournament selection and offspring in one cycle
config.num_sync_offspring = 0;%config.deme;    % length of cycle/synchronisation step

% type of metrics to apply; if necessary
config.metrics = {'linearMC'};          % list metrics to apply in cell array: see getVirtualMetrics.m for types of metrics available
config.record_metrics = 1;                  % save metrics

% additional pruning, if required
config.prune = 0;                           % after best individual is found prune the network to make more efficient
config.tolerance = [0 0 0];
config.prune_iterations = 500;

%% Run experiments
for test = 1:config.num_tests
    
    clearvars -except config test best best_indv store_error population ID node_size task
    
    warning('off','all')
    fprintf('\n Test: %d  ',test);
    fprintf('Processing genotype......... %s, task is %s \n',datestr(now, 'HH:MM:SS'),task)
    tic
    
    % tracker for sim plot
    last_best = 1;
    
    %set random seed
    rng(ID + test,'twister');
    
    % create initial population
    population = config.createFcn(config);
    
    %Assess population
%     for pop_indx = 1:config.pop_size
%         tic
%         population(pop_indx).core_indx = pop_indx;
%         population(pop_indx) = config.testFcn(population(pop_indx),config);
%         fprintf('\n i = %d, error = %.4f, took: %.4f\n',pop_indx,getError(config.error_to_check,population(pop_indx)),toc);
%     end
%     
%     
%     % find and print best individual
%     [best(test,1),best_indv(test,1)] = min(getError(config.error_to_check,population));
%     fprintf('\n Starting loop... Best error = %.4f\n',best(test,1));
%     
%     % store error that will be used as fitness in the GA
%     store_error(test,1,:) = getError(config.error_to_check,population);
%     
    %% start GA
    for gen = 2:config.total_gens
        
        % redefine seed - some functions/scripts may reset the seed
        rng(ID + test + gen,'twister');
        
        % reshape stored error to compare errors
        cmp_error = reshape(store_error(test,gen-1,:),1,size(store_error,3));

        % Tournment selection - pick two individuals
        equal = 1;
        while(equal) % find pair who are within deme range
            indv1 = randi([1 config.pop_size]);
            indv2 = indv1 + randi([1 config.deme]);
            if indv2 > config.pop_size
                indv2 = indv2 - config.pop_size; %loop around population ring if too big
            end
            if indv1 ~= indv2
                equal = 0;
            end
        end
        
        % Assess fitness of both and assign winner/loser
        if cmp_error(indv1) < cmp_error(indv2)
            winner=indv1; loser = indv2;
        else
            winner=indv2; loser = indv1;
        end
        
        % Infection and mutation. Place offspring in loser position
        population(loser) = config.recFcn(population(winner),population(loser),config);
        population(loser) = config.mutFcn(population(loser),config);
        
        %% Evaluate and update fitness
        [population(loser)] = config.testFcn(population(loser),config);
        
        %update errors
        store_error(test,gen,:) =  store_error(test,gen-1,:);
        store_error(test,gen,loser) = getError(config.error_to_check,population(loser));%population(loser).val_error;
        [best(test,gen),best_indv(test,gen)] = min(store_error(test,gen,:));
        
        % print info
        if (mod(gen,config.gen_print) == 0)
           fprintf('Gen %d, time taken: %.4f sec(s)\n  Winner: %.4f, Loser: %.4f, Best Error: %.4f \n',gen,toc/config.gen_print,getError(config.error_to_check,population(winner)),getError(config.error_to_check,population(loser)),getError('test',population(best_indv(test,gen))));
            tic;
            if best_indv(gen) ~= last_best
                % plot reservoir structure, task simulations etc.
                %plotReservoirDetails(population,best_indv(test,:),gen,loser,config);
            end
            
            last_best = best_indv(gen);
        end
        
        
        %save data
        if mod(gen,config.save_gen) == 0
            saveData(population,test,store_error,config)
            %saveTable(population,store_error,config)
        end
        
        % if found best, stop
        if best(test,gen) == 0
            return;
        end
    end
    
    % apply metrics to final population
    if config.record_metrics
        parfor pop_indx = 1:config.pop_size
            population(pop_indx).core_indx = config.test + pop_indx;
            metrics(pop_indx,:) = getMetrics(population(pop_indx),config);
            population(pop_indx).metrics = metrics(pop_indx,:);
        end
    end
    
    if config.prune
        rng(1,'twister')
        pop_error = [population.train_error]+ [population.val_error]+[population.test_error];
        
        [~,indx] = sort(pop_error);
        
        parfor p = 1:4
            [Pruned_best_individual{p},old_W_fitness(p),old_Win_fitness(p),full_W{p}] = maxPruning(@testReservoir,{'train_error','val_error','test_error'},population(indx(p))...
                ,[population(indx(p)).train_error population(indx(p)).val_error population(indx(p)).test_error]...
                ,config.tolerance,config);
        end
        % plotReservoirDetails(population,best_indv,gen,loser,config)
    end
end
end

function saveData(population,test,store_error,config)
config.figure_array =[];
if iscell(config.res_type)
    res_type ='Heterotic';
else
    res_type = config.res_type;
end

run = num2str(config.test);
run = run(end-1:end);

save(strcat(res_type,'_',config.dataset,'_run',run,'_gens',num2str(config.total_gens),'_',num2str(config.total_units),'nodes_het',num2str(config.hetero_exchange_interaction),'_graph_',config.graph_type,'.mat'),...
    'population','store_error','config','-v7.3');
end

function saveDataTable(database,tests,material,config)
% config.figure_array =[];
% save(strcat('Framework_substrate_',material,'_',config.res_type,'_run',num2str(tests),'_gens',num2str(config.total_gens),'_',num2str(config.num_reservoirs),'Nres_'),...
%     'database_history','database','config','quality','-v7.3');

T_database = database;
metrics = reshape([database.behaviours],length(database(1).behaviours),length(database))';

% remove unwanted fields
fields = {'xy','input_weights','input_widths','last_state'};
T_database = rmfield(T_database,fields);

T = repmat({material},1,length(database)); %struct('material',
T = num2cell(T);
[T_database.material]=  T{:};

% metrics
T = num2cell(metrics(:,1));
[T_database.KR]=  T{:};
T = num2cell(metrics(:,3));
[T_database.GR]=  T{:};
T = num2cell(metrics(:,2));
[T_database.MC]=  T{:};                                                                        % metrics
                                                                                                                                                 
Table = struct2table(T_database);     
writetable(Table,strcat('task_',material,'_run',num2str(tests),'_size',num2str(config.num_nodes{1}),'.csv'),'Delimiter',',','QuoteStrings',true);

end
