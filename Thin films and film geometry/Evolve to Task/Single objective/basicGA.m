
%% Evolve substrate for a specific task
% This script can be used to evolve any reservoir directly to a task. It
% uses the steady-state Microbial Genetic Algorithm to evolve the best
% solution.

% Author: M. Dale
% Date: 03/07/19
% clear %vars -except population
% close all
function basicGA(ID,res_size)
% add all subfolders to the path --> make all functions in subdirectories available
addpath(genpath('/mnt/lustre/users/md596/Magnets/MassMM')); % make sure filepath is correct

rng(ID,'twister');

%% Setup
config.parallel = 0;                        % use parallel toolbox

% type of network to evolve
config.res_type = 'MM';            % state type of reservoir(s) to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network with multiple functions), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more. Place reservoirs in cell ({}) for heterotic systems.
config.num_nodes = [res_size];                   % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each.
config = selectReservoirType(config);         % collect function pointers for the selected reservoir type

%% Evolutionary parameters
config.num_tests = ID;                         % num of tests/runs

config.num_parents = 1;
config.num_offspring = 3;
config.pop_size = config.num_offspring + config.num_parents;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

config.total_gens = 2000;                    % number of generations to evolve
config.mut_rate = 0.1;                       % mutation rate
config.deme_percent = 0.1;                   % speciation percentage; determines interbreeding distance on a ring.
config.deme = round(config.pop_size*config.deme_percent);
config.rec_rate = 0.5;                       % recombination rate
config.error_to_check = 'train&val&test';

%% Task parameters
config.dataset = 'narma_10';          % Task to evolve for
%config.figure_array = [figure figure];

% get any additional params. This might include:
% details on reservoir structure, extra task variables, etc.
config = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% general params
config.gen_print = 1;                       % after 'gen_print' generations print task performance and show any plots
config.start_time = datestr(now, 'HH:MM:SS');
config.save_gen = 25;                       % save data at generation = save_gen

% type of metrics to apply; if necessary
config.metrics = {'KR','GR','linearMC'};          % list metrics to apply in cell array: see getVirtualMetrics.m for types of metrics available
config.record_metrics = 0;                  % save metrics

% additional pruning, if required
config.prune = 0;                           % after best individual is found prune the network to make more efficient
config.tolerance = [0 0 0];
config.prune_iterations = 500;

%% Run experiments
for test = ID%1:config.num_tests
    
    clearvars -except config ID test best best_indv store_error population
    
    warning('off','all')
    fprintf('\n Test: %d  ',test);
    fprintf('Processing genotype......... %s \n',datestr(now, 'HH:MM:SS'))
    tic
    
    % tracker for sim plot
    last_best = 1;
    config.test = test;
    
    %set random seed
    rng(test,'twister');
    
    % create initial population
    population = config.createFcn(config);
    
    %Assess population
    if config.parallel % use parallel toolbox - faster
        parfor pop_indx = 1:config.pop_size
            warning('off','all')
            population(pop_indx) = config.testFcn(population(pop_indx),config);
            fprintf('\n i = %d, error = %.4f\n',pop_indx,getError(config.error_to_check,population(pop_indx)));
        end
    else
        for pop_indx = 1:config.pop_size
            tic
            population(pop_indx) = config.testFcn(population(pop_indx),config);
            fprintf('\n i = %d, error = %.4f, took: %.4f',pop_indx,getError(config.error_to_check,population(pop_indx)),toc);
        end
    end
    
    % find and print best individual
    [best(test,1,:),best_indv(test,1,:)] = sort(getError(config.error_to_check,population));
    fprintf('\n Starting loop... Best error = %.4f\n',best(test,1));
    
    % store error that will be used as fitness in the GA
    store_error(test,1,:) = getError(config.error_to_check,population);
    
    % plot best
    %plotReservoirDetails(population,best_indv(test,1,1),1,best_indv(test,1,1),config);
    
    %% start GA
    for gen = 2:config.total_gens
        
        % redefine seed - some functions/scripts may reset the seed
        rng(gen,'twister');
        
        % reshape stored error to compare errors
        %cmp_error = reshape(store_error(test,gen-1,:),1,size(store_error,3));
        
        % random parent with mutation and crossover. create offspring
        %indv = [];
        
        indv1 = best_indv(test,gen-1,1);
        %indv2 = best_indv(test,gen-1,2);
        if config.parallel
            parfor o = 1:config.num_offspring
                %offspring(o) = config.recFcn(population(indv1),population(indv2),config);
                offspring(o) = config.mutFcn(population(indv1),config);
                % Evaluate and update fitness
                offspring(o).core_indx = o;
                offspring(o) = config.testFcn(offspring(o),config);
            end
        else
            for o = 1:config.num_offspring
                %offspring(o) = config.recFcn(population(indv1),population(indv2),config);
                offspring(o) = config.mutFcn(population(indv1),config);
                % Evaluate and update fitness
                offspring(o) = config.testFcn(offspring(o),config);
            end
        end
        population = [population(best_indv(test,gen-1,1:config.num_parents)) offspring];
        
        %update errors
        store_error(test,gen,:) =  store_error(test,gen-1,:);
        store_error(test,gen,:) = getError(config.error_to_check,population);%population(loser).val_error;
        [best(test,gen,:),best_indv(test,gen,:)] = sort(reshape(store_error(test,gen,:),1,config.pop_size));
        
        % print info
        if (mod(gen,config.gen_print) == 0)
            fprintf('Gen %d, time taken: %.4f sec(s)\n  Best indv: %.4f, Old bestt: %.4f, Best Error: %.4f \n',gen,toc/config.gen_print,getError('test',population(best_indv(test,gen,1))),getError('test',population(best_indv(test,gen-1,1))),best(test,gen,1));
            tic;
            if best(test,gen,1) ~= last_best
                % plot reservoir structure, task simulations etc.
                %plotReservoirDetails(population,best_indv(test,:,1),gen,best_indv(test,gen-1,1),config);
                drawnow
            end
            
            last_best = best(test,gen,1);
        end
        
        
        %save data
        if mod(gen,config.save_gen) == 0
            saveData(population,store_error,best,best_indv,config)
        end
        
        % if found best, stop
        if best(test,gen) == 0
            fprintf('Found solution. Gen %d, \n  Best indv: %d, Error: %.4f \n',gen,best_indv(test,gen,1),best(test,gen,1));
            plotReservoirDetails(population,best_indv(test,:,1),gen,best_indv(test,gen,1),config);
            return;
        end
    end
    
    % apply metrics to final population
    if config.record_metrics
        parfor pop_indx = 1:config.pop_size
            metrics(pop_indx,:) = getMetrics(population(pop_indx),config);
        end
    end
    
    if config.prune
        rng(1,'twister')
        pop_error = getError(config.error_to_check,population);
        
        [~,indx] = sort(pop_error);
        
        parfor p = 1:4
            [Pruned_best_individual{p}] = maxPruning(@testReservoir,{'train_error','val_error','test_error'},population(indx(p))...
                ,[population(indx(p)).train_error population(indx(p)).val_error population(indx(p)).test_error]...
                ,config.tolerance,config); %,old_W_fitness(p),old_Win_fitness(p),full_W{p}
        end
        % plotReservoirDetails(population,best_indv,gen,loser,config)
    end
end

    function saveData(population,store_error,best,best_indv,config)
        config.figure_array =[];
        if iscell(config.res_type)
            res_type ='Heterotic';
        else
            res_type = config.res_type;
        end
        
        save(strcat(res_type,'_',config.dataset,'_run',num2str(config.num_tests),'_gens',num2str(config.total_gens),'_',num2str(sum(config.num_reservoirs)),'Nres.mat'),...
            'population','store_error','best','best_indv','config','-v7.3');
    end

end
