%% Evolve best filter for the BSA algorithm
function [best,encodedData,decodedData] = filterGA(data,config)

%close all
rng(1,'twister')

figure1 = figure;

pop_size = 50;
num_gens = config.num_gens;
recRate = 0.5;
mutRate = 0.1;
deme_percent = 0.1;
deme = round(pop_size*deme_percent);



%% create pop of filters
for i = 1:pop_size
    population(i).order = randi([1 config.max_order]);
    population(i).passband = rand;
    population(i).threshold = rand;
    population(i).time_period = randi([1 config.max_period]);
    population(i).input_length = randi([1 population(i).time_period]);
    
    % evaluate
    [population(i).mse,~,~,population(i).f] = evaluateFilter(population(i),data);
    fprintf('Pop %d, Error: %.4f \n',i,population(i).mse);
     
    set(0,'currentFigure',figure1)
    plot(population(i).f); 
    drawnow
end

% store error that will be used as fitness in the GA
store_error(1,:) = [population.mse];
[best(1),best_indv(1)] = min(store_error(1,:));

%% run GA
for gen = 2:num_gens
    
    % reshape stored error to compare errors
    cmp_error = store_error(gen-1,:);
    
    % Tournment selection - pick two individuals
    equal = 1;
    while(equal) % find pair who are within deme range
        indv1 = randi([1 pop_size]);
        indv2 = indv1 + randi([1 deme]);
        if indv2 > pop_size
            indv2 = indv2 - pop_size; %loop around population ring if too big
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
    population(loser) = recFcn(population(winner),population(loser),recRate);
    population(loser) = mutFcn(population(loser),mutRate,config);
    
    %% Evaluate and update fitness
    [population(loser).mse,~,~,population(loser).f] = evaluateFilter(population(loser),data);
           
    % store error
    store_error(gen,:) =  store_error(gen-1,:);
    store_error(gen,loser) = population(loser).mse;
    [best(gen),best_indv(gen)] = min(store_error(gen,:));
    
    if mod(gen,100) == 0
        fprintf('Gen %d, Winner: %.4f, Loser: %.4f, Best Error: %.4f \n',gen,population(winner).mse,population(loser).mse,best(gen));
        set(0,'currentFigure',figure1)
        plot(population(best_indv(gen)).f);
        drawnow
    end
end


[mse,encodedData,decodedData,population(best_indv(gen)).f,mul_data] = evaluateFilter(population(best_indv(gen)),data);

best = population(best_indv(gen));

figure
subplot(1,3,1)
plot(data);
title('Scaled data');

subplot(1,3,2)
%stem(encodedData)
plot(encodedData,'r')
title('Spike trains');

subplot(1,3,3)
hold on
plot(mul_data,'b');
plot(decodedData,'r')
title('Reconstructed data');
hold off 

end


function l = recFcn(w,l,recRate)

W= w.order(:);
L = l.order(:);
pos = randperm(length(L),sum(rand(length(L),1) < recRate));         
L(pos) = W(pos);
l.order = reshape(L,size(l.order));

W= w.passband(:);
L = l.passband(:);
pos = randperm(length(L),sum(rand(length(L),1) < recRate));         
L(pos) = W(pos);
l.passband = reshape(L,size(l.passband));

W= w.threshold(:);
L = l.threshold(:);
pos = randperm(length(L),sum(rand(length(L),1) < recRate));         
L(pos) = W(pos);
l.threshold = reshape(L,size(l.threshold));

W= w.input_length;
L = l.input_length;
pos = randperm(length(L),sum(rand(length(L),1) < recRate));         
L(pos) = W(pos);
l.input_length = reshape(L,size(l.input_length));

W= w.time_period;
L = l.time_period;
pos = randperm(length(L),sum(rand(length(L),1) < recRate));         
L(pos) = W(pos);
l.time_period = reshape(L,size(l.time_period));


% check input period
if l.input_length > l.time_period
    l.input_length = l.time_period;
end
    
end

function filter = mutFcn(filter,mutRate,config)

order = filter.order(:);
pos = randperm(length(order),sum(rand(length(order),1) < mutRate));
order(pos) = round(mutate(order(pos),[1, config.max_order],config));
filter.order = reshape(order,size(filter.order));

passband = filter.passband(:);
pos = randperm(length(passband),sum(rand(length(passband),1) < mutRate));
passband(pos) = mutate(passband(pos),[0, 1],config);
filter.passband = reshape(passband,size(filter.passband));

threshold = filter.threshold(:);
pos = randperm(length(threshold),sum(rand(length(threshold),1) < mutRate));
threshold(pos) = mutate(threshold(pos),[0, 1],config);
filter.threshold = reshape(threshold,size(filter.threshold));

time_period = filter.time_period;
pos =  randperm(length(time_period),sum(rand(length(time_period),1) < mutRate));
time_period(pos) = round(mutate(time_period(pos),[1 config.max_period],config));
filter.time_period = reshape(time_period,size(filter.time_period));

input_length = filter.input_length;
pos =  randperm(length(input_length),sum(rand(length(input_length),1) < mutRate));
input_length(pos) = round(mutate(input_length(pos),[filter.time_period filter.time_period],config));
filter.input_length = reshape(input_length,size(filter.input_length));

% check input period
%if filter.input_length > filter.time_period
    filter.input_length = filter.time_period;
%end

end

function value = mutate(value,range,config)

if range(1)~=range(2)
    switch(config.mutate_type)
        case 'gaussian'
            for i = 1:length(value)
                flag = 1;
                while(flag)
                    t_value = value(i) + (range(1) + (range(2)-range(1))*randn);
                    
                    % check within range
                    if (t_value <= range(2)) && (t_value >= range(1))
                        flag = 0;
                    end
                end
                value(i) = t_value;
            end
        case 'uniform'
            value = range(1) + (range(2)-range(1))*rand;
    end
end
end
% figure
% plot(normData);
% title('Scaled data');
%
% figure
% %stem(encodedData)
% stem(encodedData,'r')
% title('Spike trains');
%
% figure
% hold on
% plot(normData,'b');
% plot(decodedData,'r')
% title('Reconstructed data');
% hold off
