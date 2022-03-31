% vowels best: T(560) : load('/mnt/lustre/users/md596/Magnets/multilayer/Evolve to Task/Single objective/STO_japanese_vowels_run17_gens1000_64nodes_het0_graph_Ring.mat')
% Diabetes best: T(98): load('/mnt/lustre/users/md596/Magnets/multilayer/Evolve to Task/Single objective/STO_diabetes_run18_gens1000_64nodes_het0_graph_basicLattice.mat')
% Narma-10 best: T(228): load('/mnt/lustre/users/md596/Magnets/multilayer/Evolve to Task/Single objective/STO_narma_10_run16_gens1000_64nodes_het1_graph_basicLattice.mat')

% plot confusion matrix
best_individual = population(6);

%states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
[~,states,test_sequence] = evalReservoir(best_individual,config);
ypred = states*best_individual.output_weights;

%% plot classification accuracy
% apply softmax
a = [];
for i = 1:size(ypred,1)
    a(i,:) = exp(ypred(i,:))/sum(exp(ypred(i,:)));
end
[~,ypred] = max(a,[],2);

% for conf mat
ypred_onehot = zeros(size(a));
for i = 1:size(ypred,1)
ypred_onehot(i,ypred(i)) = 1;
end

[~,y] = max(config.test_output_sequence,[],2);

% conf chart
figure
cm= confusionchart(y,ypred)
%cm.Title = 'PIMA Indians diabetes Classification';
cm.Title = 'Japanese Vowels: Accuracy 97.6%';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';

% conf mat
figure
plotconfusion(config.test_output_sequence',ypred_onehot')
 set(gca,'FontSize',16,'FontName','Arial')
xlabel('Target Speaker')
ylabel('Predicted Speaker')
title('')

%% plot namra performance
figure
plot(config.test_output_sequence(config.wash_out+1:end),'k')
hold on
plot(ypred,'r')
hold off
legend({'Target','Predicted'})
xlim([100 500])
xlabel('Time')
ylabel('Amplitude')
set(gca,'FontSize',16,'FontName','Arial')
title('NRMSE = 0.3644')

%% plot transient dynamics
load('/mnt/lustre/users/md596/Magnets/multilayer/Evolve to Task/Single objective/STO_narma_10_run16_gens1000_64nodes_het1_graph_basicLattice.mat')
best_individual = population(6);
final_states = collectSTOStates(best_individual,[zeros(100,1); 1; zeros(200,1)],config,[zeros(100,1); 1; zeros(200,1)]);
subplot(1,2,1)
plot(final_states(:,1:end-2))

load('/mnt/lustre/users/md596/Magnets/multilayer/Evolve to Task/Single objective/STO_diabetes_run18_gens1000_64nodes_het0_graph_basicLattice.mat')
best_individual = population(2);
final_states = collectSTOStates(best_individual,[zeros(100,8); ones(1,8); zeros(200,8)],config,[zeros(100,1); 1; zeros(200,1)]);
subplot(1,2,2)
plot(final_states(:,51:end-2))