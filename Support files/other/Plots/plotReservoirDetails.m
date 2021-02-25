function F = plotReservoirDetails(population,best_indv,gen,loser,config)

% individual to print - maybe cell if using MAPelites
if iscell(population(best_indv(gen)))
    best_individual = population{best_indv(gen)};
    loser_individual = population{loser};
else
    best_individual = population(best_indv(gen));
    loser_individual = population(loser);
end

desktop     = com.mathworks.mde.desk.MLDesktop.getInstance;
cw          = desktop.getClient('Command Window');
xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
h_cw        = handle(xCmdWndView,'CallbackProperties');
set(h_cw, 'KeyPressedCallback', @CmdKeyCallback);

CmdKeyCallback('reset');
fprintf('Press any key to skip simulation \n')

set(0,'currentFigure',config.figure_array(1))

% plot task specific details
switch(config.dataset)
    
    case {'spiral','spiral_multi_class'}
        set(gcf,'Position',[138 543 1359 435])
        [x,y] = meshgrid(floor(min(min(config.train_input_sequence))):0.1:max(max(config.train_input_sequence)),...
            floor(min(min(config.train_input_sequence))):0.1:max(max(config.train_input_sequence)));
        [X2,Y2] = meshgrid(floor(min(min(config.train_input_sequence))):0.01:max(max(config.train_input_sequence)),...
            floor(min(min(config.train_input_sequence))):.01:max(max(config.train_input_sequence)));
        
        subplot(1,3,1)
        % plot boundaries
        boundary_states = config.assessFcn(best_individual,[x(:),y(:)],config,config.train_output_sequence);
        boundary_sequence = boundary_states*best_individual.output_weights;
        [a,boundary_class] =max(boundary_sequence,[],2);
        
        if best_individual.n_output_units == 2
            boundary_class(boundary_class~=2) = 0;
            outData = interp2(reshape(x,size(x,1),size(x,2)), reshape(y,size(y,1),size(y,2)), reshape(boundary_class,size(x,1),size(x,2)), X2, Y2, 'linear')-1;
            imagesc(outData*0.5);
            set(gca,'YDir','normal')
        else
            outData = interp2(reshape(x,size(x,1),size(x,2)), reshape(y,size(y,1),size(y,2)), reshape(boundary_class,size(x,1),size(x,2)), X2, Y2, 'linear');
            imagesc(outData);
            set(gca,'YDir','normal')
        end
        
        % plot classes
        train_states = config.assessFcn(best_individual,config.train_input_sequence,config,config.train_output_sequence);
        train_sequence = train_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(train_sequence,config.train_output_sequence,config);
        
        j = 1;
        b_max = max(max(config.train_input_sequence));
        b_min = floor(min(min(config.train_input_sequence)));
        if b_min < 0
            n = b_max*100;
            mag = 100;
        else
            n = 0;
            mag = b_max*100;
        end
        
        hold on
        s = scatter(mag.*rot90(config.train_input_sequence(config.wash_out+1:end,1),j)+n,mag.*rot90(config.train_input_sequence(config.wash_out+1:end,2),j)+n,30,rot90(system_output,j),'filled');
        xticks(1:(n*2+n*0.25)/9:n*2+ n*0.25)
        xticklabels(round(b_min:b_max/4:b_max))
        yticks(1:(n*2+n*0.25)/9:n*2+ n*0.25)
        yticklabels(round(b_min:b_max/4:b_max))
        
        s.LineWidth = 1;
        s.MarkerEdgeColor = 'k';
        hold off
        
        subplot(1,3,2)
        % plot boundaries
        if best_individual.n_output_units == 2
            imagesc(outData*0.5);
            set(gca,'YDir','normal')
        else
            imagesc(outData);
            set(gca,'YDir','normal')
        end
        
        % plot classes
        val_states = config.assessFcn(best_individual,config.val_input_sequence,config,config.val_output_sequence);
        val_sequence = val_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(val_sequence,config.val_output_sequence,config);
        
        j = 1;
        
        hold on
        s = scatter(mag.*rot90(config.val_input_sequence(config.wash_out+1:end,1),j)+n,mag.*rot90(config.val_input_sequence(config.wash_out+1:end,2),j)+n,30,rot90(system_output,j),'filled');
        xticks(1:(n*2+n*0.25)/9:n*2+ n*0.25)
        xticklabels(round(b_min:b_max/4:b_max))
        yticks(1:(n*2+n*0.25)/9:n*2+ n*0.25)
        yticklabels(round(b_min:b_max/4:b_max))
        
        s.LineWidth = 1;
        s.MarkerEdgeColor = 'k';
        hold off
        
        subplot(1,3,3)
        % plot boundaries
        if best_individual.n_output_units == 2
            imagesc(outData*0.5);
            set(gca,'YDir','normal')
        else
            imagesc(outData);
            set(gca,'YDir','normal')
        end
        
        % plot classes
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(test_sequence,config.test_output_sequence,config);
        
        j = 1;
        
        hold on
        s = scatter(mag.*rot90(config.test_input_sequence(config.wash_out+1:end,1),j)+n,mag.*rot90(config.test_input_sequence(config.wash_out+1:end,2),j)+n,30,rot90(system_output,j),'filled');
        xticks(1:(n*2+n*0.25)/9:n*2+ n*0.25)
        xticklabels(round(b_min:b_max/4:b_max))
        yticks(1:(n*2+n*0.25)/9:n*2+ n*0.25)
        yticklabels(round(b_min:b_max/4:b_max))
        
        s.LineWidth = 1;
        s.MarkerEdgeColor = 'k';
        hold off
        drawnow
        
        
    case 'multi_signal'
        test_states = config.assessFcn(best_individual,config.train_input_sequence,config,config.train_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        
        for i = 1:4
            subplot(3,2,i)
            plot(config.train_output_sequence(config.wash_out+1:end,i),'--')
            hold on
            plot(test_sequence(:,i))
            hold off
            legend({'Target','Predicted'})
        end
        
        subplot(3,2,[5 6])
        xdft = fft(test_states);
        xdft = xdft(1:length(test_states)/2+1);
        freq = 0:config.Fs/length(test_states):config.Fs/2;
        plot(freq,abs(xdft));
        xlabel('Hz');
        
%         xdft = fft(test_states);
%         plot(abs(test_states))

%     case {'laser','narma_10','narma_30', 'sunspot'}
%         
%         test_states = config.assessFcn(best_individual,config.train_input_sequence,config,config.train_output_sequence);
%         test_sequence = test_states*best_individual.output_weights;
%         plotresponse(num2cell(config.train_output_sequence(config.wash_out+1:end,:)'),num2cell(test_sequence'))
%         
    case 'lorenz'
        test_states = config.assessFcn(best_individual,config.train_input_sequence,config,config.train_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        
        subplot(1,2,1)
        plot(config.train_output_sequence(config.wash_out+1:end,:))
        hold on
        plot(test_sequence)
        hold off
        
        subplot(1,2,2)
        plot3(test_sequence(:,1),test_sequence(:,2),test_sequence(:,3))
        
    case 'autoencoder'
        
        test_states = config.assessFcn(best_individual,config.train_input_sequence,config,config.train_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        
        rand_input = randn(size(config.train_input_sequence));
        rand_states = config.assessFcn(best_individual,rand_input,config,config.train_output_sequence);
        rand_sequence = rand_states*best_individual.output_weights;
        if size(config.train_input_sequence,2) == 3072
            for i = 1:6
                subplot(4,6,i)
                imshow(reshape(config.train_input_sequence(i*10,:),sqrt(size(config.train_output_sequence,2)/3),sqrt(size(config.train_output_sequence,2)/3),3));
                title('Actual')
                
                subplot(4,6,6+i)
                imshow(reshape(test_sequence(i*10,:),sqrt(size(config.train_output_sequence,2)/3),sqrt(size(config.train_output_sequence,2)/3),3));
                title('Recovered')
                
                subplot(4,6,12+i)
                imshow(reshape(rand_input(i*10,:),sqrt(size(config.train_output_sequence,2)/3),sqrt(size(config.train_output_sequence,2)/3),3));
                title('Random Input')
                
                subplot(4,6,18+i)
                imshow(reshape(rand_sequence(i*10,:),sqrt(size(config.train_output_sequence,2)/3),sqrt(size(config.train_output_sequence,2)/3),3));
                title('Output: Random Input')
            end
        else
            for i = 1:6
                subplot(4,6,i)
                imshow(reshape(config.train_input_sequence(i*10,:),sqrt(size(config.train_output_sequence,2)),sqrt(size(config.train_output_sequence,2))));
                title('Actual')
                
                subplot(4,6,6+i)
                imshow(reshape(test_sequence(i*10,:),sqrt(size(config.train_output_sequence,2)),sqrt(size(config.train_output_sequence,2))));
                title('Recovered')
                
                subplot(4,6,12+i)
                imshow(reshape(rand_input(i*10,:),sqrt(size(config.train_output_sequence,2)),sqrt(size(config.train_output_sequence,2))));
                title('Random Input')
                
                subplot(4,6,18+i)
                imshow(reshape(rand_sequence(i*10,:),sqrt(size(config.train_output_sequence,2)),sqrt(size(config.train_output_sequence,2))));
                title('Output: Random Input')
            end
        end
        
    case 'image_painting'
        
        test_states = config.assessFcn(best_individual,config.train_input_sequence,config,config.train_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        
        subplot(1,2,1)
        imshow(reshape(config.train_output_sequence,sqrt(size(config.train_output_sequence,1)),sqrt(size(config.train_output_sequence,1)),3));
        title('Actual')
        
        subplot(1,2,2)
        imshow(reshape(test_sequence,sqrt(size(test_sequence,1)),sqrt(size(test_sequence,1)),3));
        title('Predicted')
        
    case 'pole_balance'
        
        config.run_sim = 1;
        config.testFcn(best_individual,config);
        config.run_sim = 0;
        
    case 'attractor'
        
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        
        subplot(1,3,1)
        plot(config.test_output_sequence(config.wash_out+1:end,:),'r')
        hold on
        plot(test_sequence,'b')
        plot(config.test_input_sequence(config.wash_out+1:end,:),'g')
        hold off
        legend({'Target','Output','Input'})
        title('I/O data');
        drawnow
        
        subplot(1,3,2)
        X = config.test_output_sequence(config.wash_out+1:end,:);
        T = test_sequence;
        switch(size(X,2))
            case 3
                plot3(X(:,1),X(:,2),X(:,3),'r');
                hold on
                plot3(T(:,1),T(:,2),T(:,3),'b');
                hold off
                xlabel('X'); ylabel('Y'); zlabel('Z');
            case 2
                plot(X(:,1),X(:,2),'r');
                hold on
                plot(T(:,1),T(:,2),'b');
                hold off
                xlabel('X'); ylabel('Y');
            case 1
                plot(X(:,1),'r');
                hold on
                plot(T(:,1),'b');
                hold off
                xlabel('X');
        end
        
        %axis equal;
        legend('Target','output')
        grid;
        title('Attractor');
        
        subplot(1,3,3)
        plot(test_states)
        title('Reservoir states');
        
        drawnow
        
    case 'robot'
        
        config.run_sim = 1;
        config.testFcn(best_individual,config);
        config.run_sim = 0;
        
    case 'CPPN'
        
        subplot(1,3,1)
        G1 = digraph(best_individual.W{1});
        [X_grid,Y_grid] = ndgrid(linspace(-1,1,sqrt(size(G1.Nodes,1))));
        
        p = plot(G1,'XData',X_grid(:),'YData',Y_grid(:));
        p.EdgeCData = G1.Edges.Weight;
        colormap(gca,bluewhitered);
        colorbar
        title('Best')
        
        subplot(1,3,2)
        G2 = digraph(loser_individual.W{1});
        [X_grid,Y_grid] = ndgrid(linspace(-1,1,sqrt(size(G2.Nodes,1))));
        
        p = plot(G2,'XData',X_grid(:),'YData',Y_grid(:));
        p.EdgeCData = G2.Edges.Weight;
        colormap(gca,bluewhitered);
        colorbar
        title('loser')
        
        subplot(1,3,3)
        imagesc(best_individual.W{1})
        colormap(gca,bluewhitered);
        set(gca,'YDir','normal')
        
        drawnow
        %return;
        
    case 'image_gaussian'
        states = config.assessFcn(best_individual,config.test_input_sequence,config);
        output = states*best_individual.output_weights;
        
        for image_indx = 1:3
            subplot(3,2,(image_indx*2)-1)
            imagesc(reshape(output(image_indx,:),sqrt(size(output,2)),sqrt(size(output,2))));
            xlabel('Reservoir Output')
            subplot(3,2,image_indx*2)
            imagesc(reshape(config.test_output_sequence(image_indx,:),sqrt(size(output,2)),sqrt(size(output,2))));
            xlabel('Target Output')
        end
        
    case 'franke_fcn'
        
        states = config.assessFcn(best_individual,config.test_input_sequence,config);
        output = states*best_individual.output_weights;
        
        subplot(1,2,1)
        x = config.test_input_sequence(config.wash_out+1:end,1);
        y = config.test_input_sequence(config.wash_out+1:end,2);
        z = config.test_output_sequence(config.wash_out+1:end);
        
        xv = linspace(min(x), max(x), 20);
        yv = linspace(min(y), max(y), 20);
        [X,Y] = meshgrid(xv, yv);
        Z = griddata(x,y,z,X,Y);
        surf(X, Y, Z);
        
        subplot(1,2,2)
        z = output;
        xv = linspace(min(x), max(x), 20);
        yv = linspace(min(y), max(y), 20);
        [X,Y] = meshgrid(xv, yv);
        Z = griddata(x,y,z,X,Y);
        surf(X, Y, Z);
        
        drawnow
        
    case {'boston_housing'} % for regression problems
        
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        plotregression(config.test_output_sequence,test_sequence)
        drawnow
        
        
    case {'MNIST','hand_digits'}
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(test_sequence,config.test_output_sequence,config);
        
        title(strcat('Error: ', num2str(err)))
        n = size(config.test_input_sequence,1)/10;
        for i = 1:n
            subplot(ceil(sqrt(n)),ceil(sqrt(n)),i)
            imshow(reshape(config.test_input_sequence(10*i,:),sqrt(size(config.test_input_sequence,2)),sqrt(size(config.test_input_sequence,2))))
            if system_output(10*i) == desired_output(10*i)
                xlabel(strcat('Predicted: ', num2str(system_output(10*i))-1),'Color','g')
            else
                xlabel(strcat('Predicted: ', num2str(system_output(10*i))-1),'Color','r')
            end
        end
        
        drawnow
        
    case 'cifar10'
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(test_sequence,config.test_output_sequence,config);
        
        title(strcat('Error: ', num2str(err)))
        n = size(config.test_input_sequence,1)/10;
        
        labels = load('batches.meta.mat');
        
        for i = 1:n
            subplot(ceil(sqrt(n)),ceil(sqrt(n)),i)
            imshow(reshape(config.test_input_sequence(10*i,:),sqrt(size(config.test_input_sequence,2)/3),sqrt(size(config.test_input_sequence,2)/3),3))
            if system_output(10*i) == desired_output(10*i)
                xlabel(strcat('Predicted: ', labels.label_names{system_output(10*i)}),'Color','g')
            else
                xlabel(strcat('Predicted: ', labels.label_names{system_output(10*i)}),'Color','r')
            end
        end
        
        drawnow
        
    case 'signal_classification'
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(test_sequence,config.test_output_sequence,config);
        
        
        plot(desired_output,'b')
        hold on
        plot(system_output,'r')
        hold off
        legend('Target','Output')
        
    case {'MSO1','MSO2','MSO3','MSO4','MSO5','MSO6','MSO7','MSO8','MSO9','MSO10','MSO11','MSO12'}
        
        test_states = config.assessFcn(best_individual,config.test_input_sequence,config,config.test_output_sequence);
        test_sequence = test_states*best_individual.output_weights;
        [err,system_output,desired_output] = calculateError(test_sequence,config.test_output_sequence,config);
        
        %plot(config.test_input_sequence,[0.9290 0.6940 0.1250],'--')
        %hold on
        subplot(2,2,[1 2])
        plot(desired_output,'--')
        hold on
        plot(system_output,'-')
        plot(config.test_input_sequence(config.wash_out+1:end,:),'-')
        hold off
        legend('Target','Output','Input')
        
        subplot(2,2,3)
        Fs = 2000; % sample frequency
        %test = test_states(:,sum(test_states,1)>0);
        test = test_states;
        L = size(test,1);
        xdft = fft(test);
        P2 = abs(xdft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        plot(f,P1)
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')

        subplot(2,2,4)
        %phase=angle(P1);
        phase = fftshift(xdft);
        ly = length(phase);
        f1 = (-ly/2:ly/2-1)/ly*Fs;
        tol = 1e-6;
        phase(abs(phase) < tol) = 0;
        theta = angle(phase);
        plot(f1,theta/pi)
        %plot(f1,abs(phase))
        grid on
        xlabel('Frequency(Hz)');
        ylabel('Phase (rad)')
%         xdft = fft(test_states);
%         xdft = xdft(1:length(test_states)/2+1);
%         freq = 0:Fs/length(test_states):10000/2;
%         plot(freq,abs(xdft));
%         xlabel('Hz');
%         p = plot([desired_output system_output]);
%         p(1).Color = [0.9290 0.6940 0.1250];
%         p(2).Color =[0 0.4470 0.7410];
%         p(3).Color =[0.3010 0.7450 0.9330];
end

%draw task details
drawnow

% plot reservoir details
if ~iscell(config.res_type)
    set(0,'currentFigure',config.figure_array(2))
    switch(config.res_type)
        case 'Graph'
            plotGridNeuron(config.figure_array(2),population,best_indv(gen),loser,config)
            
        case {'2D_CA','GOL'}
            colormap('bone')
            
            if config.run_sim
                
                set(gcf,'position',[-1005 458 981 335])
                [states,~,extra_states] = config.assessFcn(population(best_indv(gen)),config.test_input_sequence,config);
                for i = 1:size(states,1)
                    subplot(1,2,1)
                    t_state = states(i,1:end-population(best_indv(gen)).n_input_units);
                    imagesc(reshape(t_state,sqrt(size(t_state,2)),sqrt(size(t_state,2))));
                    %title(strcat('n = ',num2str(i)))
                    subplot(1,2,2)
                    t_state = extra_states(i,1:end);
                    imagesc(reshape(t_state,sqrt(size(t_state,2)),sqrt(size(t_state,2))));
                    
                    drawnow
                    if config.film
                        F(i) = getframe(gcf);
                    else
                        F =[];
                    end
                    if CmdKeyCallback()
                        i = size(states,1);
                    end
                    pause(0.05);
                end
            end
            
        case {'Ising'}
            
            if config.run_sim
                colormap('bone')
                set(gcf,'position',[-737 236 713 557])
                [states] = config.assessFcn(population(best_indv(gen)),config.test_input_sequence,config);
                for i = 1:size(states,1)
                    t_state = states(i,1:end-population(best_indv(gen)).n_input_units);
                    imagesc(reshape(t_state,sqrt(size(t_state,2)),sqrt(size(t_state,2))));
                    colorbar
                    caxis([-1 1])
                    drawnow
                    if config.film
                        F(i) = getframe(gcf);
                    else
                        F =[];
                    end
                    if CmdKeyCallback()
                        i = size(states,1);
                    end
                    pause(0.05);
                end
            end
            
        case 'BZ'
            
            config.plot_states = 1;
            
            states = config.assessFcn(population(best_indv(gen)),config.test_input_sequence,config);
            
            if config.run_sim
                set(0,'currentFigure',config.figure_array(1))
                states = config.assessFcn(population(best_indv(gen)),config.test_input_sequence,config);
                
                for i = 1:size(states,1)
                    p = reshape(states(i,1:end-population(best_indv(gen)).n_input_units),sqrt(population(best_indv(gen)).nodes),sqrt(population(best_indv(gen)).nodes),3);
                    image(uint8(255*hsv2rgb(p)));
                    drawnow;
                    if CmdKeyCallback()
                        i = size(states,1);
                    end
                    if config.film
                        F(i) = getframe;
                    else
                        F =[];
                    end
                end
            end
            
        case {'RoR','Pipeline','Ensemble'}
            
            plotRoR(config.figure_array(2),best_individual,loser_individual,config);
            
            if size(config.num_nodes,2) > 1
                
                G = digraph((best_individual.W_switch.*best_individual.W_scaling));
                
                set(0,'currentFigure',config.figure_array(1))
                hold off
                title('Reservoir structure')
               
                plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',abs(G.Edges.Weight))
                drawnow
            end
            
        case 'RoRmin'
            
            ax1 = subplot(1,3,1);
             imagesc((best_individual.input_weights.*best_individual.input_scaling)');
            colormap(ax1,bluewhitered)
            colorbar
            xlabel('Input mapping')
            ax2 = subplot(1,3,2);
            imagesc(best_individual.W.*best_individual.W_scaling);
            colormap(ax2,bluewhitered)
            colorbar
            xlabel('Internal weights')
            ax3 = subplot(1,3,3);
            imagesc(best_individual.output_weights);
            colormap(ax3,bluewhitered)
            colorbar
            xlabel('Output mapping')
            
        case {'RBN','elementary_CA'}
            plotRBN(best_individual,config)
            
        case 'Wave'
            
            config.wave_sim_speed = 1;
            
            %% plot input locations
            subplot(1,2,1)
            indx=1;
            % write input values for each location
            input = best_individual.input_scaling(indx)*(best_individual.input_weights{indx}*[config.test_input_sequence repmat(best_individual.bias_node,size(config.test_input_sequence,1),1)]')';
            
            % change input widths
            node_grid_size = sqrt(best_individual.nodes(indx));
            for n = 1:size(input,1)
                b_max = reshape(input(n,:),node_grid_size,node_grid_size);
                f_pos = find(b_max);
                input_matrix_2d = b_max;
                for p = 1:length(f_pos)
                    t = zeros(size(b_max));
                    t(f_pos(p)) = b_max(f_pos(p));
                    [t] = adjustInputShape(t,best_individual.input_widths{indx}(f_pos(p)));
                    input_matrix_2d = input_matrix_2d + t;
                end
                input(n,:) = input_matrix_2d(:);
            end
            
            imagesc(reshape(max(abs(input)),node_grid_size,node_grid_size))
            title('Input location')
            colormap(bluewhitered)
            colorbar
            
            %% plot states
            subplot(1,2,2)
            % run wave reservoir on task
            switch(config.dataset)
                case 'robot'
                    [~,states] = robot(best_individual,config);
                case 'pole_balance'
                    [~,states]= poleBalance(best_individual,config);
                otherwise
                    states = config.assessFcn(best_individual,config.test_input_sequence,config);
            end
            
            if config.num_reservoirs
                states = states(:,1:config.num_nodes(1)+best_individual.n_input_units);
            end
            
            %plot
            set(0,'currentFigure',config.figure_array(2))
            if config.add_input_states
                h=surf(reshape(states(1,1:end-best_individual.n_input_units),node_grid_size,node_grid_size));
            else
                h=surf(reshape(states(1,1:end),node_grid_size,node_grid_size));
            end
            
            
            change_scale = 0;
            
            i = 2;
            colormap(gca,'bone'); %
            %set(gca,'visible','off')
            set(gca,'XColor', 'none','YColor','none')
            shading interp
            %lighting phong;
            %material shiny;
            %lightangle(-45,30)
            while(i < size(states,1))
                if mod(i,config.wave_sim_speed) == 0
                    if config.add_input_states
                        newH = reshape(states(i,1:end-best_individual.n_input_units),node_grid_size,node_grid_size);
                        set(h,'zdata',newH,'facealpha',0.65);
                        
                        if change_scale
                            set(gca, 'xDir', 'reverse',...
                                'camerapositionmode','manual','cameraposition',[1 1 max((states(i,1:end-best_individual.n_input_units)))]);
                            axis([1 node_grid_size 1 node_grid_size min((states(i,1:end-best_individual.n_input_units))) max((states(i,1:end-best_individual.n_input_units)))]);
                        else
                            set(gca, 'xDir', 'reverse',...
                                'camerapositionmode','manual','cameraposition',[1 1 max(max(states(:,1:end-best_individual.n_input_units)))]);
                            axis([1 node_grid_size 1 node_grid_size min(min(states(:,1:end-best_individual.n_input_units))) max(max(states(:,1:end-best_individual.n_input_units)))]);
                        end
                    else
                        newH = reshape(states(i,1:end),node_grid_size,node_grid_size);
                        set(h,'zdata',newH,'facealpha',0.65);
                        if change_scale
                            set(gca, 'xDir', 'reverse',...
                                'camerapositionmode','manual','cameraposition',[1 1 max(states(i,1:end))]);
                            axis([1 node_grid_size 1 node_grid_size min(states(i,1:end)) max(states(i,1:end))]);
                        else
                            set(gca, 'xDir', 'reverse',...
                                'camerapositionmode','manual','cameraposition',[1 1 max(max(states(:,1:end)))]);
                            axis([1 node_grid_size 1 node_grid_size min(min(states(:,1:end))) max(max(states(:,1:end)))]);
                        end
                    end
                    if config.film
                        F(i) = getframe;
                    else
                        F =[];
                    end
                    drawnow
                    pause(0.05)
                end
                
                if CmdKeyCallback()
                    i = size(states,1);
                end
                
                i = i +1;
            end
            
        case 'MM'
            
            config.parallel = 0;
            config.plot_states = 1;
            states = config.assessFcn(population(best_indv(gen)),config.test_input_sequence,config);
            states = states(:,1:end-best_individual.n_input_units);
            
        case 'Oregonator'
            config.plot_states = 1;
            states = config.assessFcn(population(best_indv(gen)),config.test_input_sequence,config);
    end
    
    
    if config.film
        v = VideoWriter('ReservoirPlot_video');%,'MPEG-4');
        v.Quality = 100;
        v.FrameRate = 15;
        open(v);
        writeVideo(v,F);
        close(v);
    end
    
end

% draw reservoir details
drawnow

end

function Value = CmdKeyCallback(ObjectH, EventData)

persistent KeyPressed

switch nargin
    case 0
        Value = ~isempty(KeyPressed);
    case 1
        KeyPressed = [];
    case 2
        KeyPressed = true;
    otherwise
        error('Programming error');
end
end

function plotFFT(states,config)
subplot(1,2,1)
imagesc(states(:,1:end-config.task_num_inputs)')
colorbar
subplot(1,2,2)
imagesc(real(fft(states(:,1:end-config.task_num_inputs)))')
colorbar
colormap(bluewhitered)
end