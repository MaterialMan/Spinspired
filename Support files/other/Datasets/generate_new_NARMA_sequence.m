function [inputSequence, outputSequence] = generate_new_NARMA_sequence(sequenceLength, memoryLength, lowerDist,upperDist)
%  Generates a sequence using a nonlinear autoregressive moving average
% (NARMA) model. The sequence at the beginning includes a  ramp-up
% transient, which should be deleted if necessary. The NARMA equation to be
% used must be hand-coded into this function (at bottom)
% inputs:
% sequenceLength: a natural number, indicating the length of the
% sequence to be generated
% memoryLength: a natural number indicating the dependency length
%
% outputs:
% InputSequence: array of size sequenceLength x 2. First column contains
%                uniform noise in [0,1] range, second column contains bias
%                input (all 1's)
% tmp_outputSequence: array of size sequenceLength x 1 with the NARMA output
%
% usage example:
% [a b] = generate_linear_sequence(1000,10) ;
%
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision H. Jaeger Feb 23, 2007

%set default distribution
if (nargin<3)
    lowerDist = 0;
    upperDist = 0.5;
end

%previously used to initiate random seed
%rng(1,'twister');
washout =0;

%%%% create input (set to between a distribution given by default or user
% use the input sequence to drive a NARMA equation
inputSequence = (upperDist-lowerDist)*rand(sequenceLength+washout,1)+lowerDist;

for n = 1:length(memoryLength)
    
    tmp_outputSequence = zeros(sequenceLength+washout,1);
    
    switch(memoryLength(n))
        %% Settings according to: A Comparative Study of Reservoir Computing for Temporal Signal Processing
        case 20
            for i = memoryLength(n)+1 : sequenceLength+washout
                
                middleSum(i) = sum(tmp_outputSequence(i-(memoryLength(n)-1):i));
                
                %10th order NARMA variables
                tmp_outputSequence(i) = tanh(0.3*tmp_outputSequence(i-1) + ...
                    0.05*tmp_outputSequence(i-1)*middleSum(i-1) + ...
                    1.5*inputSequence(i-(memoryLength(n)))*inputSequence(i-1)+...
                    0.01);%+ noise(i+1);
                %tmp_outputSequence(i+1) = tmp_outputSequence(i+1) ;
                
            end
            
        case {10,5} %[0,0.5]
            %% NRMSE = 0.173 to 0.14 in Atiya's first experiment (NARMA 10)
            %% NRMSE = 0.4 eqivil to shift regesiter. Thus, lower requires non-linearity (NARMA 10)...
            %NRMSE ~0.15. Results of simulated system. Information processing using a single dynamical node as complex system
            
            for i = memoryLength(n)+1 : sequenceLength+washout
                
                middleSum(i) = sum(tmp_outputSequence(i-(memoryLength(n)-1):i));
                
                %10th order NARMA variables
                tmp_outputSequence(i) = 0.3*tmp_outputSequence(i-1) + ...
                    0.05*tmp_outputSequence(i-1)*middleSum(i-1) + ...
                    1.5*inputSequence(i-(memoryLength(n)))*inputSequence(i-1)+...
                    0.1 ;%+ noise(i+1);
                %tmp_outputSequence(i+1) = tmp_outputSequence(i+1) ;
                
            end
            
        case {30,40}
            for i = memoryLength(n)+1 : sequenceLength+washout
                
                middleSum(i) = sum(tmp_outputSequence(i-(memoryLength(n)-1):i));
                
                %10th order NARMA variables
                tmp_outputSequence(i) = 0.2*tmp_outputSequence(i-1) + ...
                    0.004*tmp_outputSequence(i-1)*middleSum(i-1) + ...
                    1.5*inputSequence(i-(memoryLength(n)))*inputSequence(i-1)+...
                    0.001 ;%+ noise(i+1);
                %tmp_outputSequence(i+1) = tmp_outputSequence(i+1) ;
                
            end
            
        otherwise
            
            %         inputSequence = (upperDist-lowerDist)*rand(sequenceLength+memoryLength,1)+lowerDist;
            %         tmp_outputSequence = zeros(sequenceLength+memoryLength,1);
            
            for i = memoryLength(n)+1 : sequenceLength+washout
                
                %collect delay
                for n = 1:memoryLength(n)-1
                    tmp_outputSequence(i+1)= tmp_outputSequence(i+1)+ tmp_outputSequence(i-n);
                end
                
                tmp_outputSequence(i+1) = (tmp_outputSequence(i+1)*(0.05*tmp_outputSequence(i)));
                tmp_outputSequence(i+1) = tmp_outputSequence(i+1)+(0.3*tmp_outputSequence(i))...
                    + 1.5*inputSequence(i-(memoryLength(n)-1))*inputSequence(i) + 0.1;
                                
            end                       
    end  
    
    
    outputSequence(:,n) = tmp_outputSequence(max(memoryLength)+washout+1:end);
end

inputSequence = inputSequence(max(memoryLength)+washout+1:end);
