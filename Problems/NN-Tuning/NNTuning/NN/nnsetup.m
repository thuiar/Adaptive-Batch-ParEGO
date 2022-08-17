function nn = nnsetup(X)
%NNSETUP creates a Feedforward Backpropagate Neural Network
% nn = nnsetup(architecture) returns an neural network structure with n=numel(architecture)
% layers, architecture being a n x 1 vector of layer sizes e.g. [784 100 10]

% happywhy

% INPUT: X is a size of 10*1 vector whose elements are the number of
% hidden_layers, number of neurous per hidden layer, dropout,
% log_learning_rate, logL1_weight_regularization and
% logL2_weight_regularization, respectively

%  X range:
%   options.xmin = [1 50  0   -10 -10 -10]
%   options.xmax = [3 500 0.9 0  0  0]


%% NN hyper-paramters by happywhy

    nn.hidden_layers                    = round(X(1));
    nn.hidden_neuros                    = round(X(2)).*ones(1,nn.hidden_layers);
    nn.size                             = [784, nn.hidden_neuros, 10];
    nn.n                                = numel(nn.size); % input  784 and output 10
    nn.activation_function              = 'tanh_opt';   %  Activation functions of hidden layers: 'sigm' (sigmoid) or 'tanh_opt' (optimal tanh).
    nn.learningRate                     = exp(X(4));    % learning rate
%     nn.momentum                         = 0.5;          %  Momentum
    nn.scaling_learningRate             = 1;            %  Scaling factor for the learning rate (each epoch)
    nn.weightPenaltyL2                  = exp(X(5));            %  L2 regularization
    nn.nonSparsityPenalty               = 0;            %  Non sparsity penalty
    nn.sparsityTarget                   = 0.05;         %  Sparsity target
    nn.inputZeroMaskedFraction          = 0;            %  Used for Denoising AutoEncoders
    nn.dropoutFraction                  = X(3);            %  Dropout level (http://www.cs.toronto.edu/~hinton/absps/dropout.pdf)
    nn.testing                          = 0;            %  Internal variable. nntest sets this to one.
%     nn.output                           = 'softmax';       %  output unit 'sigm' (=logistic), 'softmax' and 'linear'
    
    nn.numepoches = round(X(6));                    %  Number of full sweeps through data
    nn.batchsize = round(X(7));                     %  Take a mean gradient step over this many samples
    nn.outputpara = round(X(8));                      % output chosen parameter
    nn.momentum  =X(9) ;          %  Momentum
    nn.active_func = round(X(10)); %  Activation functions of hidden layers: 'sigm' (sigmoid) or 'tanh_opt' (optimal tanh).
 
    %  output unit 'sigm' (=logistic), 'softmax' and 'linear'
    if nn.outputpara ==0
        nn.output                           = 'softmax';
    elseif nn.outputpara ==1
        nn.output                           = 'sigm'; 
    elseif nn.outputpara ==2
        nn.output                           = 'linear';
    end

    if nn.active_func == 0
         nn.activation_function              = 'tanh_opt';
    elseif nn.active_func == 1
         nn.activation_function              = 'sigm';
    end
    

    for i = 2 : nn.n   
        % weights and weight momentum
        nn.W{i - 1} = (rand(nn.size(i), nn.size(i - 1)+1) - 0.5) * 2 * 4 * sqrt(6 / (nn.size(i) + nn.size(i - 1)));
        nn.vW{i - 1} = zeros(size(nn.W{i - 1}));
        
        % average activations (for use with sparsity)
        nn.p{i}     = zeros(1, nn.size(i));   
    end
end
