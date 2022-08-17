function [er,predict_time] = EvaluateModelError(X,plot)
% This function is used to obtain the training error of NN

% INPUT: X is a size of 10*1 vector whose elements are the number of
% hidden_layers, number of neurous per hidden layer, dropout,
% log_learning_rate, logL1_weight_regularization and
% logL2_weight_regularization, num_epoche,batch size, 
% the output parameter, momnetum, activation_function,respectively
%  

%  X range:
%  xmin = [1  50    0     -10  -10  3  800    0   0.4   0]
%  xmax = [3  500  0.9  0     0      6  1200  2   0.6   1]

%----------------------------Reference -----------------------------------
% Palm R B. Prediction as a candidate for learning deep hierarchical models of 
% data. Technical University of Denmark, 2012, 5. The source code is from:
% https://github.com/rasmusbergpalm/DeepLearnToolbox
%-------------------------------------------------------------------------
[row,~] = size(X); % the size of the popdec

load mnist_uint8;

train_x = double(train_x) / 255;
test_x  = double(test_x)  / 255;
train_y = double(train_y);
test_y  = double(test_y);

% normalize
[train_x, data_mu, data_sigma] = zscore(train_x);
test_x = normalize(test_x, data_mu, data_sigma);

%% NN sigmoid activation and plotting of validation and training error
% split training data into training and validation data
vx   = train_x(1:10000,:);
tx = train_x(10001:end,:);
vy   = train_y(1:10000,:);
ty = train_y(10001:end,:);

rand('state',0);
nn  = nnsetup(X);     
opts.plot = plot;

start_time = clock;
nn = nntrain(nn, tx, ty, opts, vx, vy);  %  nntrain takes validation set as last two arguments (optionally)

[er, bad] = nntest(nn, test_x, test_y);
end_time = clock;
predict_time = etime(end_time,start_time)/100;
end

