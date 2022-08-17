function [er,predict_time] = EvaluateModelError(X,plot)
% This function is used to obtain the training error of NN

% INPUT: X is a size of 5*1 vector whose elements are the number of
% hidden_layers, number of neurous per hidden layer, dropout,
% log_learning_rate and logL2_weight_regularization, respectively
%  

%  X range:
%   options.xmin = [1 50  -10 -10 -10]
%   options.xmax = [3 500  0  0  0]

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

