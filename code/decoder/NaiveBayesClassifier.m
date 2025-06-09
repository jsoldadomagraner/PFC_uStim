% Fit Naive Bayes Classifier

function [Yhat,Err,p] = NaiveBayesClassifier(X,Y,NBmethod,CVmethod,folds, ...
    binWidth)

[~,trials] = size(X);
class      = unique(Y);
c          = length(class);

% remove low FR channels (more adequate for NBP)
Xfr  = mean(X,2)/(binWidth/1000);
idx  = Xfr<2;    % spikes/s.
X    = X(~idx,:);

[N,~] = size(X);

switch CVmethod
    case 'Kfold'
        cvpart = cvpartition(Y,"KFold",folds);
    case 'LOOCV'
        folds  = trials;
        cvpart = cvpartition(trials,"Leaveout");
end

Yhat = zeros(length(Y));
Err  = zeros(1,folds);
p    = zeros(length(Y),c);

for f = 1:folds
    
    idxtrain = cvpart.training(f);
    idxtest  = cvpart.test(f);
    
    Xtrain = X(:,idxtrain);
    Xtest  = X(:,idxtest);
    
    Ytrain = Y(idxtrain);
    Ytest  = Y(idxtest);
    
    testn = length(Ytest);
    Lpost = zeros(testn,c);
    
    switch NBmethod
        case 'Gaussian'
            
            S = zeros(N,N);
            m = cell(1,c);
            for i=1:c
                % Sample mean and cov from training data for each class
                idxtraini = Ytrain==i;
                m{i} = mean(Xtrain(:,idxtraini),2);
                Si   = cov(Xtrain(:,idxtraini)');
                ni   = length(idxtraini);
                S    = S + (ni/testn)*Si;
            end
            
            for i=1:c
                % Log posterior of the test data for each class
                for tr = 1:testn
                    x = Xtest(:,tr);
                    Lpost(tr,i) = log(1/c) + m{i}'*(S\x) - 0.5* m{i}'*(S\m{i});
                end
            end              
            
        case 'Poisson'
            
            for i=1:c
                % Sample mean from training data for each class
                idxtraini = Ytrain==i;
                mi = mean(Xtrain(:,idxtraini),2);
                
                % Log posterior of the test data for each class
                for tr = 1:testn
                    x = Xtest(:,tr);
                    Lpost(tr,i) = log(1/c) + ...
                        nansum(x.*log(mi)-mi-log(factorial(x)));
                end
            end          
    end
    
    % test set label predictions
    [~,Yhatf]  = max(Lpost,[],2);
    Yhat(idxtest) = Yhatf;
    
    % log probability for each class
    p(idxtest,:) = Lpost;
    
    % Compute error fraction
    Err(f) = sum((Ytest-Yhat(idxtest))~=0)/testn;
    
end

end


