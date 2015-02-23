%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metrics = checkPredictions(exp_data,model_data)
% Calculates all associated metrics to model predictions.
%
% INPUT:    exp_data and model_data: cell objects of mx2 (single KO) or
%           mx3 (double KO or epistasis) with the first row (or the first
%           two rows in the double KO or epistasis cases) with the gene
%           names, and the second (or third in the double KO or epistasis 
%           case) with the fitness value (Wx = mu_KOx/mu_WT in the single &
%           double KO cases, and Exy = Wxy-Wx*Wy in the epistasis case).
% OUTPUT:   metrics, a structure containing:
%           *CM:   Confusion Matrix ([TP FP;FN TN])
%           *sens: Sensitivity (TP/(TP+FN))
%           *spec: Specificity (TN/(TN+FP))
%           *accu: Accuracy ((TP+TN)/(TP+FP+FN+TN))
%           *PPV:  Positive Predictive Value (TP/(TP+FP))
%           *NPV:  Negative Predictive Value (TN/(TN+FN))
%           *ROC:  ROC curve (TP/(TP+FN) Vs. FP/(FP+TN))
%           *MCC:  Mathews correlation coefficient:
%                  (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
%           *GM:   Geometric mean (sqrt(sens*spec))
%           *F1:   F1-score (2*sens*spec/(sens+spec))
%
% Benjamín J. Sánchez. Last edited: 2014-11-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metrics = checkPredictions(exp_data,model_data)

%Calculation of TP, TN, FP & FN:
[m,n]         = size(model_data);
metrics.T     = zeros(100,1);
metrics.CM    = cell(100,1);
metrics.sens  = zeros(100,1);
metrics.spec  = zeros(100,1);
metrics.accu  = zeros(100,1);
metrics.PPV   = zeros(100,1);
metrics.NPV   = zeros(100,1);
metrics.ROC   = zeros(100,2);
metrics.MCC   = zeros(100,1);
metrics.GM    = zeros(100,1);
metrics.F1    = zeros(100,1);

for i = 1:100
    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;
    T  = i/100;  %Threshold
    for j = 1:m
        %If epistasis is evaluated, invert values:
        if n == 3
            model = -model_data{j,n};
            exp   = -exp_data{j,n};
        else
            model = model_data{j,n};
            exp   = exp_data{j,n};
        end
        %Update metrics:
        [TP,TN,FP,FN] = compare_prediction(model,exp,T,TP,TN,FP,FN);
    end
    
    %Calculate scores:
    metrics.T(i)     = T;
    metrics.CM{i,1}  = [TP FP;FN TN];
    metrics.sens(i)  = TP/(TP+FN);
    metrics.spec(i)  = TN/(TN+FP);
    metrics.accu(i)  = (TP+TN)/(TP+FP+FN+TN);
    metrics.PPV(i)   = TP/(TP+FP);
    metrics.NPV(i)   = TN/(TN+FN);
    metrics.ROC(i,:) = [FP/(FP+TN) TP/(TP+FN)];
    metrics.MCC(i)   = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    metrics.GM(i)    = sqrt(TP/(TP+FN)*TN/(TN+FP));
    metrics.F1(i)    = 2*TP/(TP+FN)*TN/(TN+FP)/(TP/(TP+FN)+TN/(TN+FP));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TP,TN,FP,FN] = compare_prediction(model,exp,T,TP,TN,FP,FN)
%Compares the model prediction with the experimental results and decides if
%it is a TP, TN, FP or FN. Returns all metrics updated.

if model >= T && exp >= 0.95
    TP = TP+1;
elseif model >= T && exp < 0.95
    FP = FP+1;
elseif model < T && exp >= 0.95
    FN = FN+1;
elseif model < T && exp < 0.95
    TN = TN+1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%