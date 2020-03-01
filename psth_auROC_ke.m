function [PSTH_auROC] = psth_auROC_ke(baseline, test)
%This function computes a normalized psth using the area under the
%operating characteristic curve. This results in a normalized psth with
%possible values ranging from -1 to 1, with 0 being no difference between baseline and test.
%Here the baseline actually is the other trial type. For example, the
%baseline can be Right trial, the test the left trial;The timebine should
%be the same between the baseline and test


%The baseline can either be input to the function (for events such as taste
%events in which the baseline is taken from outside the psth window.
% if ~isempty(varargin)
%     baseline = varargin{1};
% else
%This creates a matrix of the baseline activity
% baseline = psth(:,1:baseline_time*1000/bin);
% end

totalbins = size(baseline,2); 

%this finds the maximum value of the full psth to be used 
maxval = max(max([baseline; test]));

%This sets the number of increments for the critical value to take for the
%Sliding_Probability
increments = 100;

%This computes the probability vector that the baseline activity is greater
%than the criterion value which ranges from 0 to maxval.
for j = 1:totalbins
    Baseline_Prob = Sliding_Probability(baseline(:,j), maxval, increments);
    
    Baseline_Prob = fliplr(Baseline_Prob);
    
    %Baseline_Prob(end + 1) = 1;
    
    Baseline_Prob = [0 Baseline_Prob 1];

%This computes the normalized psth for each bin of the psth matrix


    %this computes the auROC for the currrent bin
    CurrentBin_Prob = fliplr(Sliding_Probability(test(:,j), maxval, increments));
    
    CurrentBin_Prob = [0 CurrentBin_Prob 1]; %#ok<AGROW>
   
    %This uses a trapezoidal approximation to compute the integral. This
    %is the method used by the perfcurve function which computes ROC stats.
    PSTH_auROC(j) = (trapz(Baseline_Prob,CurrentBin_Prob)-0.5)*2;% scale auROC from 0-1 to -1 to 1;
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Prob_Vector] = Sliding_Probability(vector, maxval, increments)
%This computes a sliding probability of whether the vector values are
%greater than a criterion value which moves from zero to maxval in
%increments of the input variable, increments
%   Detailed explanation goes here

critvalues = 0:maxval/increments:maxval;

if maxval == 0
    
    Prob_Vector = zeros(1, length(critvalues));
    return
end

for i = 1:length(critvalues)
    
    Prob_Vector(i) = sum(vector > critvalues(i))/length(vector);
    
end


end
end

