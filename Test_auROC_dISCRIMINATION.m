%% TEST
% function adapted by Ke and Roberto from the original auRoc function used
% in the lab

function p = Test_auROC_dISCRIMINATION(A,B)


% This function is used to caluculate the p (probability value) that firing
% rate in the two vector A and B


%this finds the maximum value of the full psth (A) to be used 
maxval = max(max(A));


%This sets the number of increments for the critical value to take for the
%Sliding_Probability
increments = 1000;


Baseline_Prob = Sliding_Probability(B, maxval, increments); % vecor B (left goal)
Baseline_Prob = fliplr(Baseline_Prob);

Test_Prob = Sliding_Probability(A, maxval, increments);     % vecor A (right goal)
Test_Prob = fliplr(Test_Prob);

Baseline_Prob = [0 Baseline_Prob 1];
Test_Prob     = [0 Test_Prob 1];


PSTH_auROC = trapz(Baseline_Prob,Test_Prob); % problem if is negative


% mainen
p = 2*(PSTH_auROC-0.5);



%%
%figure(2)
%plot(Baseline_Prob,Test_Prob);xlim([0 1]);ylim([0 1]);hold on;plot([0 1],[0 1]);



