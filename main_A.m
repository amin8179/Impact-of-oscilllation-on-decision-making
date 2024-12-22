% This script simulates network A over multiple trials, varying evidence qualities and phase differences.

clear all
num_trials=3000;   % Number of trials
fs=1000;           % Sampling frequency

decision_boundary=0.026;  % Decision threshold

a_values1 = [0.6:0.1/10:0.8];    % Amplitude of stimulus inputs to column 1
a_values2 = [0.6:(-0.1)/10:0.4]; % Amplitude of stimulus inputs to column 2

ev=(a_values1-a_values2)./(a_values1+a_values2);  % Evidence qualities (input amplitudes)


% Loop over different evidence qualities and phase differences

for i=1:length(a_values1)

% Simulation of network A
[acc(i),decision_reaction_times,error_reaction_times,decision_trial_indices,error_trial_indices,non_decision_trial_indices{i},avg_decision_rts(i),avg_error_rts(i),x2_thr{i},x17_thr{i},x2_thre{i},x17_thre{i},stc_rtc(i),std_error_rts(i), o1_e1,o2_e1,input_1,input_2]=fsm_A(num_trials,fs,decision_boundary,a_values1(i),0.07,12,0,a_values2(i));

end

% End of simulation script for network A
