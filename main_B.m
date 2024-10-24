% This script simulates network B over multiple trials, varying evidence qualities and phase differences.

clear all

num_trials = 3000;             % Number of trials
fs = 1000;                     % Sampling frequency
decision_boundary1 = 0.15;     % Decision threshold
a_values = 0.1:0.025:0.35;     % Evidence qualities (input amplitudes)
phi = 0:pi/4:pi;               % Phase differences for the sine wave inputs

% Loop over different evidence qualities and phase differences
for i = 1:length(a_values)
    for j = 1:length(phi)
        % Simulation of network B with current parameters
        [acc1(i, j), decision_reaction_times1, error_reaction_times1, decision_trial_indices1, error_trial_indices1, non_decision_trial_indices, avg_decision_rts1(i, j), avg_error_rts1(i, j), x2_thrr{i, j}, x17_thrr{i, j}, x2_three{i, j}, x17_three{i, j}, stc_rtc1(i, j), std_error_rts1(i, j), oo1_e1, oo2_e1] = fsm_B(num_trials, fs, decision_boundary1, a_values(i), 0.1, 12, phi(j));
    end
end

% End of simulation script for network B
