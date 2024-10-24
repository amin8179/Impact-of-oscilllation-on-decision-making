% This script simulates network A over multiple trials, varying evidence qualities and phase differences.

clear all

num_trials = 3000;             % Number of trials
fs = 1000;                     % Sampling frequency
decision_boundary = 0.005;     % Decision threshold
a_values = 0.1:0.025:0.35;     % Evidence qualities (input amplitudes)
phi = 0:pi/4:pi;               % Phase differences for the sine wave inputs

% Loop over different evidence qualities and phase differences
for i = 1:length(a_values)
    for j = 1:length(phi)
        % Simulation of network B with current parameters
        [acc(i, j), decision_reaction_times, error_reaction_times, decision_trial_indices, error_trial_indices, non_decision_trial_indices, avg_decision_rts(i, j), avg_error_rts(i, j), x2_thrr{i, j}, x17_thrr{i, j}, x2_three{i, j}, x17_three{i, j}, stc_rtc(i, j), std_error_rts(i, j), oo1_e1, oo2_e1] = fsm_A(num_trials, fs, decision_boundary, a_values(i), 0.1, 12, phi(j));
    end
end

% End of simulation script for network A
