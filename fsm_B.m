% Simulation of network B
function [acc,decision_reaction_times,error_reaction_times,decision_trial_indicesa,error_trial_indicesa,non_decision_trial_indicesa,avg_decision_rts,avg_error_rts,x2_thr,x17_thr,x2_thre,x17_thre,std_decision_rts,std_error_rts, o1_e1,o2_e1,input_1,input_2]=fsm_B(num_trials,fs,decision_boundary,a1,a2,f2,phi,a3)
  % This function simulates a neural decision-making model over multiple trials.
    % Inputs:
    %   num_trials           - Number of simulation trials
    %   fs                   - Sampling frequency
    %   decision_boundary    - Threshold for decision-making
    %   a1                   - Amplitude of the Gaussian bump (input stimulus) to column 1
    %   a2                   - Amplitude for the sine wave inputs (external stimuli)
    %   a3                   - Amplitude of the Gaussian bump (input stimulus) to column 2
    %   f2                   - Frequency for the sine wave inputs
    %   phi                  - Phase shift for one of the sine wave inputs
    %
    % Outputs:
    %   acc                          - Accuracy percentage across trials
    %   decision_reaction_times      - Reaction times for decision trials
    %   error_reaction_times         - Reaction times for error trials
    %   decision_trial_indicesa      - Indices of decision trials
    %   error_trial_indicesa         - Indices of error trials
    %   non_decision_trial_indicesa  - Indices of non-decision trials
    %   avg_decision_rts             - Average reaction time for decision trials
    %   avg_error_rts                - Average reaction time for error trials
    %   x2_thr, x17_thr              - Threshold values at reaction time for decision trials
    %   x2_thre, x17_thre            - Threshold values at reaction time for error trials
    %   std_decision_rts             - Standard deviation of decision reaction times
    %   std_error_rts                - Standard deviation of error reaction times
    %   o1_e1, o2_e1                 - Outputs from specific neurons in the model
    %   input_1, input_2             - Input stimuli used in the simulation

    Ae = 3.25;           % Excitatory synaptic gain
    Ai = 22;             % Inhibitory synaptic gain
    a1e = 10;           % 1/Time constant for excitatory population
    a1i = 100;             % 1/Time constant for inhibitory population
    cL21 = 0.85;            % Coupling from population 2 to 1
    cL12 = 0.85;            % Coupling from population 1 to 2
    cLi21 = 50;          % Inhibitory coupling from 2 to 1
    cLi12 = 50;          % Inhibitory coupling from 1 to 2
    c1 = 0.15;            % Coupling parameters (intrinsic connections)
    c2 = 0.1;
    c3 = 32;
    c4 = 55;
    c5 = 50;
    c6 = 0.1;
    c7 = 31.5;
    noise_amp = 2;      % Amplitude of noise





    tampa = 5/1000;      %Time constant
    f = 3;               %End time of stimulation    
decision_boundaries1=0.08;     % Decision boundary for classification
  
  
    % Time vector
    t = linspace(0, 4, 4*fs);

    % Parameters for Gaussian bumps
    b = 1;               % Start time of stimulation 
    c = 0.005;           % Width of the Gaussian bump

    X0 = zeros(1, 30);   % Initial state vector for the model (30 state variables)

    % Time step for Euler's method
    dt = t(2) - t(1);

    % Initialize variables to store reaction times and accuracies
    avg_decision_rts = [];
    avg_error_rts = [];
    std_decision_rts = [];
    std_error_rts = [];

    % Generate external sine wave inputs with phase shift
    y1 = a2 .* sin(2 .* pi .* f2 .* t);  % Sine wave inputs to column 1
    y2 = a2 .* sin(2 .* pi .* f2 .* t+ phi);       % Sine wave inputs to column 2

    % Loop over different 'a' values 
  a_values = (a1);
  a1_values = (a3);


    for a_idx = 1:length(a_values)
        a = a_values(a_idx);  % Amplitude for the first Gaussian bump

        % Construct the first Gaussian bump (stimulus) for population 1
        rise1 = a * exp(-((t - b).^2) / (2 * c^2));
        rise1(t > b) = 0;
        sustain1 = zeros(size(t));
        sustain1((t >= b) & (t <= f)) = a;
        decay1 = a * exp(-((t - f).^2) / (2 * c^2));
        decay1(t < f) = 0;
        bump1 = max([rise1; sustain1; decay1]);

        % Construct the second Gaussian bump (stimulus) for population 2
        amplitude_bump2 = a1_values(a_idx);
        rise2 = amplitude_bump2 * exp(-((t - b).^2) / (2 * c^2));
        rise2(t > b) = 0;
        sustain2 = zeros(size(t));
        sustain2((t >= b) & (t <= f)) = amplitude_bump2;
        decay2 = amplitude_bump2 * exp(-((t - f).^2) / (2 * c^2));
        decay2(t < f) = 0;
        bump2 = max([rise2; sustain2; decay2]);

        % Initialize lists to store trial indices and decision times
        decision_trial_indices = [];
        error_trial_indices = [];
        non_decision_trial_indices = [];
        decision_reaction_times = [];
        error_reaction_times = [];
        x2_thr = [];
        x17_thr = [];
        x2_thre = [];
        x17_thre = [];

        % Run simulations for the specified number of trials
        for trial = 1:num_trials
            rng(trial);  % Seed the random number generator for reproducibility

            % Initialize the state variable matrix for this trial
            X = zeros(length(t), length(X0));
            X(1, :) = X0;

            % Simulate the model using Euler's method
            for i = 1:length(t) - 1
                X(i + 1, :) = X(i, :) + dt * model33(X(i, :), i, bump1, bump2, y1, y2, noise_amp, tampa, Ae, Ai, a1e, a1i, cL21, cL12, cLi21, cLi12, c1, c2, c3, c4, c5, c6, c7);
            end

            % Store inputs and outputs for analysis
            input_1(trial, :) = bump1 + X(:, 29)' + y1;
            input_2(trial, :) = bump2 + X(:, 30)' + y2;
            o1_e1(trial, :) = X(:, 1)';       % Output from excitatory neuron e1 in population 1
            o1_e2(trial, :) = X(:, 2)';       % Output from excitatory neuron e2 in population 1
            o1_i1(trial, :) = X(:, 3)';       % Output from inhibitory neuron i1 in population 1
            o2_e1(trial, :) = X(:, 15)';      % Output from excitatory neuron e1 in population 2
            o2_e2(trial, :) = X(:, 16)';      % Output from excitatory neuron e2 in population 2
            o2_i1(trial, :) = X(:, 17)';      % Output from inhibitory neuron i1 in population 2

            % Calculate Local Field Potentials (LFPs) and EEG signals (not used in further analysis)
            output1_LFP(trial, :) = abs((X(:, 1)' + X(:, 6)')) + abs((X(:, 2)' + X(:, 7)') - X(:, 3)') + abs(X(:, 4)' - X(:, 5)');
            output2_LFP(trial, :) = abs((X(:, 16)' + X(:, 21)')) + abs((X(:, 17)' + X(:, 22)') - X(:, 18)') + abs(X(:, 19)' - X(:, 20)');
            output1_EEG(trial, :) = X(:, 2)' - X(:, 3)' + X(:, 7)';
            output2_EEG(trial, :) = X(:, 17)' - X(:, 18)' + X(:, 22)';

        X1=X(:, 1)';  %EPSP of neuron e1 in column 1
        X2=X(:, 15)';  %EPSP of neuron e1 in column 2

        X1=X1';
        X2=X2';


        % Analyze segments of x2 and x17 within the specified time window
        start_idx = find(t >= 1, 1, 'first');
        end_idx = find(t <= f, 1, 'last');
        x2 = X1(start_idx:end_idx + 1,1);
        x17 = X2(start_idx:end_idx + 1,1);

        % Classify the trials

        % Classification logic as per Python code

            if max(x2) >= max(x17) && abs(max(abs(x2)) - max(abs(x17))) >= decision_boundaries1
                decision_trial_indices = [decision_trial_indices, trial];
            elseif max(x17) >= max(x2) && abs(max(abs(x17)) - max(abs(x2))) >= decision_boundaries1
                error_trial_indices = [error_trial_indices, trial];
            else
                non_decision_trial_indices = [non_decision_trial_indices, trial];
            end
        end

        % Store trial indices for output
        decision_trial_indicesa = decision_trial_indices;
        error_trial_indicesa = error_trial_indices;
        non_decision_trial_indicesa = non_decision_trial_indices;

        % Calculate decision times for correct trials
        for trial_idx = 1:length(decision_trial_indices)
            trial = decision_trial_indices(trial_idx);
            rng(trial);

            % Re-simulate the trial to get accurate decision times
            X = zeros(length(t), length(X0));
            X(1, :) = X0;
            for i = 1:length(t) - 1
                X(i + 1, :) = X(i, :) + dt * model33(X(i, :), i, bump1, bump2, y1, y2, noise_amp, tampa, Ae, Ai, a1e, a1i, cL21, cL12, cLi21, cLi12, c1, c2, c3, c4, c5, c6, c7);
            end
            x22 = X(:, 1);  %EPSP of neuron e1 in column 1
            x177 = X(:, 15);%EPSP of neuron e1 in column 2 

            x2 = x22;
            x17 = x177;

            % Find the decision time when x2 crosses the decision boundary
            rt_indices = find(x2 > decision_boundary, 1, 'first');
            if ~isempty(rt_indices)
                first_rt = t(rt_indices);
                x2_th = x2(rt_indices);
                x17_th = x17(rt_indices);

                decision_reaction_times = [decision_reaction_times, first_rt];
                x2_thr = [x2_thr, x2_th];
                x17_thr = [x17_thr, x17_th];
            end
        end

        % Calculate decision times for error trials
        for trial_idx = 1:length(error_trial_indices)
            trial = error_trial_indices(trial_idx);
            rng(trial);

            % Re-simulate the trial to get accurate decision times
            X = zeros(length(t), length(X0));
            X(1, :) = X0;
            for i = 1:length(t) - 1
                X(i + 1, :) = X(i, :) + dt * model33(X(i, :), i, bump1, bump2, y1, y2, noise_amp, tampa, Ae, Ai, a1e, a1i, cL21, cL12, cLi21, cLi12, c1, c2, c3, c4, c5, c6, c7);
            end
            x22 = X(:, 1);
            x177 = X(:, 15);

            x2 = x22;
            x17 = x177;

            % Find the decision time when x17 crosses the decision boundary
            rt_indices = find(x17 > decision_boundary, 1, 'first');
            if ~isempty(rt_indices)
                first_rt = t(rt_indices);
                x2_th = x2(rt_indices);
                x17_th = x17(rt_indices);

                error_reaction_times = [error_reaction_times, first_rt];
                x2_thre = [x2_thre, x2_th];
                x17_thre = [x17_thre, x17_th];
            end
        end

        % Calculate average decision times and standard deviations
        if ~isempty(decision_reaction_times)
            avg_decision_rt = mean(decision_reaction_times, 'omitnan');
            std_decision_rt = std(decision_reaction_times, 'omitnan') / sqrt(length(decision_reaction_times));
        else
            avg_decision_rt = NaN;
            std_decision_rt = NaN;
        end

        if ~isempty(error_reaction_times)
            avg_error_rt = mean(error_reaction_times, 'omitnan');
            std_error_rt = std(error_reaction_times, 'omitnan') / sqrt(length(error_reaction_times));
        else
            avg_error_rt = NaN;
            std_error_rt = NaN;
        end

        % Store the decision times and accuracies
        avg_decision_rts(a_idx) = avg_decision_rt;
        avg_error_rts(a_idx) = avg_error_rt;
        std_decision_rts(a_idx) = std_decision_rt;
        std_error_rts(a_idx) = std_error_rt;

        % Calculate accuracy percentage
        acc(a_idx) = (length(decision_trial_indices)) / (length(decision_trial_indices) + length(error_trial_indices)) * 100;
    end
end

% Note: The functions 'model33' is assumed to be defined elsewhere in your codebase.
