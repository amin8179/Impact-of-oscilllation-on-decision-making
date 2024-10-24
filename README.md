**README for Neural Decision-Making Model Simulations**


Overview
This project simulates neural decision-making processes using two different network models (Network A and Network B). Each model is simulated over multiple trials, varying evidence qualities (input stimulus amplitudes) and phase differences between inputs. The scripts employ Euler's method to solve systems of differential equations that represent the dynamics of neuron populations. The goal is to simulate decision and error trials and analyze reaction times, accuracy, and other neural outputs.

File Descriptions

1. fsm_A.m 
This script simulates Network A. The function fsm_A runs a neural decision-making model over multiple trials, simulating the activity of two interconnected neural populations. The inputs include Gaussian bumps and sine waves with various phase shifts. Key outputs include decision and error reaction times, accuracy, and other relevant neural signals.

Inputs:
num_trials: Number of trials for the simulation.
fs: Sampling frequency.
decision_boundary: Threshold for decision-making.
a1: Amplitude for the first Gaussian bump (input stimulus).
a2: Amplitude for the sine wave inputs (external stimuli).
f2: Frequency for the sine wave inputs.
phi: Phase shift for one of the sine wave inputs.
Outputs:
Various metrics like reaction times, decision indices, average reaction times, accuracy, etc.


2. fsm_B.m 
This script simulates Network B. It follows a similar structure to fsm_A.m, with differences in parameters and boundary conditions. It calculates decision-making performance metrics like reaction times, accuracy, and neural output.

Inputs:
Same structure as fsm_A.m.
Outputs:
Similar to fsm_A.m, providing reaction times, decision indices, and accuracy metrics.


3. main_A.m 
This script is used to run simulations for Network A over multiple trials. It loops over different evidence qualities (input amplitudes) and phase differences to explore how these parameters affect decision-making in the model.

Key parameters:
num_trials = 3000: Number of trials per simulation.
fs = 1000: Sampling frequency.
decision_boundary1 = 0.005: Threshold for decision-making.
a_values = 0.1:0.025:0.35: Range of evidence quality values.
phi = 0:pi/4:pi: Phase differences for the sine wave inputs.


4. main_B.m 
This script runs simulations for Network B. Like main_A.m, it loops over different evidence qualities and phase differences but applies to the second model. It helps assess the impact of changing parameters on decision-making.

Key parameters:
Similar structure to main_A.m but for Network B.




5. model33.m 
This file contains the model33 function, which defines the system of differential equations that represent the dynamics of neural populations. The equations describe how each neuron’s activity evolves based on external inputs, coupling between neurons, and internal noise.

Key Components:
State variables (x1, x2, etc.): Represent different neural populations.
Inputs (bump1, bump2, y1, y2): External inputs to the model (Gaussian bumps and sine waves).
Noise: Random noise affecting the neural dynamics.
Differential Equations: Represent the change in neural activity over time.




6. movingAverage.m (​(movingAverage))
This file defines the movingAverage function, which calculates the moving average of the input data over a specified window size. It smooths the data for easier trend analysis during classification of trials.

Inputs:

data: The input data vector.
windowSize: Size of the moving average window.
Outputs:

mav: The moving average of the input data.



How to Run the Simulations
Ensure all files are in the same directory.
Edit simulation parameters in main_A.m or main_B.m to fit your needs. Modify variables like num_trials, fs, a_values, and phi to explore different setups.
Run main_A.m or main_B.m from the MATLAB environment to execute the simulations. The script will automatically call fsm_A.m or fsm_B.m to handle the trial-by-trial simulations.
Analyze outputs: After running, the reaction times, accuracy metrics, and other outputs will be stored in the specified variables for further analysis.



Dependencies

MATLAB R2018b or later.
The model33 and movingAverage functions must be present in the working directory.



Notes
Ensure that the functions model33.m and movingAverage.m are available, as they are essential for simulating neural activity and processing the outputs.
Adjust the decision_boundary and noise parameters to suit different decision-making scenarios.
The moving average window size (k) is set to 100 by default but can be adjusted to better fit your data.
By using these files, you can simulate how neural populations respond to varying inputs, analyze decision and error trials, and explore neural decision-making mechanisms.











