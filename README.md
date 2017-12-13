# electrokinetic_flow_one_electrolyte_drop_model
This repository provide a set of Matlab code for solving the one-electrolyte drop model.

## An Academic Project at New Jersey Institute of Technology

## Prerequistite

This project requires Matlab 2015R or later version to run the code for solving the problem.

### Code

The matlab code consists a set of reusable function files and four scripts. The scripts are used to solve for background fluid velocity with respect to different magnitudes of zeta potential and applied electric field in the time-dependent one-electrolyte drop model and in the steady state.

time_indep_lt0_problem.m is looking for an improved initial condition for initial guess of the steady state problem.
time_indep_gt0_problem.m is looking for a steady state solution of the system.
time_dep_lt0_problem.m is looking for an improved initial condition for running the simulation with t>0 successfully.
time_dep_gt0_problem.m is simulating the time dependent system with applied electric field for t>0.
 
### Run

You can run the script file by swithcing current editor to the specific script file and clicking the 'run' button, or type in the script's file name in the command window and press enter. 

### Adjusting parameters

The parameters for the problems can be adjusted within the corresponding script.
