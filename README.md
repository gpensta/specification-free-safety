# Specification-Free-Safety

Objective: find a set W which always allow the system to recover. 

Example: pendulum.

# Install

You just need to add CORA 2021 toolbox in the path of your MatLab. 

https://tumcps.github.io/CORA/



Files:

    * pendulum.m: plot of sets. 
    * pendulum_evolution.m: plot of trajectories, sequence of inputs. 
    * heatmaps.m : plot of heatmaps of the state space.

# Bug 

In the file : pendulum_exact_integration.m, get_max_w() does not find the right W whereas it is OK for Euler: get_max_w_euler(), I do not know why yet... 
