# Csto
#Customizing the Simulation Environment

The simulation environment can be easily customized by modifying the parameters in the main script.

1. Changing the Start and Goal Locations

The start and goal positions of the path planning problem are defined in main.
By updating these variables, users can generate different map configurations and planning tasks:

start_point = [x_s, y_s, z_s];
end_point   = [x_g, y_g, z_g];


Different choices of start and goal points may lead to different polytope chains and homotopy classes.

2. Modifying Obstacle Configurations

The obstacle layout of the environment is specified by the obstacle data structure in obstacle data.
By changing the obstacle parameters (e.g., position, size, or shape), users can construct different simulation scenarios, including cluttered or narrow-passage environments.

After updating the obstacle data, the convex decomposition and polytope-chain construction are automatically recomputed.

3. Regenerating the Environment

Once the start/goal points or obstacle data are modified, simply rerun the main script.
The algorithm will automatically:

reconstruct the convex free-space decomposition,

search for a valid polytope chain,

and compute an obstacle-free trajectory using convex optimization.
