NetWalk (based on the Random walk with restarts algorithm proposed by Kohler, et al., 2008)
iteratively simulates random transitions of a walker from a node to a randomly selected neighbour
node and where at any time step the walk can be restarted depending on a predefined probability.
Random walk with restarts is slightly different than PageRank with priors in the way that it
normalizes the link weights. The convergence is decided by either having a probability difference
less than 10e-6 between two consecutive time steps or achieving the limit of the number of
iterations, set to 50 (though in practice less than 20 iterations are typically sufficient to satisfy
the first criterion). There is no method specific parameter for NetWalk.
Therefore an example run for NetWalk is as follows:

R --slave --args test_proteins.txt test_interactions.txt output.txt < random_walk.r

R --slave --args nodes.txt edges.txt output.txt < random_walk.r

