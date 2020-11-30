# Probabilistic_DynamicOutputFeedback

Like other code that you may have perused, I have done a decent job with commenting these scripts and functions.  You'll want to start with the main invoking script, which is:

# decentdynoutputGAnested_main.m

This work is detailed in Chapter 7 of my dissertation.  There, you can even find algorithm pseudo-code, which describes what all of this code is doing.  

Similar to what you may have seen before, I create the system models (Euler-Bernoulli beam finite element models) along with a nominal interconnection stiffness element.  I then form the composite system and evaluate its open-loop performance, and establish a performance reduction that I want the probabilistic robust controllers to accomplish.  All of this is done between lines 1 and 141.

On lines 143 to 155 I define my stochastic cost functional parameters and create a data structure for passing this on later.  

Lines 157 to 178 is done for computational efficiency.  During my optimization, I do not want the size function to be called repeatedly, as this would slow my code down.

initialpopulationloadingsupercomputer() is called on line 187.  You may want to take a look at this script that I wrote.  This script actually loads the results from the Loop_at_a_Time_DKIt repository that I had shared before.  The script creates my initial conditions for optimization: they are "seeds" from loop at a time synthesis via D/K iteration.

Now, we move into the function controllerdiagonalvectorization.m.  We will also move into controllervectorcompression.m.  Both of these functions and associated transformations were created for computational efficiency and so that the opimization variables (the controllers, in this case) could be put into a form that is required for genetic algorithm-based optimization.

# controllerdiagonalvectorization.m
Here, I transform the controllers into their complex modal form.  I then separate the real and complex parts from one another.  I also identify the controller poles that are purely real, as I must keep track of the real and complex eigenvectors as I am rearranging these systems.  Once again, if interested, I believe that I have some decent comments in this code if you want to peer under the hood to see my thought process.  Finally, we end up with the partitioned real/complex parts of the modal controller as it is passed to the main invoking script.  

# controllervectorcompression.m
This function simply removes redunancy in the mathematics.  Complex-conjugate eigenvalues appear in pairs --- this means that if I evolve the real or complex part of one pair, I am doing it to both of them.  Thus, for every congugate pair (along with their associated eigenvectors) I only need to evolve the positive or negative conjugate.  I then vectorize the entire controller.  That is, the controllers are in strictly-proper state variable form with some A,B, and C matrices, which I had worked to separate into real and complex portions previously.  At this step, we throw out the redundancy and vectorize the system as we ready ourselves for optimization.

This brings us to line 254 of decentdynoutputGAnested_main.m.  All problem data has been transformed for efficiency and into the form that is needed for parallel, stochastic optimization using genetic algorithms.  As you have worked through this script, you likely saw how I was packaging problem data into structures.  Now all of these structures are passed to passvariablestonestedfunction.m function.  This begins on line 269.

On line 326, I call the genetic algorithm function, which is a part of Matlab's Optimization Toolbox.  Just prior to this, though, you may want to see my custom gaoptions() script (line 320).  This is where I make important algorithm choices and also set the problem up so that it uses all available computer cores on a particular supercomputing node.  My call to ga involves an anonymous function call, @decentdynoutputfeedbackcostfcn_nested.  This function is defined by lines 330--445.  It may or may not be obvious, but this is a Monte Carlo evaluation that evaluates a combination of the probability of stability/performance (note that performance is a subset of stability).  The interconnection element is randomized at each step.  Thus, we culminate with a parallelized stochastic optimization problem that parameterizes mu synthesized controllers into vectors that represent their real/complex modal form for computational efficiency.  This step in the research was very successful, leading to the creation of a novel controller synthesis methods that yields probabilistically-robust, decentralized dynamic output feedback controllers for plants that possess random interconnection uncertainty.  As with all of my code, I have created it to work with plants of arbitrary dimension and input/output numbers.  










