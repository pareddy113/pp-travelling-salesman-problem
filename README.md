# pp-travelling-salesman-problem
Use C and MPI to compute solution the Travelling Salesman Problem using parallel search with branch and bound algorithm.


**Input:**
usa13509.tsp from: [URL](http://elib.zib.de/pub/mp−testdata/tsp/tsplib/tsp/index.html)


**Terminology:**
Travelling Salesman problem
NP Hard
Branch & Bound search algorithm
Held-Karp algorithm


**RUN:**

**_.pbs_** files are the job files that are submitted to the cluster where we specify the number of nodes, name of the executable file, and all other configurations for the job to run.

**_.pl_** files are the script files that automate the whole process by executing the c program and submitting the .pbs job to the cluster.

**_.err_** files give any error after executing the job.

**_.out_** files are the output files after the job being executed.

Run the .pl file to start the program execution on the server.
