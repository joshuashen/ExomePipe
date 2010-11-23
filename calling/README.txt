job.sh is the master script. Any parameters need to be specified in this script.
Presently indel calling is failing with Java Heap out of space error. I've allocated 1G of memory for the Heap. Allocating more than 1G of memory for the Heap is failing. Tried testing using qsub but the job is never triggered. Hence I am not able to test beyond indel calling. 
