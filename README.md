Gabe Appleton

# Python Process Pools on Multiple Hosts

## Abstract

The goal for this project is to reimplement something like Python's `multiprocessing.Pool` API on top of OpenMPI. This would enable significantly faster prototyping of parallel code. It would also enable for easier porting of existing parallel code to a cluster environment. As the local interface would have the same types of built-in restrictions that an OpenMPI reimplementation would have, it is possible to develop a mostly-drop-in replacement that enables you to quickly scale to an entire cluster. The only thing that would need changing is the synchronization primitives, and every effort will be taken to preserve the existing API where possible.

## Implementation Strategies

The goal is essentially as follows. Across N systems, N+1 total processes will be created. The first is a coordinator process on the master node, and the remainder are worker daemons. The master node is responsible for remembering the state of job requests, resubmitting them if necessary, and spawning worker processes to satisfy job requests. In the event of a non-master daemon failure, the system should be expected to continue by resubmitting the lost jobs to another node. The master node will intentionally prioritize its own daemon when possible to ensure that during low-outstanding-job periods computation is kept as local as possible.

Note that there is some chance that OpenMPI will not end up being suitable for this task, but I believe it will be.

Shared objects are divided into two categories. Read-only shared objects will be managed by the coordinator daemon, who will keep an ID'd cache of objects, as well as keep track of which job is using them. If no active job is using them for some configurable time, the object is purged from cache. The master node also has the option to invalidate a cache entry by updating the corresponding object.

Shared writable objects will require some degree of mutual exclusion. It may be that OpenMPI primitives are sufficient for this purpose. If not, the such synchronization primitives will be implemented via the daemons keeping a Hybrid Logical Clock or similar to enforce causal consistency.

## Benchmarking Proposal

Benchmarks will be done by comparing two metrics against the current `multiprocessing.Pool` implementation. The first is naive wall-clock time, and the second is an approximation of core time. It is unclear at the moment whether this will be by estimation or by observation.

Each API will be compared using various sizes of the Heat Equation model from Section 36.3. This is chosen in part because it requires some form of shared object support. Since this is a particular weak-point of the current API, and is a somewhat large hurdle for cross-host computation, this will be a good comparison point. Improvements beyond the cost of communication between hosts will indicate real success.

## Anticipated Challenges

### Python State

The biggest challenge will likely be managing Python state between hosts. Not only can I not assume that each host has the same set of libraries installed, but it is not obvious how to determine what libraries are needed for a given job. I can think of two possible workarounds, but further research is required. The workarounds are as follows:

1. If a job fails due to an ImportError, determine the missing library and attempt to install it as a user, then restart the job
2. If the above fails or is not feasible, rerun job on the master node and raise a warning in the logs
3. If the above fails, kill the program and return an error message asking for an administrator to check the packages

Additionally, one could require a list of packages to be submitted for checking on each daemon node. This may be best left optional, however.

### Serialization Overhead

This is a problem shared with the `multiprocessing.Pool` API, but the fact that objects need to be serialized in order to be transferred is going to cause overhead problems. It is possible that workarounds will be needed in order to accomodate this. For instance, one could imagine a small C++ extension that would grab the interpreter lock and serialize all the objects in parallel. This would mean not going with Python's built-in serialization tools, however, and it is possible that one of the several Python OpenMPI bindings effectively solves this.
