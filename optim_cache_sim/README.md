# Cache Simulator

A two-level cache simulator.

The simulator reads the provided trace files (``trace1.txt`` and ``trace2.txt``) and simulates a cache’s operations. Finally, it will calculate the miss rate, access time (latency), replacement count and the AMAT (average memory average time).

The simulator is implemented as a two-level cache, *i.e.*, a L1 cache and another L2 cache. It also simulates a main memory which is the container of all original data. The configurations of L1, L2 and main memory are listed as follows: 

|  Level   | Capacity | Associativity | Block Size |        Write Policy         | Bus Latency | Hit Latency |
| :------: | :------: | :-----------: | :--------: | :-------------------------: | :---------: | :---------: |
| L1 Cache |   32KB   |       8       |    64B     | write back + write allocate |      0      |  3 cycles   |
| L2 Cache |  256KB   |       8       |    64B     | write back + write allocate |  6 cycles   |  4 cycles   |
|  Memory  |   inf    |       -       |     -      |              -              |      0      | 100 cycles  |

The frequency of CPU and memory are set to 2Ghz. The replacement policy of cache is set to workset policy. L1 and L2 caches use prefetch algorithm to reduce cache misses.

## Files

``main.cc``: program entrance. It parses the command-line arguments, does the setup of the cache, and reads the trace files.

``cache.h``: definition of cache.

``cache.cc``: implementation of cache functions which are defined in ``cache.h``

``memory.h``: definition of lower memory.

``memory.cc``: calculates the memory latency for each time memory is accessed.

``storage.h``: the upper class of both cache and memory.

``def.h``: definition of some constants.

``01-mcf-gem5-xcg.trace``: the first trace file, with each line representing a cache access.

``02-stream-gem5-xaa.trace``: the second trace file, with each line representing a cache access.

## Usage

1. Compile the C++ code by typing ``make`` in your terminal.

2. Run the code:

   ```bash
   ./sim -i 01-mcf-gem5-xcg.trace -r 50 -l50 -p4 -s16 
   ```

   Required arguments:

   ```bash
   -i			the input trace file
   -r			the number of rounds that you want to run
   ```

   Optional arguments:

   ```bash
   -l			   workset length (default: 50)
   -p			   the number of prefetch lines (default: 4)
   -s			   prefetch stride (default: 16)
   --no_optim	disable the optimization metrics used in this simulator, i.e., workset policy and prefetch. Instead, LRU will be used as
               replacement policy. This option can be used to testify the efficacy of these optimization metrics.
   ```

   P.S. This simulator uses C function ``getopt`` to read command-line arguments. Therefore, in order to input optional short arguments, *i.e.*, ``-l``, ``-p`` or ``-s``, you should omit the blank space between the argument label and its value. That is to say, type ``-l50`` instead of ``-l 50``.

