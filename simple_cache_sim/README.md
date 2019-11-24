# Cache Simulator

A single-level cache simulator.

The simulator reads the provided trace files (``trace1.txt`` and ``trace2.txt``), simulates a cache’s operations, and calculate hits and misses.

The simulator also simulates a memory as the lower level of cache. Memory latency is taken into consideration and is set to 100 cycles. The frequency of CPU and memory are set to 2Ghz.

The replacement policy of cache is set to LRU.

## Files

``main.cc``: program entrance. It parses the command-line arguments, does the setup of the cache, and reads the trace files.

``cache.h``: definition of cache.

``cache.cc``: implementation of cache functions which are defined in ``cache.h``

``memory.h``: definition of lower memory.

``memory.cc``: calculates the memory latency for each time memory is accessed.

``storage.h``: the upper class of both cache and memory.

``def.h``: definition of some constants.

``trace1.txt``: the first trace file, with each line representing a cache access.

``trace2.txt``: the second trace file, with each line representing a cache access.

## Usage

1. Compile the C++ code by typing ``make`` in your terminal.

2. Run the code:

   ```bash
   ./sim -i trace1.txt -c 16 -s 4 -l 1 -b 8 --write_through
   ```

   Required arguments:

   ```bash
   -i		the input trace file
   -c		cache size (in KB)
   -s		associativity
   -b		block size
   -l		cache hit latency
   ```

   Optional arguments:

   ```bash
   --write_through		use "write through" to replace "write back"
   --write_allocate	enable "write allocate"
   --debug				verbose mode, print all intermediate results (you may need to 					  redirect the outputs to a text file)
   ```

   

