# RISC-V Simulator

A simple simulator for the RV64I instruction set of RISC-V ISA.

The simulator reads a executable file (elf file), and simulates the execution process of RV64I (e.g., registers and memory)

The simulator allows **single-step mode**, i.e., executing the instructions one by one, and you can check the status of the simulated registers and memory whenever you want during the execution process.

The code has been tested on ``g++`` 7.4 and ``riscv64-unknown-elf-gcc`` 9.2 (the RISC-V toolchain).

All scripts are written based on the template given by Computer Architecture Lab course, PKU.

## Files

``Simulation.cpp``: The main script, including all core functions.

``Simulation.h``: The header file, including some global variables and constants.

``Simulation``: The executable file of ``Simulation.cpp``, compiled by g++.

``quicksort.c``: The code to test the simulator.

``quicksort``: The executable file of ``quicksort.c``, compiled by RISC-V toolchain.

``RISCV-simple-greencard.pdf``: The instructions of RISC-V RV64I ISA, used as the reference document of the simulator.

## Usage

#### 1. Install the RISC-V toolchain

You can follow the instructions in this [repo](https://github.com/riscv/riscv-gnu-toolchain).

#### 2. Compile the test code

```bash
riscv64-unknown-elf-gcc -Wa,-march=rv64i -o quicksort quicksort.c
```

Note: You can directly use the provided executable file ``quicksort``. It is compiled using the exact command as above. By doing so you can skip the long process of installing RISC-V toolchain.

#### 3. Compile the simulator

```bash
g++ Simulation.cpp -o Simulation
```

#### 4. Run the simulator

```bash
./Simulation quicksort
```

