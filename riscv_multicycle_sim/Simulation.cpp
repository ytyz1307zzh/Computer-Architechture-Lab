/*
 @Date  : 11/6/2019
 @Author: Zhihan Zhang
 @mail  : zhangzhihan@pku.edu.cn
*/

#include <string.h>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include "Simulation.h"
using namespace std;


#define add 0
#define sub 1
#define mul 2
#define div 3
#define _xor 4
#define _and 5
#define _or 6
#define _not 7
#define addw 8
#define mulw 9
#define subw 10
#define divw 11
#define modw 12
#define sll 16
#define srl 17
#define mulh 18
#define lt 19
#define mod 20
#define sra 21
#define beq 22
#define bneq 23
#define ge 24
#define upPC 25
#define addiw 26
#define lui 27
#define auipc 28

// jump choices
#define no_jump 0
#define imm_jump 1
#define imm_reg_jump 2

// ALU source operators
#define alu_reg_reg 0
#define alu_reg_imm 1
#define alu_pc_4 2
#define alu_pc_imm 3

// Memory read choices
#define no_read 0
#define read_byte 1
#define read_half 2
#define read_word 3
#define read_double 4

//Memory write choices
#define no_write 0
#define write_byte 1
#define write_half 2
#define write_word 3
#define write_double 4

Elf64_Ehdr elf64_hdr;

// count instruction numbers
long long instr_num=0;
long long cycle_num = 0;

// whether to use single-step mode
bool single_step = false;

//除法与取余
bool div_flag = false;
int div_rs1 = 0, div_rs2 = 0;
int div_rdq = 0, div_rdr = 0;


void read_elf(char *filename);
void read_elf_header();
void read_elf_sections();
void read_symtable();
void read_program();

void Run();
void Fetch();
void Decode();
void Execute();
void Memory();
void Writeback();

//get the immediate of SB-type instructions
unsigned int BEQ_imm(unsigned int x)
{
    unsigned int imm12 = (x >> 31) & 1;
    unsigned int imm11 = (x >> 7) & 1;
    unsigned int imm1_4 = (x >> 8) & 0xf;
    unsigned int imm5_10 = (x >> 25) & 0x3f;
    return (imm12 << 12) | (imm11 << 11) | (imm5_10 << 5) | (imm1_4 << 1);
}

//get the immediate of UJ-type instructions
unsigned int JAL_imm(unsigned int x)
{
    unsigned int imm20 = (x >> 31) & 1;
    unsigned int imm11 = (x >> 20) & 1;
    unsigned int imm1_10 = (x >> 21) & 0x3ff;
    unsigned int imm12_19 = (x >> 12) & 0xff;
    return (imm20 << 20) | (imm12_19 << 12) | (imm11 << 11) | (imm1_10 << 1);
}

//get the immediate of S-type instructions
unsigned int SW_imm(unsigned int x)
{
    unsigned int imm0_4 = (x >> 7) & 0x1f;
	unsigned int imm5_11 = (x >> 25) & 0x7f;
    return (imm5_11 << 5) | imm0_4;
}

//获取指定位
unsigned int getbit(unsigned int instr, int start, int end)
{
    assert (start < end);
    unsigned int instr_bak = instr >> start;
    int bits = end - start;
    unsigned int tail = (1 << (bits + 1)) - 1;
    return instr_bak & tail;
}

//符号扩展
long long ext_signed(unsigned int src, int bits)
{
    if (bits > 0) {
        int head = 64 - bits;
        return ((long long)src << head) >> head;
    }
    else return (long long)src;
}

// print all register information
void print_reg()
{
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 4; j++)
            cout << reg_name[i*4+j] << ": " << reg[i*4+j] << "\t\t";
        cout << endl; 
    }
}

void print_mem()
{
    unsigned long long addr;
    cout << "Input a valid memory address: 0x";
    cin >> hex >> addr;
    if (addr >= 0 && addr < MAX - 8) {
        cout << "8 Bytes from address " << hex << addr << ": ";
        cout << "0x" << hex << setw(8) << setfill('0') << *(unsigned long long*)(memory + addr);
        cout << endl;
    }
    else cout << "Invalid address." << endl;
}

void read_elf_header()
{
	fread(&elf64_hdr, 1, sizeof(elf64_hdr), file);
	pnum = (int)elf64_hdr.e_phnum;
	padr = elf64_hdr.e_phoff;
	psize = (int)elf64_hdr.e_phentsize;
	entry = elf64_hdr.e_entry;
	snum = (int)elf64_hdr.e_shnum;
	sadr = elf64_hdr.e_shoff;
	ssize = (int)elf64_hdr.e_shentsize;
    shstrindex = (int)elf64_hdr.e_shstrndx;
}


void read_elf_sections()
{
    Elf64_Shdr* section_header = new Elf64_Shdr[snum];
    fseek(file, sadr, SEEK_SET);
    fread(section_header, sizeof(Elf64_Shdr), snum, file);
    for (int j = 0; j < snum; j++)
        if (j == shstrindex) 
            shstradr = section_header[j].sh_offset;
    for (int i = 0; i < snum; i++) {
        unsigned long long addr = shstradr + section_header[i].sh_name;
        char section_name[50] = {0};
        fseek(file, addr, SEEK_SET);
        while(char next_char = fgetc(file)) 
            section_name[strlen(section_name)] = next_char;
        if (strcmp(section_name, ".strtab") == 0)
            stradr = section_header[i].sh_addr + section_header[i].sh_offset;
        if (strcmp(section_name, ".symtab") == 0) {
            symadr = section_header[i].sh_offset;
            symnum = section_header[i].sh_size / section_header[i].sh_entsize;
            symsize = section_header[i].sh_entsize;
        }
    }
}

void read_program()
{
    Elf64_Phdr* program_header = new Elf64_Phdr[pnum];
    fseek(file, padr, SEEK_SET);
    fread(program_header, pnum, sizeof(Elf64_Phdr), file);
	for(int i = 0; i < pnum; i++) {
        fseek(file, program_header[i].p_offset, SEEK_SET);
        fread(memory + program_header[i].p_vaddr, 1, program_header[i].p_filesz, file);
	}
}


void read_symtable()
{

    Elf64_Sym* symtable = new Elf64_Sym[symnum];
    fseek(file, symadr, SEEK_SET);
    fread(symtable, symnum, sizeof(Elf64_Sym), file);
    for(int i = 0; i < symnum; i++)
	{
        char symbol_name[50] = {0};
        unsigned int addr = stradr + symtable[i].st_name;
        fseek(file, addr, SEEK_SET);

        while(char next_char = fgetc(file)) 
            symbol_name[strlen(symbol_name)] = next_char;

		if (strcmp(symbol_name, "main") == 0) {
            PC = symtable[i].st_value;
            endPC = PC + ((symtable[i].st_size) / 4) * 4;
		}
		if (strcmp(symbol_name, "__global_pointer$") == 0)
            gp = symtable[i].st_value;
    }
}

void read_elf(char *filename)
{
	file = fopen(filename, "r");
	read_elf_header();
	read_elf_sections();
    read_program();
	read_symtable();
}

void Finish()
{
    bool finish = false;
    while (!finish) {
        cout << endl;
        cout << "----------------------------------" << endl;
        cout << "Input r to print register values." << endl;
        cout << "Input m to print memory values." << endl;
        cout << "Input e to exit." << endl;
        char exit_choice;
        cin >> exit_choice;
        switch(exit_choice) {
            case 'r':
                print_reg();
                break;
            case 'm':
                print_mem();
                break;
            case 'e':
                cout << "Exiting..." << endl;
                finish = true;
                break;
        }
    }
}


int main(int argc, char** argv)
{
	assert (argc == 2);

    cout << "Run in single-step mode? (y/N) ";
    char step_choice[5] = {0};
    cin.getline(step_choice, 5);
    if (strlen(step_choice) != 0 && step_choice[0] == 'y')
        single_step = true;

    char* elf_filename = argv[1];
	read_elf(elf_filename);

	reg[3] = gp;
	reg[2] = MAX / 2;

	Run();

    cout << "Finished." << endl;
	cout << "Instruction Numbers: " << dec << instr_num << endl;
    cout << "Cycle Numbers: " << dec << cycle_num << endl;
    double cpi = (double)cycle_num / instr_num;
    cout << "CPI: " << cpi << endl;

    Finish();

	return 0;
}



void Run()
{
    
	while(PC != endPC - 4)
	{
		Fetch();
		Decode();
		Execute();
		Memory();
		Writeback();

         // single-step mode
        if (single_step) {
            bool Continue = false; // continue to next instr?
            while (!Continue) {
                Continue = false;
                cout << endl;
                cout << "------------------------------------------" << endl;
                cout << "Input c to stop single-stepping & continue." << endl;
                cout << "Input n to next step." << endl;
                cout << "Input r to print register values." << endl;
                cout << "Input m to print memory values." << endl;
                char cont_choice;
                cin >> cont_choice;
                switch (cont_choice) {
                    case 'c':
                        single_step = false;
                        Continue = true;
                        break;
                    case 'r':
                        print_reg();
                        break;
                    case 'm':
                        print_mem();
                        break;
                    case 'n':
                        Continue = true;
                        break;
                }
            }
        }

        cout << endl;
        reg[0] = 0;
        instr_num++;

	}
}

void Fetch()
{
    cout << hex << "Current PC: 0x" << PC << endl;
	unsigned int instr = *(unsigned int*)(memory + PC);
    cout << hex << "Current Instr: 0x" << instr << endl;
    FET_DEC.PC = PC;
	FET_DEC.instr = instr;
	PC = PC + 4;
	cycle_num++;
}

void Decode()
{
	unsigned int instr=FET_DEC.instr;
	int EXTbits=0;
	unsigned int EXTsrc=0;

	int ALUop, ALUSrc = 0;
    int branch = no_jump; // no jump
    int MemRead = no_read;
    int MemWrite = no_write;
    bool RegWrite = true;
    bool MemtoReg = false;

    string instr_name;

    OP = getbit(instr, 0, 6);
	rd = getbit(instr, 7, 11);
	func3 = getbit(instr, 12, 14);
	func7 = getbit(instr, 25, 31);
	
    rs1 = getbit(instr, 15, 19);
    rs2 = getbit(instr, 20, 24);


	if(OP==OP_R) {

		EXTbits = 0;
		ALUSrc = alu_reg_reg;

		if(func3 == 0) {
            if (func7 == 0) {
                ALUop = add;
                instr_name = "add";
            }
            if (func7 == 1) {
                ALUop = mul;
                instr_name = "mul";
            }
            if (func7 == 0x20) {
                ALUop = sub;
                instr_name = "sub";
            }
		}

		if(func3 == 1) {
            if (func7 == 0) {
                ALUop = sll;
                instr_name = "sll";
            }
            if (func7 == 1) {
                ALUop = mulh;
                instr_name = "mulh";
            }
		}

		if(func3 == 2) {
            ALUop = lt;
            instr_name = "slt";
		}

		if(func3 == 4) {
            if (func7 == 0) {
                ALUop = _xor;
                instr_name = "xor";
            }
            if (func7 == 1) {
                ALUop = div;
                instr_name = "div";
            }
		}

		if(func3 == 5) {
            if (func7 == 0) {
                ALUop = srl;
                instr_name = "srl";
            }
            if (func7 == 0x20) {
                ALUop = sra;
                instr_name = "sra";
            }
		}

		if(func3 == 6) {
            if (func7 == 0) {
                ALUop = _or;
                instr_name = "_or";
            }
            if (func7 == 1) {
                ALUop = mod;
                instr_name = "rem";
            }
		}

		if(func3 == 7) {
            ALUop = _and;
            instr_name = "and";
		}
	}

	if(OP==OP_I) {
        shamt = getbit(instr, 20, 25); // shift bits
        EXTsrc = getbit(instr, 20, 31);
        func7 = getbit(instr, 26, 31);
        EXTbits = 12;
		ALUSrc = alu_reg_imm;

        if (func3 == 0) {
            ALUop = add;
            instr_name = "addi";
        }
        
        if (func3 == 1) {
            EXTsrc = shamt;
            EXTbits = 6;
            ALUop = sll;
            instr_name = "slli";
        }

        if (func3 == 2) {
            ALUop = lt;
            instr_name = "slti";
        }

        if (func3 == 4) {
            ALUop = _xor;
            instr_name = "xori";
        }

        if (func3 == 5) {
            EXTbits = 6;
            if (func7 == 0) {
                ALUop = srl;
                instr_name = "srli";
                EXTsrc = shamt;
            }

            if (func7 == 0x10) {
                ALUop = sra;
                instr_name = "srai";
                EXTsrc = shamt;
            }
        }

        if (func3 == 6) {
            ALUop = _or;
            instr_name = "ori";
        }

        if (func3 == 7) {
            ALUop = _and;
            instr_name = "andi";
        }
    }

    if(OP==OP_SW){

        EXTbits = 12;
        EXTsrc = SW_imm(instr);
		ALUSrc = alu_reg_imm;
		RegWrite = false;

        if (func3 == 0) {
            MemWrite = write_byte;
            ALUop = add;
            instr_name = "sb";
        }

        if (func3 == 1) {
            MemWrite = write_half;
            ALUop = add;
            instr_name = "sh";
        }
        
        if (func3 == 2) {
            MemWrite = write_word;
            ALUop = add;
            instr_name = "sw";
        }

        if (func3 == 3) {
            MemWrite = write_double;
            ALUop = add;
            instr_name = "sd";
        }
    }

    if(OP==OP_LW) {

        EXTbits = 12;
        EXTsrc = getbit(instr, 20, 31);
		ALUSrc = alu_reg_imm;
		MemtoReg = true;

        if (func3 == 0) {
            MemRead = read_byte;
            ALUop = add;
            instr_name = "lb";
        }

        if (func3 == 1) {
            MemRead = read_half;
            ALUop = add;
            instr_name = "lh";
        }

        if (func3 == 2) {
            MemRead = read_word;
            ALUop = add;
            instr_name = "lw";
        }

        if (func3 == 3) {
            MemRead = read_double;
            ALUop = add;
            instr_name = "ld";
        }
    }

    if(OP==OP_BEQ) {

        EXTbits = 13;
        EXTsrc = BEQ_imm(instr);
		ALUSrc = alu_reg_reg;
		branch = imm_jump; // imm jump
		RegWrite = false;

        if (func3 == 0) {
            ALUop = beq;
            instr_name = "beq";
        }
        
        if (func3 == 1) {
            ALUop = bneq;
            instr_name = "bneq";
        }

        if (func3 == 4) {
            ALUop = lt;
            instr_name = "blt";
        }

        if (func3 == 5) {
            ALUop = ge;
            instr_name = "bge";
        }
    }

    if(OP==OP_JALR) {

        EXTbits = 12;
        EXTsrc = getbit(instr, 20, 31);
		ALUSrc = alu_pc_4;
		branch = imm_reg_jump; // imm + reg jump
		ALUop = add;
		instr_name = "jalr";
    }

    if (OP == OP_JAL) {

        EXTbits = 20;
        EXTsrc = JAL_imm(instr);
		ALUSrc = alu_pc_4;
		branch = imm_jump;
		ALUop = add;
		instr_name = "jal";
    }

    if (OP == OP_IW) {
        EXTbits = 12;
        shamt = getbit(instr, 20, 25);
        EXTsrc = getbit(instr, 20, 31);
		ALUSrc = alu_reg_imm;
        ALUop = addiw;
        instr_name = "addiw";
    }

    if (OP == OP_LUI) {
        EXTbits = 0;
        EXTsrc = getbit(instr, 12, 31);
        EXTsrc = (EXTsrc << 12);
		ALUSrc = alu_reg_imm;
		ALUop = lui;
		instr_name = "lui";
    }

    if (OP == OP_AUIPC) {
        EXTbits = 0;
        EXTsrc = getbit(instr, 12, 31);
        EXTsrc = (EXTsrc << 20);
		ALUSrc = alu_pc_imm;
		ALUop = auipc;
		instr_name = "auipc";
    }

    if (OP == OP_RW) {
        EXTbits = 0;
		ALUSrc = alu_reg_reg;

        if (func3 == 0) {
            if (func7 == 0) {
                ALUop = addw;
                instr_name = "addw";
            }
            
            if (func7 == 1) {
                ALUop = mulw;
                instr_name = "mulw";
            }

            if (func7 == 0x20) {
                ALUop = subw;
                instr_name = "subw";
            }
        }

        if (func3 == 4 && func7 == 1) {
            ALUop = divw;
            instr_name = "divw";
        }
        if (func3 == 6 && func7 == 1) {
            ALUop = modw;
            instr_name = "remw";
        }
    }

    cout << "Instruction Name:" + instr_name << endl;
    DEC_EX.PC = FET_DEC.PC;
    DEC_EX.Rd = rd;
    DEC_EX.rs1 = rs1;
    DEC_EX.rs2 = rs2;
    DEC_EX.Imm = ext_signed(EXTsrc, EXTbits);
    DEC_EX.Ctrl_EX_ALUOp = ALUop;
    DEC_EX.Ctrl_EX_ALUSrc = ALUSrc;
    DEC_EX.Ctrl_EX_RegDst = rd;
    DEC_EX.Ctrl_M_Branch = branch;
    DEC_EX.Ctrl_M_MemRead = MemRead;
    DEC_EX.Ctrl_M_MemWrite = MemWrite;
    DEC_EX.Ctrl_WB_MemtoReg = MemtoReg;
    DEC_EX.Ctrl_WB_RegWrite = RegWrite;
    cycle_num++;
}

void Execute()
{
	//read DEC_EX
	long long nextPC = DEC_EX.PC;
	int ALUOp = DEC_EX.Ctrl_EX_ALUOp;
	int ALUSrc = DEC_EX.Ctrl_EX_ALUSrc;
	int branch = DEC_EX.Ctrl_M_Branch;
	long long imm = DEC_EX.Imm;
	int exe_cycles = 1; // default execution cycles

    if (ALUOp != mod || ALUOp != modw)
        div_flag = false;

    // 处理跳转
	if (branch == imm_jump)
        nextPC += DEC_EX.Imm;
    if (branch == imm_reg_jump) {
        nextPC = reg[DEC_EX.rs1] + DEC_EX.Imm;
        nextPC = nextPC - (nextPC % 2);
	}

    // ALU source operators
	long long op1 = 0, op2 = 0;
    // Branch: jump or not
	bool jump = false;

    if (ALUSrc == alu_reg_reg) {
        op1 = reg[rs1];
        op2 = reg[rs2];
    }
    if (ALUSrc == alu_reg_imm) {
        op1 = reg[rs1];
        op2 = imm;
    }
    if (ALUSrc == alu_pc_4) {
        op1 = DEC_EX.PC;
        op2 = 4;
        jump = true;
    }
    if (ALUSrc == alu_pc_imm) {
        op1 = DEC_EX.PC;
        op2 = imm;
    }

	REG ALUout;

    switch(ALUOp)
    {
    case add:
        ALUout = op1 + op2;
        break;
    case sub:
        ALUout = op1 - op2;
        break;
    case mul:
        ALUout = op1 * op2;
        exe_cycles = 2;
        break;
    case div:
        ALUout = op1 / op2;
        exe_cycles = 40;
        div_flag = true;
        div_rs1 = rs1;
        div_rs2 = rs2;
        div_rdq = rd;
        break;
    case _xor:
        ALUout = op1 ^ op2;
        break;
    case _and:
        ALUout = op1 & op2;
        break;
    case _or:
        ALUout = op1 | op2;
        break;
    case sll:
        ALUout = op1 << op2;
        break;
    case srl:
        ALUout = ((unsigned long long)op1) >> op2;
        break;
    case mulh:
    {
        long long rs1_high = (reg[rs1] >> 32) & 0xffffffff;
        long long rs1_low = reg[rs1] & 0xffffffff;
        long long rs2_high = (reg[rs2] >> 32) & 0xffffffff;
        long long rs2_low = reg[rs2] & 0xffffffff;

        long long h1h2 = rs1_high * rs2_high;
        long long h1l2 = ((rs1_high * rs2_low) >> 32) & 0xffffffff;
        long long l1h2 = ((rs2_high * rs1_low) >> 32) & 0xffffffff;
        ALUout = h1h2 + h1l2 + l1h2;
        exe_cycles = 2;
        break;
    }
    case lt:
        if (op1 < op2) jump = true;
        break;
    case mod:
        ALUout = op1 % op2;
        exe_cycles = 40;
        break;
    case sra:
        ALUout = ((long long)op1) >> op2;
        break;
    case beq:
        if (op1 == op2) jump = true;
        break;
    case bneq:
        if (op1 != op2) jump = true;
        break;
    case ge:
        if (op1 >= op2) jump = true;
        break;
    case addiw:
        ALUout = (long long)((int)(op1 + (long long)op2));
        break;
    case lui:
        ALUout = (long long)op2;
        break;
    case auipc:
        ALUout = op1 + (long long)op2;
        break;
    case addw:
        ALUout = (long long)((int)op1 + (int)op2);
        break;
    case mulw:
        ALUout = (long long)((int)op1 * (int)op2);
        break;
    case subw:
        ALUout = (long long)((int)op1 - (int)op2);
        break;
    case divw:
        ALUout = (long long)((int)op1 / (int)op2);
        exe_cycles = 40;
        div_flag = true;
        div_rs1 = rs1;
        div_rs2 = rs2;
        div_rdq = rd;
        break;
    case modw:
        ALUout = (long long)((int)op1 % (int)op2);
        exe_cycles = 40;
	}

    EX_MEM.ALU_out = ALUout;
    EX_MEM.Ctrl_M_Branch = branch;
    EX_MEM.Ctrl_M_MemRead = DEC_EX.Ctrl_M_MemRead;
    EX_MEM.Ctrl_M_MemWrite = DEC_EX.Ctrl_M_MemWrite;
    EX_MEM.Ctrl_WB_MemtoReg = DEC_EX.Ctrl_WB_MemtoReg;
    EX_MEM.Ctrl_WB_RegWrite = DEC_EX.Ctrl_WB_RegWrite;
    EX_MEM.PC = nextPC;
    EX_MEM.RegDst = DEC_EX.Ctrl_EX_RegDst;
    EX_MEM.rt = DEC_EX.rs2;
    EX_MEM.jump = jump;
    cycle_num += exe_cycles;
}

void Memory()
{
    int MemRead = EX_MEM.Ctrl_M_MemRead;
    int MemWrite = EX_MEM.Ctrl_M_MemWrite;
    int jump = EX_MEM.jump;
    int rs2 = EX_MEM.rt;
    long long ALUout = EX_MEM.ALU_out;
    long long read_data = 0;

	//complete Branch instruction PC change
    if (jump) PC = EX_MEM.PC;

    if (MemRead || MemWrite) cycle_num++;

    // Memory Read
    if (MemRead == read_byte) read_data = *(memory + ALUout);
    if (MemRead == read_half) read_data = *(short *)(memory + ALUout);
    if (MemRead == read_word) read_data = *(int *)(memory + ALUout);
    if (MemRead == read_double) read_data = *(long long *)(memory + ALUout);

    // Memory Write
    if (MemWrite == write_byte)
        *(unsigned char *)(memory + ALUout) = (unsigned char)reg[rs2];
    if (MemWrite == write_half)
        *(unsigned short *)(memory + ALUout) = (unsigned short)reg[rs2];
    if (MemWrite == write_word)
        *(unsigned int *)(memory + ALUout) = (unsigned int)reg[rs2];
    if (MemWrite == write_double)
        *(unsigned long long *)(memory + ALUout) = (unsigned long long)reg[rs2];

    MEM_WB.ALU_out = ALUout;
    MEM_WB.Ctrl_WB_MemtoReg = EX_MEM.Ctrl_WB_MemtoReg;
    MEM_WB.Ctrl_WB_RegWrite = EX_MEM.Ctrl_WB_RegWrite;
    MEM_WB.read_data = read_data;
    MEM_WB.RegDst = EX_MEM.RegDst;
}

void Writeback()
{
    unsigned long long ALUout = MEM_WB.ALU_out;
    bool MemtoReg = MEM_WB.Ctrl_WB_MemtoReg;
    bool RegWrite = MEM_WB.Ctrl_WB_RegWrite;
    long long read_data = MEM_WB.read_data;
    int rd = MEM_WB.RegDst;

    // register write
	if (RegWrite) {
        cycle_num++;
        if (MemtoReg)
            reg[rd] = read_data;
        else 
            reg[rd] = ALUout;
	}
}
