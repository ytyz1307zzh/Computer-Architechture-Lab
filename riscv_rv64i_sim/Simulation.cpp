/*
 @Date  : 10/23/2019
 @Author: Zhihan Zhang
 @mail  : zhangzhihan@pku.edu.cn
*/

#include <string.h>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include "Simulation.h"
using namespace std;

Elf64_Ehdr elf64_hdr;

extern unsigned int cadr;
extern unsigned int csize;
extern unsigned int vadr;
extern unsigned long long gp;
extern unsigned int madr;
extern unsigned long long endPC;
extern unsigned int entry;

void Execute();

// count instruction numbers
long long instr_num=0;

// whether to use single-step mode
bool single_step = false;

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

void read_Elf_header()
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
        fread(memory + program_header[i].p_vaddr, 1, program_header[i].p_memsz, file);
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

        // find PCs and gp
        while(char next_char = fgetc(file)) 
            symbol_name[strlen(symbol_name)] = next_char;

        // 初始PC和终止PC
		if (strcmp(symbol_name, "main") == 0) {
            PC = symtable[i].st_value;
            endPC = PC + ((symtable[i].st_size) / 4) * 4;
		}
        // 全局数据段地址
		if (strcmp(symbol_name, "__global_pointer$") == 0)
            gp = symtable[i].st_value;
    }
}

void read_elf(char* filename)
{
	file = fopen(filename, "r");
	read_Elf_header();
	read_elf_sections();
    read_program();
	read_symtable();
}

void Run()
{
    cout << "Program begins from 0x" << hex << PC << endl;
    cout << "Program ends at 0x" << hex << endPC << endl;

	while(PC < endPC - 4) {
		//运行
		Execute();

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

	//设置全局数据段地址寄存器
	reg[3] = gp;

	reg[2] = MAX / 2; //栈基址 （sp寄存器）

	Run();

    cout << "Finished." << endl;
	cout << "Instruct Num: " << dec << instr_num << endl;

    Finish();

	return 0;
}


// execute a single instruction -- core function
void Execute()
{
    cout << hex << "Current PC: 0x" << PC << endl;
	unsigned int instr = *(unsigned int*)(memory + PC);
    cout << hex << "Current Instr: 0x" << instr << endl;

    // arguments for sign extension
	int EXTbits = 0;
	unsigned int EXTsrc = 0;
    // whether the PC need to change (otherwise +4)
    bool change_PC = false;
    // the name of this instruction, in string
    string instr_name;

    OP = getbit(instr, 0, 6);
	rd = getbit(instr, 7, 11);
	func3 = getbit(instr, 12, 14);
	func7 = getbit(instr, 25, 31);
	func7_s = getbit(instr, 26, 31);
    rs1 = getbit(instr, 15, 19);
    rs2 = getbit(instr, 20, 24);

    // R-type instructions
	if(OP == OP_R) {

		if(func3 == 0) {
            switch(func7) {
            case 0:
                instr_name = "add";
                reg[rd] = reg[rs1] + reg[rs2];
                break;
            case 1:
                instr_name = "mul";
                reg[rd] = reg[rs1] * reg[rs2];
                break;
            case 32:
                instr_name = "sub";
                reg[rd] = reg[rs1] - reg[rs2];
                break;
            default:
                cout << "Invalid Instruction" << endl;
                exit(0);
            }
		}

		if (func3 == 1) {
            switch (func7) {
            case 0:
                instr_name = "sll";
                reg[rd] = reg[rs1] << reg[rs2];
                break;
            case 1: {
                    instr_name = "mulh";
                    long long rs1_high = (reg[rs1] >> 32) & 0xffffffff;
                    long long rs1_low = reg[rs1] & 0xffffffff;
                    long long rs2_high = (reg[rs2] >> 32) & 0xffffffff;
                    long long rs2_low = reg[rs2] & 0xffffffff;

                    long long h1h2 = rs1_high * rs2_high;
                    long long h1l2 = ((rs1_high * rs2_low) >> 32) & 0xffffffff;
                    long long l1h2 = ((rs2_high * rs1_low) >> 32) & 0xffffffff;
                    reg[rd] = h1h2 + h1l2 + l1h2;
                    break;
                }
            default:
                cout << "Invalid Instruction" << endl;
                exit(0);
            }
		}

		if(func3 == 2) {
            instr_name = "slt";
            reg[rd] = (reg[rs1] < reg[rs2]) ? 1 : 0;
		}

		if(func3 == 4) {
            switch(func7) {
            case 0:
                instr_name = "xor";
                reg[rd] = reg[rs1] ^ reg[rs2];
                break;
            case 1:
                instr_name = "div";
                reg[rd] = reg[rs1] / reg[rs2];
                break;
            default:
                cout << "Invalid Instruction" << endl;
                exit(0);
            }
		}

		if(func3 == 5) {
            switch(func7) {
            case 0:
                instr_name = "srl";
                reg[rd] = (unsigned long long)reg[rs1] >> reg[rs2];
                break;
            case 32:
                instr_name = "sra";
                reg[rd] = (long long)reg[rs1] >> reg[rs2];
                break;
            default:
                cout << "Invalid Instruction" << endl;
                exit(0);
            }
		}

		if(func3 == 6) {
            switch(func7) {
            case 0:
                instr_name = "or";
                reg[rd] = reg[rs1] | reg[rs2];
                break;
            case 1:
                instr_name = "rem";
                reg[rd] = reg[rs1] % reg[rs2];
                break;
            default:
                cout << "Invalid Instruction" << endl;
                exit(0);
            }
		}

		if(func3 == 7) {
            instr_name = "and";
            reg[rd] = reg[rs1] & reg[rs2];
		}
	}

	else if(OP==OP_I) {
        EXTsrc = getbit(instr, 20, 31);
        EXTbits = 12;
        long long imm = ext_signed(EXTsrc, EXTbits);

        switch(func3) {
        case 0:
            instr_name = "addi";
            reg[rd] = reg[rs1] + imm;
            break;
        case 1:
            instr_name = "slli";
            reg[rd] = reg[rs1] << (imm & 0x3f);
            break;
        case 2:
            instr_name = "slti";
            reg[rd] = (reg[rs1] < imm) ? 1 : 0;
            break;
        case 4:
            instr_name = "xori";
            reg[rd] = reg[rs1] ^ imm;
            break;
        case 5:
            if (func7_s == 0) {
                instr_name = "srli";
                reg[rd] = (unsigned long long)reg[rs1] >> (imm & 0x3f);
                break;
            }
            else {
                instr_name = "srai";
                reg[rd] = (long long)reg[rs1] >> (imm & 0x3f);
                break;
            }
        case 6:
            instr_name = "ori";
            reg[rd] = reg[rs1] | imm;
            break;
        case 7:
            instr_name = "andi";
            reg[rd] = reg[rs1] & imm;
            break;
        default:
            cout << "Invalid Instruction" << endl;
            exit(0);
        }
    }
    else if(OP==OP_SW)
    {
        EXTbits = 12;
        EXTsrc = SW_imm(instr);
        long long offset = ext_signed(EXTsrc, EXTbits);

        switch(func3) {
        case 0:
            instr_name = "sb";
            *(memory + reg[rs1] + offset) = (char)reg[rs2];
            break;
        case 1:
            instr_name = "sh";
            *(short *)(memory + reg[rs1] + offset) = (short)reg[rs2];
            break;
        case 2:
            instr_name = "sw";
            *(int *)(memory + reg[rs1] + offset) = (int)reg[rs2];
            break;
        case 3:
            instr_name = "sd";
            *(long long *)(memory + reg[rs1] + offset) = reg[rs2];
            break;
        default:
            cout << "Invalid Instruction" << endl;
            exit(0);
        }
    }
    else if(OP==OP_LW)
    {
        EXTbits = 12;
        EXTsrc = getbit(instr, 20, 31);
        long long offset = ext_signed(EXTsrc, EXTbits);

        switch(func3) {
        case 0:
            instr_name = "lb";
            reg[rd] = ext_signed(*(memory + reg[rs1] + offset), 8);
            break;
        case 1:
            instr_name = "lh";
            reg[rd] = ext_signed(*(short *)(memory + reg[rs1] + offset), 16);
            break;
        case 2:
            instr_name = "lw";
            reg[rd] = ext_signed(*(int *)(memory + reg[rs1] + offset), 32);
            break;
        case 3:
            instr_name = "ld";
            reg[rd] = *(long long *)(memory + reg[rs1] + offset);
            break;
        default:
            cout << "Invalid Instruction" << endl;
            exit(0);
        }
    }
    else if(OP==OP_BEQ)
    {
        EXTbits = 13;
        EXTsrc = BEQ_imm(instr);
        long long offset = ext_signed(EXTsrc, EXTbits);

        switch(func3) {
        case 0:
            instr_name = "beq";
            if (reg[rs1] == reg[rs2]) {
                PC += offset;
                change_PC = true;
            }
            break;
        case 1:
            instr_name = "bne";
            if (reg[rs1] != reg[rs2]) {
                PC += offset;
                change_PC = true;
            }
            break;
        case 4:
            instr_name = "blt";
            if (reg[rs1] < reg[rs2]) {
                PC += offset;
                change_PC = true;
            }
            break;
        case 5:
            instr_name = "bge";
            if (reg[rs1] >= reg[rs2]) {
                PC += offset;
                change_PC = true;
            }
            break;
        default:
            cout << "Invalid Instruction" << endl;
            exit(0);
        }
    }

    else if(OP==OP_JALR) {
        EXTbits = 12;
        EXTsrc = getbit(instr, 20, 31);
        long long imm = ext_signed(EXTsrc, EXTbits);
		instr_name = "jalr";
        reg[rd] = PC + 4;
        PC = (reg[rs1] + imm) & ~1;
        change_PC = true;
    }

    else if (OP == OP_JAL) {
        EXTbits = 20;
        EXTsrc = JAL_imm(instr);
        long long imm = ext_signed(EXTsrc, EXTbits);
		instr_name = "jal";
        reg[rd] = PC + 4;
        PC += imm;
        change_PC = true;
    }

    else if (OP == OP_IW) {
        EXTbits = 12;
        EXTsrc = getbit(instr, 20, 31);
        long long imm = ext_signed(EXTsrc, EXTbits);
		instr_name = "addiw";
        reg[rd] = (long long)((int)reg[rs1] + (int)imm);
    }
    
    else if (OP == OP_LUI) {
        EXTbits = 0;
        EXTsrc = getbit(instr, 12, 31);
        EXTsrc = (EXTsrc << 12);
        long long offset = ext_signed(EXTsrc, EXTbits);
		instr_name = "lui";
        reg[rd] = offset;
    }

    else if (OP == OP_AUIPC) {
        EXTbits = 0;
        EXTsrc = getbit(instr, 12, 31);
        EXTsrc = (EXTsrc << 20);
        long long offset = ext_signed(EXTsrc, EXTbits);
		instr_name = "auipc";
        reg[rd] = PC + offset;
    }

    else if (OP == OP_RW) {

        if (func3 == F3_ADDW) {
            switch(func7) {
            case 0:
                instr_name = "addw";
                reg[rd] = (long long)((int)reg[rs1] + (int)reg[rs2]);
                break;
            case 1:
                instr_name = "mulw";
                reg[rd] = (long long)((int)reg[rs1] * (int)reg[rs2]);
                break;
            case 32:
                instr_name = "subw";
                reg[rd] = (long long)((int)reg[rs1] - (int)reg[rs2]);
                break;
            default:
                cout << "Invalid Instruction" << endl;
                exit(0);
            }
        }

        if (func3 == 4 & func7 == 1) {
            instr_name = "divw";
            reg[rd] = (long long)((int)reg[rs1] / (int)reg[rs2]);
        }
        
        if (func3 == 6 and func7 == 1) {
            instr_name = "remw";
            reg[rd] = (long long)((int)reg[rs1] % (int)reg[rs2]);
        }
    }

    cout << "Instr Name:"+ instr_name << endl;

    // did not jump, go to the next instruction
    if (!change_PC)
        PC += 4;
}
