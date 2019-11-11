/*
 @Date  : 11/6/2019
 @Author: Zhihan Zhang
 @mail  : zhangzhihan@pku.edu.cn
*/

#include<stdint.h>

//simulation define
#define OP_JAL 0x6f
#define OP_R 51

#define OP_JALR 0x67
#define OP_LUI 55
#define OP_AUIPC 23

#define F3_ADD 0
#define F3_MUL 0

#define F7_MSE 1
#define F7_ADD 0

#define OP_I 0x13
#define F3_ADDI 0

#define OP_SW 0x23
#define F3_SB 0

#define OP_LW 0x03
#define F3_LB 0

#define OP_BEQ 0x63
#define F3_BEQ 0

#define OP_IW 0x1b
#define F3_ADDIW 0

#define OP_RW 0x3b
#define F3_ADDW 0
#define F7_ADDW 0

#define OP_SCALL 0x73
#define F3_SCALL 0
#define F7_SCALL 0

#define MAX 1<<28

typedef unsigned long long REG;
typedef unsigned short Elf64_Half;
typedef unsigned int Elf64_Word;
typedef	int  Elf64_Sword;
typedef unsigned long long Elf64_Xword;
typedef	long long  Elf64_Sxword;
typedef unsigned long long Elf64_Addr;
typedef unsigned long long Elf64_Off;
typedef unsigned short Elf64_Section;

#define	EI_CLASS 4
#define	EI_DATA 5
#define	EI_VERSION 6
#define	EI_OSABI 7
#define	EI_ABIVERSION 8
#define	EI_PAD 9
#define	EI_NIDENT 16

#define	SHN_UNDEF 0
#define	SHN_LOPROC 0xFF00
#define	SHN_HIPROC 0xFF1F
#define	SHN_LOOS 0xFF20
#define	SHN_HIOS 0xFF3F
#define	SHN_ABS 0xFFF1
#define	SHN_COMMON 0xFFF2

typedef struct
{
	unsigned char e_ident[16]; /* ELF identification */
	Elf64_Half e_type; /* Object file type */
	Elf64_Half e_machine; /* Machine type */
	Elf64_Word e_version; /* Object file version */
	Elf64_Addr e_entry; /* Entry point address */
	Elf64_Off e_phoff; /* Program header offset */
	Elf64_Off e_shoff; /* Section header offset */
	Elf64_Word e_flags; /* Processor-specific flags */
	Elf64_Half e_ehsize; /* ELF header size */
	Elf64_Half e_phentsize; /* Size of program header entry */
	Elf64_Half e_phnum; /* Number of program header entries */
	Elf64_Half e_shentsize; /* Size of section header entry */
	Elf64_Half e_shnum; /* Number of section header entries */
	Elf64_Half e_shstrndx; /* Section name string table index */
} Elf64_Ehdr;

typedef struct
{
	Elf64_Word sh_name; /* Section name */
	Elf64_Word sh_type; /* Section type */
	Elf64_Xword sh_flags; /* Section attributes */
	Elf64_Addr sh_addr; /* Virtual address in memory */
	Elf64_Off sh_offset; /* Offset in file */
	Elf64_Xword sh_size; /* Size of section */
	Elf64_Word sh_link; /* Link to other section */
	Elf64_Word sh_info; /* Miscellaneous information */
	Elf64_Xword sh_addralign; /* Address alignment boundary */
	Elf64_Xword sh_entsize; /* Size of entries, if section has table */
} Elf64_Shdr;

typedef struct
{
	Elf64_Word st_name; /* Symbol name */
	unsigned char st_info; /* Type and Binding attributes */
	unsigned char st_other; /* Reserved */
	Elf64_Half st_shndx; /* Section table index */
	Elf64_Addr st_value; /* Symbol value */
	Elf64_Xword st_size; /* Size of object (e.g., common) */
} Elf64_Sym;


typedef struct
{
	Elf64_Word p_type; /* Type of segment */
	Elf64_Word p_flags; /* Segment attributes */
	Elf64_Off p_offset; /* Offset in file */
	Elf64_Addr p_vaddr; /* Virtual address in memory */
	Elf64_Addr p_paddr; /* Reserved */
	Elf64_Xword p_filesz; /* Size of segment in file */
	Elf64_Xword p_memsz; /* Size of segment in memory */
	Elf64_Xword p_align; /* Alignment of segment */
} Elf64_Phdr;


//代码段在解释文件中的偏移地址
unsigned int cadr=0;

//代码段的长度
unsigned int csize=0;

//代码段在内存中的虚拟地址
unsigned int vadr=0;

//全局数据段在内存的地址
unsigned long long gp=0;

//main函数在内存中地址
unsigned int madr=0;

//程序结束时的PC
unsigned long long endPC=0;

//程序的入口地址
unsigned int entry=0;

FILE *file = NULL;

//寄存器堆
REG reg[32]={0};
//PC
long long PC=0;

//regNames
char reg_name[32][5] = {"r0", "ra", "sp", "gp",
						"tp", "t0", "t1", "t2",
						"s0", "s1", "a0", "a1",
						"a2", "a3", "a4", "a5",
						"a6", "a7", "s2", "s3",
						"s4", "s5", "s6", "s7",
						"s8", "s9", "s10", "s11",
						"t3", "t4", "t5", "t6"};

//Program headers
unsigned long long padr=0;
unsigned int psize=0;
unsigned int pnum=0;

//Section Headers
unsigned long long sadr=0;
unsigned int ssize=0;
unsigned int snum=0;

//Symbol table
unsigned int symnum=0;
unsigned long long symadr=0;
unsigned int symsize=0;

//用于指示 包含节名称的字符串是第几个节（从零开始计数）
unsigned int shstrindex=0;
unsigned int shstradr=0;

//字符串表在文件中地址，其内容包括.symtab和.debug节中的符号表
unsigned int stradr=0;

unsigned int OP=0;
unsigned int func3=0,func7=0, func7_s = 0;
int shamt=0;
int rs1=0,rs2=0,rd=0;

//主存
char memory[MAX] = {0};

struct FETDEC{
	unsigned int instr;
	long long PC;
}FET_DEC, FET_DEC_prev;


struct DECEX{
	int Rd;
	long long PC;
	long long Imm;
	int rs1, rs2;

	int Ctrl_EX_ALUSrc;
	int Ctrl_EX_ALUOp;
	int Ctrl_EX_RegDst;

	int Ctrl_M_Branch;
	int Ctrl_M_MemWrite;
	int Ctrl_M_MemRead;

	bool Ctrl_WB_RegWrite;
	bool Ctrl_WB_MemtoReg;

}DEC_EX, DEC_EX_prev;

struct EXMEM{
	long long PC;
	int RegDst;
	REG ALU_out;
	int jump;
	int rt;

	int Ctrl_M_Branch;
	int Ctrl_M_MemWrite;
	int Ctrl_M_MemRead;

	bool Ctrl_WB_RegWrite;
	bool Ctrl_WB_MemtoReg;

}EX_MEM, EX_MEM_prev;
struct MEMWB{
	long long read_data;
	REG ALU_out;
	int RegDst;

	bool Ctrl_WB_RegWrite;
	bool Ctrl_WB_MemtoReg;

}MEM_WB, MEM_WB_prev;
