#ifndef CACHE_DEF_H_
#define CACHE_DEF_H_

#define TRUE 1
#define FALSE 0
#define INF 0x3f3f3f3f
#define L1_CACHE 1
#define L2_CACHE 2

const double cpu_freq = 2e9;
const int m_hit_latency = 100;
const int m_bus_latency = 0;
const int l1_hit_latency = 3;
const int l1_bus_latency = 0;
const int l2_hit_latency = 4;
const int l2_bus_latency = 6;

const int l1_size = 32 * 1024;
const int l2_size = 256 * 1024;
const int l1_asso = 8;
const int l2_asso = 8;
const int l1_block = 64;
const int l2_block = 64;
// use write-back method
const int l1_write_through = FALSE;
const int l2_write_through = FALSE;
const int l1_write_allocate = TRUE;
const int l2_write_allocate = TRUE;

#endif //CACHE_DEF_H_
