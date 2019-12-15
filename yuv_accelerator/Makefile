CC = g++
SIMD = -msse2 -mmmx -mavx2
NOSIMD = -mno-sse2 -mno-mmx
OPTIM = -O2
PROC = yuv yuv_mmx yuv_sse yuv_avx


all: $(PROC)

yuv_mmx: 
	$(CC) -o yuv_mmx $(OPTIM) $(SIMD) yuv_mmx.cpp

yuv_nosimd: 
	$(CC) -o yuv $(OPTIM) $(NOSIMD) yuv.cpp

yuv_sse: 
	$(CC) -o yuv_sse $(OPTIM) $(SIMD) yuv_sse.cpp

yuv_avx:
	$(CC) -o yuv_avx $(OPTIM) $(SIMD) yuv_avx.cpp

clean:
	-rm $(PROC) 
