#include <iostream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <mmintrin.h>
using namespace std;

const int IMG_WIDTH = 1920;
const int IMG_HEIGHT = 1080;
const int IMG_SIZE = 1920 * 1080;
const char* YUV_FILE1 = "demo/dem1.yuv";
const char* YUV_FILE2 = "demo/dem2.yuv";
const char* OUT_FILE1 = "demo/mmx_fade.yuv";
const char* OUT_FILE2 = "demo/mmx_overlay.yuv";

unsigned char yuv_y[2][IMG_SIZE];
unsigned char yuv_u[2][IMG_SIZE / 4];
unsigned char yuv_v[2][IMG_SIZE / 4];
short new_yuv_y[IMG_SIZE];
short new_yuv_u[IMG_SIZE / 4];
short new_yuv_v[IMG_SIZE / 4];
short rgb_r[2][IMG_SIZE];
short rgb_g[2][IMG_SIZE];
short rgb_b[2][IMG_SIZE];
short new_rgb_r[IMG_SIZE];
short new_rgb_g[IMG_SIZE];
short new_rgb_b[IMG_SIZE];
unsigned char raw_yuv_y[IMG_SIZE];
unsigned char raw_yuv_u[IMG_SIZE / 4];
unsigned char raw_yuv_v[IMG_SIZE / 4];

__m64 *dstR, *dstG, *dstB;
__m64 vec_A, vec_R, vec_G, vec_B, vec_R1, vec_R2, vec_G1, vec_G2, vec_B1, vec_B2;
short *srcR, *srcG, *srcB;


void read_yuv(const char* filename, int file_id)
{
	FILE* file = fopen(filename, "rb");
	int index = file_id - 1;  // index in the array  = 0 / 1
	fread(yuv_y[index], 1, IMG_SIZE, file);
	fread(yuv_u[index], 1, IMG_SIZE / 4, file);
	fread(yuv_v[index], 1, IMG_SIZE / 4, file);
	fclose(file);
}


void write_yuv(FILE* file)
{
	for(int i = 0; i < IMG_SIZE; i++)
		raw_yuv_y[i] = new_yuv_y[i];

	for(int i = 0; i < IMG_SIZE / 4; i++) {
		raw_yuv_u[i] = new_yuv_u[i];
		raw_yuv_v[i] = new_yuv_v[i];
	}

	fwrite(raw_yuv_y, 1, IMG_SIZE, file);
	fwrite(raw_yuv_u, 1, IMG_SIZE / 4, file);
	fwrite(raw_yuv_v, 1, IMG_SIZE / 4, file);

	return;
}

void alpha_blending_mmx(int alpha)
{
	srcR = new_rgb_r;
	srcG = new_rgb_g;
	srcB = new_rgb_b;

	vec_A = _mm_set1_pi16((short)alpha);  	// Set 64-bit vec_A with 16-bit integers
	vec_A = _m_psllwi(vec_A, 4);			// Shift packed 16-bit integers in vec_A left by 4 (multiply by 2^4)

	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j += 4) {
			dstR = (__m64*)srcR;
			dstG = (__m64*)srcG;
			dstB = (__m64*)srcB;

			int row_y = i * IMG_HEIGHT + j;
			// R' = A * R / 256
			// Set packed 16-bit integers in 64-bit vec_R with the supplied values.
			// MMX can handle 4 16-bit integers at a time
			vec_R = _mm_set_pi16(rgb_r[0][row_y + 3], rgb_r[0][row_y + 2], rgb_r[0][row_y + 1], rgb_r[0][row_y]);
			vec_R = _m_psllwi(vec_R, 4);		// Shift packed 16-bit integers in vec_R left by 4 (multiply by 2^4)
			*dstR = _m_pmulhw(vec_R, vec_A);  	// Multiply the packed 16-bit integers in vec_R and vec_A, producing intermediate 32-bit int, 
												// and store the high 16 bits of the intermediate integers in dstR. (divided by 2^16)
												// 2^4 * 2^4 / 2^16 = 1 / 2^8 = 1/256
			// G' = A * G / 256
			vec_G = _mm_set_pi16(rgb_g[0][row_y + 3], rgb_g[0][row_y + 2], rgb_g[0][row_y + 1], rgb_g[0][row_y]);
			vec_G = _m_psllwi(vec_G, 4);
			*dstG = _m_pmulhw(vec_G, vec_A);
			// B' = A * B / 256
			vec_B = _mm_set_pi16(rgb_b[0][row_y + 3], rgb_b[0][row_y + 2], rgb_b[0][row_y + 1], rgb_b[0][row_y]);
			vec_B = _m_psllwi(vec_B, 4);
			*dstB = _m_pmulhw(vec_B, vec_A);

			srcR += 4;
			srcG += 4;
			srcB += 4;
		}
}

void image_overlay_mmx(int alpha)
{
	srcR = new_rgb_r;
	srcG = new_rgb_g;
	srcB = new_rgb_b;

	// same operation to alpha as above
	vec_A = _mm_set1_pi16((short)alpha);
	vec_A = _m_psllwi(vec_A, 4);

	for(int i = 0;i < IMG_WIDTH;i++)
		for(int j = 0;j < IMG_HEIGHT;j += 4) {
			dstR = (__m64*)srcR;
			dstG = (__m64*)srcG;
			dstB = (__m64*)srcB;

			int row_y = i * IMG_HEIGHT + j;

			// R' = (A * R1 + (256 - A) * R2) / 256 = A * (R1 - R2) / 256 + R2
			vec_R1 = _mm_set_pi16(rgb_r[0][row_y + 3], rgb_r[0][row_y + 2], rgb_r[0][row_y + 1], rgb_r[0][row_y]); 	// load image 1
			vec_R2 = _mm_set_pi16(rgb_r[1][row_y + 3], rgb_r[1][row_y + 2], rgb_r[1][row_y + 1], rgb_r[1][row_y]);	// load image 2
			vec_R1 = _m_psubsw(vec_R1, vec_R2);		// subtract vec_R2 from vec_R1
			vec_R1 = _m_psllwi(vec_R1, 4);			// shift vec_R1 left by 4
			vec_R1 = _m_pmulhw(vec_R1, vec_A);		// multiply and take higher 16-bit result
			*dstR = _m_paddsw(vec_R1, vec_R2);		// add vec_R1 and vec_R2

			// G' = (A * G1 + (256 - A) * G2) / 256 = A * (G1 - G2) / 256 + G2
			vec_G1 = _mm_set_pi16(rgb_g[0][row_y + 3], rgb_g[0][row_y + 2], rgb_g[0][row_y + 1], rgb_g[0][row_y]);
			vec_G2 = _mm_set_pi16(rgb_g[1][row_y + 3], rgb_g[1][row_y + 2], rgb_g[1][row_y + 1], rgb_g[1][row_y]);
			vec_G1 = _m_psubsw(vec_G1, vec_G2);
			vec_G1 = _m_psllwi(vec_G1, 4);
			vec_G1 = _m_pmulhw(vec_G1, vec_A);
			*dstG = _m_paddsw(vec_G1, vec_G2);

			// B' = (A * B1 + (256 - A) * B2) / 256 = A * (B1 - B2) / 256 + B2
			vec_B1 = _mm_set_pi16(rgb_b[0][row_y + 3], rgb_b[0][row_y + 2], rgb_b[0][row_y + 1], rgb_b[0][row_y]);
			vec_B2 = _mm_set_pi16(rgb_b[1][row_y + 3], rgb_b[1][row_y + 2], rgb_b[1][row_y + 1], rgb_b[1][row_y]);
			vec_B1 = _m_psubsw(vec_B1, vec_B2);
			vec_B1 = _m_psllwi(vec_B1, 4);
			vec_B1 = _m_pmulhw(vec_B1, vec_A);
			*dstB = _m_paddsw(vec_B1, vec_B2);
			
			srcR += 4;
			srcG += 4;
			srcB += 4;
		}
}


void yuv2rgb(int index)
{
	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j++) {

			int row_y = i * IMG_HEIGHT + j;
			int row_uv = (i / 2) * (IMG_HEIGHT / 2) + j / 2;
			// R = Y + 1.140 * V
			rgb_r[index][row_y] = yuv_y[index][row_y] + 1.140 * (yuv_v[index][row_uv] - 128);
			// G = Y - 0.394 * U - 0.581 * V
			rgb_g[index][row_y] = yuv_y[index][row_y] - 0.394 * (yuv_u[index][row_uv] - 128) 
									- 0.581 * (yuv_v[index][row_uv] - 128);
			// B = Y + 2.032 * U
			rgb_b[index][row_y] = yuv_y[index][row_y] + 2.032 * (yuv_u[index][row_uv] - 128);
		}
}


void rgb2yuv()
{
	for(int i = 0; i < IMG_WIDTH; i++)
		for(int j = 0; j < IMG_HEIGHT; j++) {

			int row_y = i * IMG_HEIGHT + j;
			// Y = 0.299 * R + 0.587 * G + 0.114 * B
			new_yuv_y[row_y] = 0.299 * new_rgb_r[row_y] + 0.587 * new_rgb_g[row_y] + 0.114 * new_rgb_b[row_y];

			if(i % 2 == 1 && j % 2 == 1) {
				int row_uv = (IMG_HEIGHT / 2) * (i / 2) + (j / 2);
				// U = -0.147 * R - 0.289 * G + 0.436 * B
				new_yuv_u[row_uv] = - 0.147 * new_rgb_r[row_y] - 0.289 * new_rgb_g[row_y] + 0.436 * new_rgb_b[row_y] + 128;
				// V = 0.615 * R - 0.515 * G - 0.100 * B
				new_yuv_v[row_uv] = 0.615 * new_rgb_r[row_y] - 0.515 * new_rgb_g[row_y] - 0.100 * new_rgb_b[row_y] + 128;
			}
		}
}


void process_picture(int task)
{
	clock_t start_clock = clock();
	FILE* fout = NULL;

	if (task == 1) {
		fout = fopen(OUT_FILE1, "wb");
		yuv2rgb(0);
		for(int alpha = 1; alpha < 255; alpha += 3) {
			alpha_blending_mmx(alpha);
			rgb2yuv();
			write_yuv(fout);
		}
	}
	else if (task == 2) {
		fout = fopen(OUT_FILE2, "wb");
		yuv2rgb(0);
		yuv2rgb(1);
		for(int alpha = 1; alpha < 255; alpha += 3) {
			image_overlay_mmx(alpha);
			rgb2yuv();
			write_yuv(fout);
		}
	}

	fclose(fout);
	clock_t total_time = clock() - start_clock;

	if(task == 1)
		cout << "Task1: Image Fading" << endl;
	else if (task == 2)
		cout << "Task2: Image Overlay" << endl;

	cout << "Time Elapse: " << fixed << setprecision(3) << (double)total_time / (double)CLOCKS_PER_SEC << "s" << endl;
}


int main()
{	
	read_yuv(YUV_FILE1, 1);
	read_yuv(YUV_FILE2, 2);
	process_picture(1);
	process_picture(2);
	return 0;
}
