#include <iostream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <immintrin.h>
using namespace std;

const int IMG_WIDTH = 1920;
const int IMG_HEIGHT = 1080;
const int IMG_SIZE = 1920 * 1080;
const char* YUV_FILE1 = "demo/dem1.yuv";
const char* YUV_FILE2 = "demo/dem2.yuv";
const char* OUT_FILE1 = "demo/sse_fade.yuv";
const char* OUT_FILE2 = "demo/sse_overlay.yuv";

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

__m128i *dstR, *dstG, *dstB, *dstY, *dstU, *dstV;
__m128i vec_Y, vec_U, vec_V, vec_IMM, vec_OFF, vec_A, vec_A1, vec_A2,
		vec_R, vec_G, vec_B, vec_R1, vec_R2, vec_G1, vec_G2, vec_B1, vec_B2, 
		result1, result2, result3;
short *srcY, *srcU, *srcV, *srcR, *srcG, *srcB;


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


void alpha_blending_sse(int alpha)
{
	srcR = new_rgb_r;
	srcG = new_rgb_g;
	srcB = new_rgb_b;

	vec_A = _mm_set1_epi16((short)alpha);	// set 128-bit vector with 16-bit integers

	for(int i = 0; i < IMG_HEIGHT; i++)
		for(int j = 0; j < IMG_WIDTH; j += 8) {
			dstR = (__m128i*)srcR;
			dstG = (__m128i*)srcG;
			dstB = (__m128i*)srcB;

			int row_y = i * IMG_WIDTH + j;

			// R' = A * R / 256
			vec_R = _mm_set_epi16(rgb_r[0][row_y + 7], rgb_r[0][row_y + 6], 
								  rgb_r[0][row_y + 5], rgb_r[0][row_y + 4],
								  rgb_r[0][row_y + 3], rgb_r[0][row_y + 2], 
								  rgb_r[0][row_y + 1], rgb_r[0][row_y]);

			vec_R = _mm_mullo_epi16(vec_R, vec_A);		// multiply vec_R and vec_A, select low 16 bits
			vec_R = _mm_srli_epi16(vec_R, 8);			// vec_R = vec_R / 256
			*dstR = vec_R;

			// G' = A * G / 256
			vec_G = _mm_set_epi16(rgb_g[0][row_y + 7], rgb_g[0][row_y + 6], 
								  rgb_g[0][row_y + 5], rgb_g[0][row_y + 4],
								  rgb_g[0][row_y + 3], rgb_g[0][row_y + 2], 
								  rgb_g[0][row_y + 1], rgb_g[0][row_y]);

			vec_G = _mm_mullo_epi16(vec_G, vec_A);
			*dstG = _mm_srli_epi16(vec_G, 8);

			// B' = A * B / 256
			vec_B = _mm_set_epi16(rgb_b[0][row_y + 7], rgb_b[0][row_y + 6], 
			    				  rgb_b[0][row_y + 5], rgb_b[0][row_y + 4],
								  rgb_b[0][row_y + 3], rgb_b[0][row_y + 2], 
								  rgb_b[0][row_y + 1], rgb_b[0][row_y]);

			vec_B = _mm_mullo_epi16(vec_B, vec_A);
			*dstB = _mm_srli_epi16(vec_B, 8);

			srcR += 8;
			srcG += 8;
			srcB += 8;
		}
}


void image_overlay_sse(int alpha)
{
	srcR = new_rgb_r;
	srcG = new_rgb_g;
	srcB = new_rgb_b;
	
	vec_A1 = _mm_set1_epi16((short)alpha);		// set 128-bit vector with 16-bit integers
	vec_A2 = _mm_set1_epi16((short)256 - alpha);

	for(int i = 0; i < IMG_HEIGHT; i++)
		for(int j = 0; j < IMG_WIDTH; j += 8) {
			dstR = (__m128i*)srcR;
			dstG = (__m128i*)srcG;
			dstB = (__m128i*)srcB;

			int row_y = i * IMG_WIDTH + j;

			// R' = (A * R1 + (256 - A) * R2) / 256
			vec_R1 = _mm_set_epi16(rgb_r[0][row_y + 7], rgb_r[0][row_y + 6], 
								   rgb_r[0][row_y + 5], rgb_r[0][row_y + 4],
								   rgb_r[0][row_y + 3], rgb_r[0][row_y + 2], 
								   rgb_r[0][row_y + 1], rgb_r[0][row_y]);

			vec_R2 = _mm_set_epi16(rgb_r[1][row_y + 7], rgb_r[1][row_y + 6], 
								   rgb_r[1][row_y + 5], rgb_r[1][row_y + 4],
								   rgb_r[1][row_y + 3], rgb_r[1][row_y + 2], 
								   rgb_r[1][row_y + 1], rgb_r[1][row_y]);

			vec_R1 = _mm_mullo_epi16(vec_R1, vec_A1);	// A * R1, get lower 16 bits
			vec_R2 = _mm_mullo_epi16(vec_R2, vec_A2);	// (256 - A) * R2, get lower 16 bits
			vec_R = _mm_add_epi16(vec_R1, vec_R2);		// vec_R = A * R1 + (256 - A) * R2
			*dstR = _mm_srli_epi16(vec_R, 8);			// vec_R = vec_R / 256

			// G' = (A * G1 + (256 - A) * G2) / 256
			vec_G1 = _mm_set_epi16(rgb_g[0][row_y + 7], rgb_g[0][row_y + 6], 
								   rgb_g[0][row_y + 5], rgb_g[0][row_y + 4],
								   rgb_g[0][row_y + 3], rgb_g[0][row_y + 2], 
								   rgb_g[0][row_y + 1], rgb_g[0][row_y]);

			vec_G2 = _mm_set_epi16(rgb_g[1][row_y + 7], rgb_g[1][row_y + 6], 
								   rgb_g[1][row_y + 5], rgb_g[1][row_y + 4],
								   rgb_g[1][row_y + 3], rgb_g[1][row_y + 2], 
								   rgb_g[1][row_y + 1], rgb_g[1][row_y]);

			vec_G1 = _mm_mullo_epi16(vec_G1, vec_A1);
			vec_G2 = _mm_mullo_epi16(vec_G2, vec_A2);
			vec_G = _mm_add_epi16(vec_G1, vec_G2);
			*dstG = _mm_srli_epi16(vec_G, 8);

			// B' = (A * B1 + (256 - A) * B2) / 256
			vec_B1 = _mm_set_epi16(rgb_b[0][row_y + 7], rgb_b[0][row_y + 6], 
								   rgb_b[0][row_y + 5], rgb_b[0][row_y + 4],
								   rgb_b[0][row_y + 3], rgb_b[0][row_y + 2], 
								   rgb_b[0][row_y + 1], rgb_b[0][row_y]);

			vec_B2 = _mm_set_epi16(rgb_b[1][row_y + 7], rgb_b[1][row_y + 6], 
								   rgb_b[1][row_y + 5], rgb_b[1][row_y + 4],
								   rgb_b[1][row_y + 3], rgb_b[1][row_y + 2], 
								   rgb_b[1][row_y + 1], rgb_b[1][row_y]);

			vec_B1 = _mm_mullo_epi16(vec_B1, vec_A1);
			vec_B2 = _mm_mullo_epi16(vec_B2, vec_A2);
			vec_B = _mm_add_epi16(vec_B1, vec_B2);
			*dstB = _mm_srli_epi16(vec_B, 8);

			srcR += 8;
			srcG += 8;
			srcB += 8;
		}
}


void yuv2rgb_sse(int index)
{
	srcR = rgb_r[index];
	srcG = rgb_g[index];
	srcB = rgb_b[index];

	for(int i = 0; i < IMG_HEIGHT; i++)
		for(int j = 0; j < IMG_WIDTH; j += 8) {
			dstR = (__m128i*)srcR;
			dstG = (__m128i*)srcG;
			dstB = (__m128i*)srcB;

			int row_y = i * IMG_WIDTH + j;
			int row_uv = (i / 2) * (IMG_WIDTH / 2) + j / 2;

			// set 128-bit vector with 16-bit integers
			vec_Y = _mm_set_epi16(yuv_y[index][row_y + 7], yuv_y[index][row_y + 6], 
								  yuv_y[index][row_y + 5], yuv_y[index][row_y + 4],
								  yuv_y[index][row_y + 3], yuv_y[index][row_y + 2], 
								  yuv_y[index][row_y + 1], yuv_y[index][row_y]);

			vec_U = _mm_set_epi16(yuv_u[index][row_uv + 3], yuv_u[index][row_uv + 3], 
								  yuv_u[index][row_uv + 2], yuv_u[index][row_uv + 2],
								  yuv_u[index][row_uv + 1], yuv_u[index][row_uv + 1], 
								  yuv_u[index][row_uv], yuv_u[index][row_uv]);
								  
			vec_V = _mm_set_epi16(yuv_v[index][row_uv + 3], yuv_v[index][row_uv + 3], 
								  yuv_v[index][row_uv + 2], yuv_v[index][row_uv + 2],
								  yuv_v[index][row_uv + 1], yuv_v[index][row_uv + 1], 
								  yuv_v[index][row_uv], yuv_v[index][row_uv]);

			vec_IMM = _mm_set1_epi16(128);
			vec_U = _mm_sub_epi16(vec_U, vec_IMM);		// U - 128 (two operands must be packed 128-bit vectors)
			vec_V = _mm_sub_epi16(vec_V, vec_IMM);		// V - 128

			// R = Y + 1.140 * V
			vec_IMM = _mm_set1_epi16((1 << 16) * 0.140);
			vec_R = _mm_mulhi_epi16(vec_V, vec_IMM);
			vec_R = _mm_add_epi16(vec_R, vec_V);			// 0.140 * V + V = 1.140 * V
			*dstR = _mm_add_epi16(vec_Y, vec_R);

			// G = Y - 0.394 * U - 0.581 * V
			vec_IMM = _mm_set1_epi16((1 << 16) * 0.394);
			result1 = _mm_mulhi_epi16(vec_U, vec_IMM);
			vec_IMM = _mm_set1_epi16((1 << 16) * 0.081);
			result2 = _mm_mulhi_epi16(vec_V, vec_IMM);
			vec_IMM = _mm_srai_epi16(vec_V, 1);
			result2 = _mm_add_epi16(result2, vec_IMM);		// 0.081 * V + V / 2 = 0.581 * V
			vec_G = _mm_sub_epi16(vec_Y, result1);
			*dstG = _mm_sub_epi16(vec_G, result2);

			// B = Y + 2.032 * U
			vec_IMM = _mm_set1_epi16((1 << 16) * 0.032);
			vec_B = _mm_mulhi_epi16(vec_U, vec_IMM);
			vec_B = _mm_add_epi16(vec_B, vec_U);
			vec_B = _mm_add_epi16(vec_B, vec_U);			// 0.032 * U + U + U = 2.032 * U
			*dstB = _mm_add_epi16(vec_Y, vec_B);

			srcR += 8;
			srcG += 8;
			srcB += 8;
		}
}


void rgb2yuv_sse()
{
	srcY = new_yuv_y;
	srcU = new_yuv_u;
	srcV = new_yuv_v;
	
	for(int i = 0; i < IMG_HEIGHT; i++)
		for(int j = 0; j < IMG_WIDTH; j += 8) {
			dstY = (__m128i*)srcY;

			int row_y = i * IMG_WIDTH + j;
			
			// Y = 0.299 * R + 0.587 * G + 0.114 * B
			vec_R = _mm_set_epi16(new_rgb_r[row_y + 7], new_rgb_r[row_y + 6], 
								  new_rgb_r[row_y + 5], new_rgb_r[row_y + 4],
								  new_rgb_r[row_y + 3], new_rgb_r[row_y + 2], 
								  new_rgb_r[row_y + 1], new_rgb_r[row_y]);
			vec_G = _mm_set_epi16(new_rgb_g[row_y + 7], new_rgb_g[row_y + 6], 
								  new_rgb_g[row_y + 5], new_rgb_g[row_y + 4],
								  new_rgb_g[row_y + 3], new_rgb_g[row_y + 2], 
								  new_rgb_g[row_y + 1], new_rgb_g[row_y]);
			vec_B = _mm_set_epi16(new_rgb_b[row_y + 7], new_rgb_b[row_y + 6], 
								  new_rgb_b[row_y + 5], new_rgb_b[row_y + 4],
								  new_rgb_b[row_y + 3], new_rgb_b[row_y + 2], 
								  new_rgb_b[row_y + 1], new_rgb_b[row_y]);

			vec_IMM = _mm_set1_epi16((1 << 16) * 0.299);
			result1 = _mm_mulhi_epi16(vec_IMM, vec_R);

			vec_IMM = _mm_set1_epi16((1 << 16) * 0.087);
			result2 = _mm_mulhi_epi16(vec_IMM, vec_G);
			vec_IMM = _mm_srai_epi16(vec_G, 1);
			result2 = _mm_add_epi16(result2, vec_IMM);			// 0.087 * G + G / 2 = 0.587 * G

			vec_IMM = _mm_set1_epi16((1 << 16) * 0.114);
			result3 = _mm_mulhi_epi16(vec_IMM, vec_B);

			vec_Y = _mm_add_epi16(result1, result2);
			vec_Y = _mm_add_epi16(vec_Y, result3);
			*dstY = vec_Y;										// add up R channel, G channel and B channel

			if(i % 2 == 0 && j % 16 == 0) {
				dstU = (__m128i*)srcU;
				dstV = (__m128i*)srcV;

				vec_R = _mm_set_epi16(new_rgb_r[row_y + 14], new_rgb_r[row_y + 12], 
									  new_rgb_r[row_y + 10], new_rgb_r[row_y + 8],
									  new_rgb_r[row_y + 6], new_rgb_r[row_y + 4], 
									  new_rgb_r[row_y + 2], new_rgb_r[row_y]);
				vec_G = _mm_set_epi16(new_rgb_g[row_y + 14], new_rgb_g[row_y + 12], 
									  new_rgb_g[row_y + 10], new_rgb_g[row_y + 8],
									  new_rgb_g[row_y + 6], new_rgb_g[row_y + 4], 
									  new_rgb_g[row_y + 2], new_rgb_g[row_y]);
				vec_B = _mm_set_epi16(new_rgb_b[row_y + 14], new_rgb_b[row_y + 12], 
									  new_rgb_b[row_y + 10], new_rgb_b[row_y + 8],
									  new_rgb_b[row_y + 6], new_rgb_b[row_y + 4], 
									  new_rgb_b[row_y + 2], new_rgb_b[row_y]);

				vec_OFF = _mm_set1_epi16(128);

				// U = - 0.147 * R - 0.289 * G + 0.436 * B
				vec_IMM = _mm_set1_epi16((1 << 16) * 0.147);
				result1 = _mm_mulhi_epi16(vec_IMM, vec_R);

				vec_IMM = _mm_set1_epi16((1 << 16) * 0.289);
				result2 = _mm_mulhi_epi16(vec_IMM, vec_G);

				vec_IMM = _mm_set1_epi16((1 << 16) * 0.436);
				result3 = _mm_mulhi_epi16(vec_IMM, vec_B);

				vec_U = _mm_sub_epi16(result3, result2);
				vec_U = _mm_sub_epi16(vec_U, result1);
				*dstU = _mm_add_epi16(vec_U, vec_OFF);

				// V = 0.615 * R - 0.515 * G - 0.100 * B
				vec_IMM = _mm_set1_epi16((1 << 16) * 0.115);			// 0.115 * R + R / 2 = 0.615 * R
				result1 = _mm_mulhi_epi16(vec_IMM, vec_R);
				vec_IMM = _mm_srai_epi16(vec_R, 1);
				result1 = _mm_add_epi16(result1, vec_IMM);

				vec_IMM = _mm_set1_epi16((1 << 16) * 0.015);			// 0.015 * G + G / 2 = 0.515 * G
				result2 = _mm_mulhi_epi16(vec_IMM, vec_G);
				vec_IMM = _mm_srai_epi16(vec_G, 1);
				result2 = _mm_add_epi16(result2, vec_IMM);

				vec_IMM = _mm_set1_epi16((1 << 16) * 0.100);
				result3 = _mm_mulhi_epi16(vec_IMM, vec_B);

				vec_V = _mm_sub_epi16(result1, result2);
				vec_V = _mm_sub_epi16(vec_V, result3);
				*dstV = _mm_add_epi16(vec_V, vec_OFF);

				srcU += 8;
				srcV += 8;
			}
			srcY += 8;
		}
}

void process_picture(int task)
{
	clock_t start_clock = clock();
	FILE* fout = NULL;

	if (task == 1) {
		fout = fopen(OUT_FILE1, "wb");
		yuv2rgb_sse(0);
		for(int alpha = 1; alpha < 255; alpha += 3) {
			alpha_blending_sse(alpha);
			rgb2yuv_sse();
			write_yuv(fout);
		}
	}
	else if (task == 2) {
		fout = fopen(OUT_FILE2, "wb");
		yuv2rgb_sse(0);
		yuv2rgb_sse(1);
		for(int alpha = 1; alpha < 255; alpha += 3) {
			image_overlay_sse(alpha);
			rgb2yuv_sse();
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
