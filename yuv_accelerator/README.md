# SIMD Image Accelerator

An accelerator for YUV image processing using SIMD instructions.

The code mainly does two things: fade in & fade out of a single YUV image, a gradient change of image overlay of two YUV images.

For the first task, the code reads a YUV image and turns it into a RGB format image:

```mathematica
R = Y + 1.140 * V
G = Y - 0.394 * U - 0.581 * V
B = Y + 2.032 * U
```

Then, add an additional brightness dimension, alpha, to the RGB image, making it in ARGB format:

```mathematica
R' = A * R / 256
G' = A * G / 256
B' = A * B / 256
```

The ARGB frames are created in stride 3. That is, 84 distinct frames will be created, each having an alpha value 1, 4, 7, ..., 253.

Lastly, turn the ARGB frames back into one YUV file:

```mathematica
Y = 0.299 * R + 0.587 * G + 0.114 * B
U = -0.147 * R - 0.289 * G + 0.436 * B = 0.492 * (B - Y)
V = 0.615 * R - 0.515 * G - 0.100 * B = 0.877 * (R - Y)
```

For the second task, the only difference is the intermediate step, *i.e.*, we use the idea of alpha blending to perform image overlay:

```mathematica
R' = (A * R1 + (256 - A) * R2) / 256
G' = (A * G1 + (256 - A) * G2) / 256
B' = (A * B1 + (256 - A) * B2) / 256
```

In this project, I first implement the above tasks in raw C++ code. Then, I re-write part of the code using SIMD instructions, namely MMX, SSE and AVX instruction sets. These instruction sets are designed mainly for accelerating vector-level operation, *i.e.*, data are processed in vectors instead of individual integers or floats.

For comparison, the code will compute the time elapse of different implementations. Compiled under g++'s -O2 option, the SIMD instructions largely improves the performance of the code.

**P.S.** The MMX instructions are included in ``<mmintrin.h>``, while SSE and AVX instructions are included in ``<immintrin.h>``. These instruction sets require the ``-msse2 -mmmx -mavx2`` options of g++.

## Requirements

OS: Ubuntu. The code was tested on Ubuntu 16.04 LTS.

A g++ compiler. The code was tested on g++ 5.4.0.

A YUV player. The code was tested using the latest version (7:2.8.15) of ffmpeg, which can be easily installed in Ubuntu:

```bash
sudo apt-get install ffmpeg
```

## Files

``demo/``: directory that contains two source YUV images (named as ``dem1.yuv`` and ``dem2.yuv``). The generated YUV images are also expected to store in this directory.

``yuv.cpp``: code in raw C++ without any SIMD instructions.

``yuv``: executable file of ``yuv.cpp``, compiled by g++.

``yuv_mmx.cpp``: code using MMX instructions to accelerate image processing. Specifically, alpha blending and image overlay are re-written using MMX. Since MMX are not compatible with float operations, rgb2yuv and yuv2rgb functions remain the same.

``yuv_mmx``: executable file of ``yuv_mmx.cpp``, compiled by g++.

``yuv_sse.cpp``: code using SSE instructions to accelerate image processing. All core functions are re-written.

``yuv_sse``: executable file of ``yuv_sse.cpp``, compiled by g++.

``yuv_avx.cpp``: code using AVX instructions to accelerate image processing. All core functions are re-written.

``yuv_avx``: executable file of ``yuv_avx.cpp``, compiled by g++.

## Usage

1. Compile C++ code by typing ``make`` in your terminal. (Optional, you can directly use the provided executable files)

2. Run the executable file. (Here I use the AVX code for example)

   ```bash
   ./yuv_avx
   ```

   The program will print how much time it takes to finish the two tasks, respectively. Then you should find two YUV files generated in ``demo/`` directory. They should have names ``avx_fade.yuv`` and ``avx_overlay.yuv``.

3. Use ``ffmpeg`` to display the YUV image.

   ```bash
   ffplay -f rawvideo -video_size 1920x1080 demo/avx_fade.yuv
   ffplay -f rawvideo -video_size 1920x1080 demo/avx_overlay.yuv
   ```

   
