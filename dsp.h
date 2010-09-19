//      dsp.h
//      
//      Copyright 2010 Ilya <ilya@laptop>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#ifndef _DSP_H_
#define _DSP_H_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#include "types.h"

void process_fftw_out_into_harmonics(fftw_complex * out, struct Harmonic * harmonics, int sample_size, int sample_msec);
int find_main_harmonic(struct Harmonic * harmonics, struct Harmonic * dest, int sample_size);
int lookup_for_harmonics(struct Harmonic * harmonics, struct Harmonic * spectre, int current, int sample_size);
void generate_chord(struct Harmonic * spectre, double * out_buffer, struct Chord chord, int sample_rate, int sample_size, int pos, int offset);
void normalize_out_buffer(double * out_buffer, char * sample_buffer, int sample_size);
double compute_freq(int note, int octave);
void push_raw_le_data_into_fftw(char * source, fftw_complex * in, int sample_size);


#endif /* _DSP_H_ */
