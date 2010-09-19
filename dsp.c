//      dsp.c
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


#include "dsp.h"

double compute_freq(int note, int octave) {
	return 27.5 * pow(2, (octave * 12 + note) / 12.0L);
}

void push_raw_le_data_into_fftw(char * source, fftw_complex * in, int sample_size) {
	int j;
	char tmp_buf[2];

	for (j = 0; j < sample_size; ++j) {
		memcpy(&tmp_buf, (char *)(source + 4 * j), 2 * sizeof(char));
		memcpy(&tmp_buf, (char *)(source + 4 * j + 2), 2 * sizeof(char));
		in[j][0] = 0xff * tmp_buf[1] + tmp_buf[0] * ((tmp_buf[1] < 0) ? -1 : 1); // Little endian :)
		in[j][1] = 0;	
	}
}

void process_fftw_out_into_harmonics(fftw_complex * out, struct Harmonic * harmonics, int sample_size, int sample_msec) {
	int j;
	
	for (j = 0; j < sample_size / 2; ++j) {
		harmonics[j].ampl = sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]) / sample_size;
		harmonics[j].phase = atan2(out[j][1], out[j][0]);
		harmonics[j].freq = j / ((double)sample_msec / 1000);
	}
}


int find_main_harmonic(struct Harmonic * harmonics, struct Harmonic * dest, int sample_size) {
	int j, current = 0;
	
	memcpy(dest, &harmonics[0], sizeof(struct Harmonic));
	for (j = 1; j < sample_size / 2; ++j)
		if ( (harmonics[j].ampl > dest->ampl) && (current = j) ) 
			memcpy(dest, &harmonics[j], sizeof(struct Harmonic));
	
	return current;
}


int lookup_for_harmonics(struct Harmonic * harmonics, struct Harmonic * spectre, int current, int sample_size) {
	int j, ch, pos = 1;	
	
	current *= 2;
	while (pos < 10) {
		memcpy(&spectre[pos], &harmonics[current - 10], sizeof(struct Harmonic));
		ch = current - 10;
		// Looking window contains 20 frequencies around n * f (needed for accurate syntheses)
		for (j = current - 9; j < current + 10; ++j)
			if ( (harmonics[j].ampl > spectre[pos].ampl) && (ch = j) ) 
				memcpy(&spectre[pos], &harmonics[j], sizeof(struct Harmonic));
		current = ch * 2;
		// Skip harmonics beyond FFT's resolution
		if (current + 10 > sample_size / 2)
			break;
		++pos;
	}
	return pos;
}

void generate_chord(struct Harmonic * spectre, double * out_buffer, struct Chord chord, int sample_rate, int sample_size, int pos, int offset) {
	int samples[5];
	int j, k, mp, phase;
	double phases[ 3 * 10 ];
	double amplitude, sample;
	const double dpi = 6.2831852;
	
	for (j = 0; j < pos; ++j) {
		mp = j * 3;
		phases[mp] = dpi / (float)sample_rate * compute_freq(chord.note[0].note, chord.note[0].octave) * (j + 1);
		phases[mp + 1] = dpi / (float)sample_rate * compute_freq(chord.note[1].note, chord.note[1].octave) * (j + 1);
		phases[mp + 2] = dpi / (float)sample_rate * compute_freq(chord.note[2].note, chord.note[2].octave) * (j + 1);
	}
	// Generate and upmix chord directly into out buffer
	for (j = 0; j < sample_size; ++j) {
		sample = 0;
		phase = offset + j;
		for (k = 0; k < pos; ++k) {
			mp = k * 3;
			sample += (spectre[k].ampl / 3) * (sin(phase * phases[mp]) + sin(phase * phases[mp + 1]) + sin(phase * phases[mp + 2]));	
		}
		out_buffer[2 * j] = out_buffer[2 * j + 1] = sample;
	}
}

void normalize_out_buffer(double * out_buffer, char * sample_buffer, int sample_size) {
	int i, sample;
	double t, max = 0;
	double multiplier;
	
	for (i = 0; i < sample_size; ++i) {
		t = fabs(out_buffer[2 * i]);
		if (t > max)
			max = t;
	}
	multiplier = (0xffff / 2 - 0xff) / max;

	for (i = 0; i < sample_size; ++i) {
		out_buffer[2 * i] *= multiplier;
		out_buffer[2 * i + 1] *= multiplier;
		sample = (int)out_buffer[2 * i];
		sample_buffer[4 * i] = sample_buffer[4 * i + 2] = sample & 0xff;
		sample_buffer[4 * i + 1] = sample_buffer[4 * i + 3] = (sample >> 8) & 0xff;		
	}
}

