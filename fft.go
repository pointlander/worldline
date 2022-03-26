// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"
	"math/rand"

	"github.com/mjibson/go-dsp/fft"
)

// MakeLoops make worldline loops
func MakeLoopsFFT(N, Loops int) []Worldline {
	loops := make([]Worldline, Loops)
	for loop := 0; loop < Loops; loop++ {
		rnd := rand.New(rand.NewSource(int64(loop + 1)))
		y := make([][]complex128, 0, d)
		for i := 0; i < d; i++ {
			y = append(y, make([]complex128, 0, N))
			y[i] = append(y[i], 0)
		}
		for i := 0; i < d; i++ {
			for j := 1; j < N; j++ {
				y[i] = append(y[i], complex(rnd.NormFloat64(), rnd.NormFloat64()))
			}
		}

		yt := make([][]complex128, 0, d)
		for i := 0; i < d; i++ {
			yt = append(yt, fft.FFT(y[i]))
		}

		var min [d]float64
		var max [d]float64
		for i := 0; i < d; i++ {
			min[i] = math.MaxFloat64
			max[i] = -math.MaxFloat64
		}
		for i := 0; i < N; i++ {
			for j := 0; j < d; j++ {
				r := real(yt[j][i])
				if r < min[j] {
					min[j] = r
				}
				if r > max[j] {
					max[j] = r
				}
			}
		}

		var norm [d]float64
		for i := 0; i < d; i++ {
			norm[i] = math.Abs(max[i] - min[i])
		}

		for i := 0; i < N; i++ {
			var point [d]float64
			for j := 0; j < d; j++ {
				point[j] = real(yt[j][i]) / norm[j]
			}
			loops[loop].Line = append(loops[loop].Line, point)
		}
	}

	return loops
}
