// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/mjibson/go-dsp/fft"
)

const (
	// d is the number of dimensions
	d = 2
	// N is the number of points in the Worldline
	N = 1024
)

func main() {
	rnd := rand.New(rand.NewSource(1))
	x := make([][]complex128, 0, d)
	for i := 0; i < d; i++ {
		x = append(x, make([]complex128, 0, N))
		x[i] = append(x[i], 0)
	}
	factor := math.Sqrt(2)
	for i := 0; i < d; i++ {
		for j := 1; j < N; j++ {
			x[i] = append(x[i], complex(rnd.Float64()*factor, rnd.Float64()*factor))
		}
	}
	xfft := make([][]complex128, 0, d)
	for i := 0; i < d; i++ {
		xfft = append(xfft, fft.FFT(x[i]))
	}

	sum := make([]float64, 2)
	manhattan := make([]float64, 2)
	length := 0.0
	first := []float64{real(xfft[0][0]), real(xfft[1][0])}
	last := []float64{first[0], first[1]}
	for i := 1; i < N; i++ {
		var diff [2]float64
		for j := 0; j < d; j++ {
			diff[j] = math.Abs(last[j] - real(xfft[j][i]))
			manhattan[j] += diff[j]
			sum[j] += real(xfft[j][i])
			last[j] = real(xfft[j][i])
		}
		length += math.Sqrt(diff[0]*diff[0] + diff[1]*diff[1])
		fmt.Println(real(xfft[0][i]), real(xfft[1][i]))
	}
	fmt.Printf("\n")
	var diff [2]float64
	for j := 0; j < d; j++ {
		diff[j] = math.Abs(last[j] - first[j])
		manhattan[j] += diff[j]
	}
	length += math.Sqrt(diff[0]*diff[0] + diff[1]*diff[1])
	fmt.Println(sum[0], sum[1])
	fmt.Println(manhattan[0], manhattan[1], length)
}
