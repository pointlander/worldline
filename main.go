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

	sum := make([]float64, d)
	manhattan := make([]float64, d)
	length := 0.0
	first := make([]float64, d)
	last := make([]float64, d)
	for i := 0; i < d; i++ {
		first[i] = real(xfft[i][0])
		last[i] = real(xfft[i][0])
	}
	var min [d]float64
	var max [d]float64
	for i := 0; i < d; i++ {
		min[i] = math.MaxFloat64
		max[i] = -math.MaxFloat64
	}
	for i := 1; i < N; i++ {
		var diff [d]float64
		for j := 0; j < d; j++ {
			r := real(xfft[j][i])
			diff[j] = math.Abs(last[j]) - r
			manhattan[j] += diff[j]
			sum[j] += r
			last[j] = r
			if r < min[j] {
				min[j] = r
			}
			if r > max[j] {
				max[j] = r
			}
		}
		l := 0.0
		for j := 0; j < d; j++ {
			l += diff[j] * diff[j]
		}
		length += math.Sqrt(l)
		fmt.Println(real(xfft[0][i]), real(xfft[1][i]))
	}
	fmt.Printf("\n")
	var diff [d]float64
	for j := 0; j < d; j++ {
		diff[j] = math.Abs(last[j] - first[j])
		manhattan[j] += diff[j]
	}
	l := 0.0
	for j := 0; j < d; j++ {
		l += diff[j] * diff[j]
	}
	length += math.Sqrt(l)
	fmt.Println(sum[0], sum[1])
	fmt.Println(manhattan[0], manhattan[1], length)

	var norm [d]float64
	for i := 0; i < d; i++ {
		norm[i] = math.Abs(max[i]-min[i]) / 2
	}
	for i := 0; i < N; i++ {
		fmt.Println(real(xfft[0][i])/norm[0], real(xfft[1][i])/norm[1])
	}
}
