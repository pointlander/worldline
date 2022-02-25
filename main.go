// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"

	"github.com/mjibson/go-dsp/fft"
)

func main() {
	rnd := rand.New(rand.NewSource(1))
	x, y := make([]complex128, 0, 1024), make([]complex128, 0, 1024)
	x = append(x, 0)
	y = append(y, 0)
	factor := math.Sqrt(2)
	for i := 1; i < 1024; i++ {
		x = append(x, complex(rnd.Float64()*factor, rnd.Float64()*factor))
		y = append(y, complex(rnd.Float64()*factor, rnd.Float64()*factor))
	}
	xx := fft.FFT(x)
	yy := fft.FFT(y)

	sumxx, sumyy := 0.0, 0.0
	lengthxx, lengthyy := 0.0, 0.0
	length := 0.0
	firstxx, firstyy := cmplx.Abs(xx[0]), cmplx.Abs(yy[0])
	lastxx, lastyy := firstxx, firstyy
	for i := 1; i < 1024; i++ {
		x, y := cmplx.Abs(xx[i]), cmplx.Abs(yy[i])
		diffxx := math.Abs(lastxx-x) / 1024
		diffyy := math.Abs(lastyy-y) / 1024
		lengthxx += diffxx
		lengthyy += diffyy
		length += math.Sqrt(diffxx*diffxx + diffyy*diffyy)
		sumxx += x
		sumyy += y
		fmt.Println(x, y)
		lastxx, lastyy = x, y
	}
	fmt.Printf("\n")
	diffxx := math.Abs(lastxx-firstxx) / 1024
	diffyy := math.Abs(lastyy-firstyy) / 1024
	lengthxx += diffxx
	lengthyy += diffyy
	length += math.Sqrt(diffxx*diffxx + diffyy*diffyy)
	fmt.Println(sumxx, sumyy)
	fmt.Println(lengthxx, lengthyy, length)
}
