// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math/cmplx"
	"math/rand"

	"github.com/mjibson/go-dsp/fft"
)

func main() {
	rnd := rand.New(rand.NewSource(1))
	x, y := make([]complex128, 0, 1024), make([]complex128, 0, 1024)
	for i := 0; i < 1024; i++ {
		x = append(x, complex(rnd.Float64()*2-1, rnd.Float64()*2-1))
		y = append(y, complex(rnd.Float64()*2-1, rnd.Float64()*2-1))
	}
	xx := fft.FFT(x)
	yy := fft.FFT(y)

	sumxx, sumyy := 0.0, 0.0
	for i := 0; i < 1024; i++ {
		sumxx += cmplx.Abs(xx[i])
		sumyy += cmplx.Abs(yy[i])
		fmt.Println(cmplx.Abs(xx[i]), cmplx.Abs(yy[i]))
	}
	fmt.Println(sumxx, sumyy)
}
