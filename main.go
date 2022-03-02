// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/mjibson/go-dsp/fft"
	"gonum.org/v1/gonum/integrate/quad"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

const (
	// d is the number of dimensions
	d = 2
	// N is the number of points in the Worldline
	N = 1024
)

func W(a float64, worldline int64, x float64) float64 {
	rnd := rand.New(rand.NewSource(worldline))
	y := make([][]complex128, 0, d)
	for i := 0; i < d; i++ {
		y = append(y, make([]complex128, 0, N))
		y[i] = append(y[i], 0)
	}
	factor := math.Sqrt(2)
	for i := 0; i < d; i++ {
		for j := 1; j < N; j++ {
			y[i] = append(y[i], complex(rnd.Float64()*factor, rnd.Float64()*factor))
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
	for i := 1; i < N; i++ {
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

	intersections := 0.0
	previous := real(yt[0][0])/norm[0] + x
	a /= 2
	for i := 1; i < N; i++ {
		x := real(yt[0][i])/norm[0] + x
		if previous < -a && x > a ||
			previous > a && x < -a {
			intersections += 2
		} else if previous < a && x > a ||
			previous > a && x < a ||
			previous < -a && x > -a ||
			previous > -a && x < -a {
			intersections++
		}
		previous = x
	}

	return intersections
}

func main() {
	points := make(plotter.XYs, 0, 10)
	for a := 1.5; a >= .1; a -= .1 {
		w := func(x float64) float64 {
			sum := 0.0
			for i := 0; i < 10; i++ {
				sum += math.Exp(-W(a, int64(i), x))
			}
			return x * sum
		}
		e := quad.Fixed(w, -2, 2, 10000, nil, 0)
		fmt.Println(a, e)
		points = append(points, plotter.XY{X: a, Y: e})
	}
	for a := .09; a >= .01; a -= .01 {
		w := func(x float64) float64 {
			sum := 0.0
			for i := 0; i < 10; i++ {
				sum += math.Exp(-W(a, int64(i), x))
			}
			return x * sum
		}
		e := quad.Fixed(w, -2, 2, 10000, nil, 0)
		fmt.Println(a, e)
		points = append(points, plotter.XY{X: a, Y: e})
	}

	p := plot.New()

	p.Title.Text = "a vs energy"
	p.X.Label.Text = "a"
	p.Y.Label.Text = "energy"

	scatter, err := plotter.NewScatter(points)
	if err != nil {
		panic(err)
	}
	scatter.GlyphStyle.Radius = vg.Length(1)
	scatter.GlyphStyle.Shape = draw.CircleGlyph{}
	p.Add(scatter)

	err = p.Save(8*vg.Inch, 8*vg.Inch, "energy.png")
	if err != nil {
		panic(err)
	}
}
