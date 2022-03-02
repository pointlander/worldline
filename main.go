// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"

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
	N = 32 * 1024
	// Loops is the number of loops
	Loops = 100
)

var (
	// CPUs is the number of CPUs
	CPUs = runtime.NumCPU()
)

func MakeLoops() [][][d]float64 {
	loops := make([][][d]float64, Loops)
	for loop := 0; loop < Loops; loop++ {
		rnd := rand.New(rand.NewSource(int64(loop + 1)))
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
				point[j] = real(yt[j][i]) / norm[0]
			}
			loops[loop] = append(loops[loop], point)
		}
	}

	return loops
}

func W(a float64, loop [][d]float64, x float64) (float64, float64) {
	intersections := 0.0
	scale := a
	a /= 2
	var previous [d]float64
	for i := 0; i < d; i++ {
		previous[i] = scale*loop[0][0] + x
	}
	length := 0.0
	for i := 1; i < N; i++ {
		var xx [2]float64
		for j := 0; j < d; j++ {
			xx[j] = scale*loop[i][j] + x
		}
		for j := 0; j < d; j++ {
			diff := xx[j] - previous[j]
			length += diff * diff
		}
		if previous[0] < -a && xx[0] > a ||
			previous[0] > a && xx[0] < -a {
			intersections += 2
		} else if previous[0] < a && xx[0] > a ||
			previous[0] > a && xx[0] < a ||
			previous[0] < -a && xx[0] > -a ||
			previous[0] > -a && xx[0] < -a {
			intersections++
		}
		previous = xx
	}

	return intersections, length / 4
}

func main() {
	fmt.Println("making loops...")
	loops := MakeLoops()
	fmt.Println("simulating...")
	type Result struct {
		Intersections float64
		Length        float64
	}
	points := make(plotter.XYs, 0, 10)
	for a := 1.0; a >= .1; a -= .1 {
		w := func(x float64) float64 {
			done := make(chan Result, 8)
			process := func(a float64, i int) {
				intersections, length := W(a, loops[i], x)
				done <- Result{
					Intersections: intersections,
					Length:        length,
				}
			}
			sum, denominator, flight, i := 0.0, 0.0, 0, 0
			for i < Loops && flight < CPUs {
				go process(a, i)
				flight++
				i++
			}
			for i < Loops {
				result := <-done
				flight--
				sum += math.Exp(-result.Intersections) * math.Exp(-result.Length)
				denominator += math.Exp(-result.Length)

				go process(a, i)
				flight++
				i++
			}
			for i := 0; i < flight; i++ {
				result := <-done
				sum += math.Exp(-result.Intersections) * math.Exp(-result.Length)
				denominator += math.Exp(-result.Length)
			}
			return x * sum / denominator
		}
		e := quad.Fixed(w, -1, 1, 1000, nil, 0)
		fmt.Println(a, e)
		points = append(points, plotter.XY{X: a, Y: e})
	}
	for a := .09; a >= .01; a -= .01 {
		w := func(x float64) float64 {
			done := make(chan Result, 8)
			process := func(a float64, i int) {
				intersections, length := W(a, loops[i], x)
				done <- Result{
					Intersections: intersections,
					Length:        length,
				}
			}
			sum, denominator, flight, i := 0.0, 0.0, 0, 0
			for i < Loops && flight < CPUs {
				go process(a, i)
				flight++
				i++
			}
			for i < Loops {
				result := <-done
				flight--
				sum += math.Exp(-result.Intersections) * math.Exp(-result.Length)
				denominator += math.Exp(-result.Length)

				go process(a, i)
				flight++
				i++
			}
			for i := 0; i < flight; i++ {
				result := <-done
				sum += math.Exp(-result.Intersections) * math.Exp(-result.Length)
				denominator += math.Exp(-result.Length)
			}
			return x * sum / denominator
		}
		e := quad.Fixed(w, -1, 1, 1000, nil, 0)
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
