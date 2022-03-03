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
	// D is the number of space time dimensions
	D = d + 1
	// N is the number of points in the Worldline
	N = 32 * 1024
	// Loops is the number of loops
	Loops = 100
	// Lambda is a plate factor
	Lambda = 1e-4
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
				point[j] = real(yt[j][i]) / norm[j]
			}
			loops[loop] = append(loops[loop], point)
		}
	}

	return loops
}

func W(a float64, loop [][d]float64, x float64) (float64, float64) {
	intersections := 0.0
	a /= 2
	length := 0.0
	for i := 0; i < N+1; i++ {
		v1, v2 := loop[(i+N-1)%N], loop[i%N]
		for j := 0; j < d; j++ {
			diff := v1[j] - v2[j]
			length += diff * diff
		}
		if x1, x2 := v1[0]+x, v2[0]+x; x1 < -a && x2 > a ||
			x1 > a && x2 < -a {
			intersections += 2
		} else if x1 < a && x2 > a ||
			x1 > a && x2 < a ||
			x1 < -a && x2 > -a ||
			x1 > -a && x2 < -a {
			intersections++
		}
	}

	return Lambda * intersections, length / 4
}

func main() {
	fmt.Println("making loops...")
	loops := MakeLoops()

	fmt.Println("simulating...")
	factor := -1 / (2 * math.Pow(4*math.Pi, D/2))
	m := 1.0
	f := func(T float64) float64 {
		return math.Exp(-math.Pow(m, 2)*T) / math.Pow(T, 1+D/2)
	}
	factor *= quad.Fixed(f, 0, math.Inf(1), 1000, nil, 0)

	type Result struct {
		Intersections float64
		Length        float64
	}
	points := make(plotter.XYs, 0, 10)
	for a := 1.0; a > 0; a -= .01 {
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
			return x * (sum/denominator - 1)
		}
		e := factor * quad.Fixed(w, -1, 1, 1000, nil, 0)
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
