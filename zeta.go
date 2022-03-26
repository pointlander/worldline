// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"

	"gonum.org/v1/gonum/mathext"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

// MakeLoopsZeta makes a loop using the zeta function
func MakeLoopsZeta(N, Loops int) []Worldline {
	loops := make([]Worldline, 1)
	x := 2.0
	for i := 0; i < N; i++ {
		var point [d]float64
		for j := 0; j < d; j++ {
			point[j] = mathext.Zeta(float64(x), float64(x))
		}
		x += .01
		loops[0].Line = append(loops[0].Line, point)
	}

	var min [d]float64
	var max [d]float64
	for i := 0; i < d; i++ {
		min[i] = math.MaxFloat64
		max[i] = -math.MaxFloat64
	}
	for i := 0; i < N; i++ {
		for j := 0; j < d; j++ {
			r := loops[0].Line[i][j]
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

	points := make(plotter.XYs, 0, N)
	for i := 0; i < N; i++ {
		for j := 0; j < d; j++ {
			loops[0].Line[i][j] /= norm[j]
		}
		points = append(points, plotter.XY{X: loops[0].Line[i][0], Y: loops[0].Line[i][1]})
	}

	p := plot.New()

	p.Title.Text = "x vs y"
	p.X.Label.Text = "x"
	p.Y.Label.Text = "y"

	scatter, err := plotter.NewScatter(points)
	if err != nil {
		panic(err)
	}
	scatter.GlyphStyle.Radius = vg.Length(1)
	scatter.GlyphStyle.Shape = draw.CircleGlyph{}
	p.Add(scatter)

	err = p.Save(8*vg.Inch, 8*vg.Inch, "zeta.png")
	if err != nil {
		panic(err)
	}

	return loops
}
