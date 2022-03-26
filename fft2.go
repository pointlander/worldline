// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"math/rand"

	"github.com/mjibson/go-dsp/fft"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

// MakeLoopsMulti make worldline loops
func MakeLoopsFFT2(N, Loops int) []Worldline {
	loops := make([]Worldline, Loops)
	rnd := rand.New(rand.NewSource(int64(1)))
	y := make([][][]complex128, 0, d)
	for i := 0; i < d; i++ {
		y = append(y, make([][]complex128, 0, Loops))
		for j := 0; j < Loops; j++ {
			y[i] = append(y[i], make([]complex128, 0, N))
		}
	}
	for i := 0; i < d; i++ {
		for j := 0; j < Loops; j++ {
			for k := 0; k < N; k++ {
				if j == 0 && k == 0 {
					y[i][j] = append(y[i][j], 0)
					continue
				}
				y[i][j] = append(y[i][j], complex(rnd.NormFloat64(), rnd.NormFloat64()))
			}
		}
	}

	yt := make([][][]complex128, 0, d)
	for i := 0; i < d; i++ {
		yt = append(yt, fft.FFT2(y[i]))
	}

	min := make([][d]float64, Loops)
	max := make([][d]float64, Loops)
	for loop := 0; loop < Loops; loop++ {
		for i := 0; i < d; i++ {
			min[loop][i] = math.MaxFloat64
			max[loop][i] = -math.MaxFloat64
		}
		for i := 0; i < N; i++ {
			if loop == 0 && i == 0 {
				continue
			}
			for j := 0; j < d; j++ {
				r := real(yt[j][loop][i])
				if r < min[loop][j] {
					min[loop][j] = r
				}
				if r > max[loop][j] {
					max[loop][j] = r
				}
			}
		}
	}

	norm := make([][d]float64, Loops)
	for loop := 0; loop < Loops; loop++ {
		for i := 0; i < d; i++ {
			norm[loop][i] = math.Abs(max[loop][i] - min[loop][i])
		}
	}

	for loop := 0; loop < Loops; loop++ {
		for i := 0; i < N; i++ {
			var point [d]float64
			for j := 0; j < d; j++ {
				if loop == 0 && i == 0 {
					yt[j][loop][i] = 0
				}
				point[j] = real(yt[j][loop][i]) / norm[loop][j]
			}
			loops[loop].Line = append(loops[loop].Line, point)
		}
	}

	if *FlagGraph {
		for loop := 0; loop < Loops; loop++ {
			points := make(plotter.XYs, 0, N)
			for i := 0; i < N; i++ {
				points = append(points, plotter.XY{X: loops[loop].Line[i][0], Y: loops[loop].Line[i][1]})
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

			err = p.Save(8*vg.Inch, 8*vg.Inch, fmt.Sprintf("fft_%d.png", loop))
			if err != nil {
				panic(err)
			}
		}
	}

	return loops
}
