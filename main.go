// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"
	"fmt"
	"image/color"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"time"

	"github.com/mjibson/go-dsp/fft"
	"gonum.org/v1/gonum/integrate/quad"
	"gonum.org/v1/gonum/mathext"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

const (
	// d is the number of dimensions
	d = 3
	// D is the number of space time dimensions
	D = d + 1
	// Lambda is a plate factor
	Lambda = 1.0
	// Genomes is the number of genomes
	Genomes = 1024
)

var (
	// CPUs is the number of CPUs
	CPUs = runtime.NumCPU()
	// FlagGA is the genetic optimization mode
	FlagGA = flag.Bool("ga", false, "ga mode")
	// FlagGraph generate graphs
	FlagGraph = flag.Bool("graph", false, "generate graphs")
	// FlagCompare comares different loop types
	FlagCompare = flag.Bool("compare", false, "compare loop types")
	// Inner is just the inner integration
	FlagInner = flag.Bool("inner", false, "inner integration loop")
)

func square(a float64) float64 {
	return a * a
}

// Worldline is a worldline
type Worldline struct {
	Line   [][d]float64
	Length float64
}

// ComputeLength computes the squared length
func (w *Worldline) ComputeLength() {
	N, length := len(w.Line), 0.0
	for i := 0; i < N+1; i++ {
		v1, v2 := w.Line[(i+N-1)%N], w.Line[i%N]
		for j := 0; j < d; j++ {
			diff := v1[j] - v2[j]
			length += diff * diff
		}
	}
	w.Length = math.Exp(-length / 4)
}

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

// MakeLoopsGA make worldline loops using genetic algorithm
func MakeLoopsGA(N, Loops int) ([]Worldline, [][]complex128) {
	target := -8.9518852409623

	factor := -1 / (2 * math.Pow(4*math.Pi, D/2))
	m := 1.0
	f := func(T float64) float64 {
		return math.Exp(-square(m)*T) / math.Pow(T, 1+D/2)
	}
	factor *= quad.Fixed(f, 0, math.Inf(1), 1000, nil, 0)

	type Genome struct {
		Genome  [][]complex128
		Fitness float64
	}
	genomes := make([]Genome, 0, Genomes)
	rnd := rand.New(rand.NewSource(int64(1)))

	new := func() Genome {
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

		return Genome{
			Genome: y,
		}
	}

	cp := func(genome Genome) Genome {
		cp := Genome{
			Genome:  make([][]complex128, len(genome.Genome)),
			Fitness: genome.Fitness,
		}
		for i := range genome.Genome {
			cp.Genome[i] = make([]complex128, len(genome.Genome[i]))
			copy(cp.Genome[i], genome.Genome[i])
		}
		return cp
	}

	getLoops := func(genome Genome) []Worldline {
		loops := make([]Worldline, 1)

		yt := make([][]complex128, 0, d)
		for i := 0; i < d; i++ {
			genome.Genome[i][0] = 0
			yt = append(yt, fft.FFT(genome.Genome[i]))
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
			loops[0].Line = append(loops[0].Line, point)
		}

		loops[0].ComputeLength()

		return loops
	}

	done := make(chan bool, 8)
	fitness := func(i int) {
		loops := getLoops(genomes[i])

		diff := factor*V(loops, 1, 1) - target
		genomes[i].Fitness = diff * diff

		done <- true
	}

	for i := 0; i < cap(genomes); i++ {
		genomes = append(genomes, new())
	}
	for i := 0; i < 256; i++ {
		start := time.Now()
		j, flight := 0, 0
		for j < len(genomes) && flight < CPUs {
			go fitness(j)
			flight++
			j++
		}
		for j < len(genomes) {
			<-done
			flight--

			go fitness(j)
			j++
			flight++
		}
		for k := 0; k < flight; k++ {
			<-done
		}

		sort.Slice(genomes, func(i, j int) bool {
			return genomes[i].Fitness < genomes[j].Fitness
		})
		fmt.Println(genomes[0].Fitness)
		if genomes[0].Fitness == 0 {
			break
		}
		genomes = genomes[:Genomes]

		for j := 0; j < 10; j++ {
			a, b := cp(genomes[rnd.Intn(10)]), cp(genomes[rnd.Intn(10)])
			x, y, z := rnd.Intn(d), rnd.Intn(N), rnd.Intn(N)
			a.Genome[x][y], b.Genome[x][z] = b.Genome[x][z], a.Genome[x][y]
			genomes = append(genomes, a, b)
		}

		for j := range genomes {
			a := cp(genomes[j])
			x, y := rnd.Intn(d), rnd.Intn(N)
			a.Genome[x][y] += complex(rnd.NormFloat64(), rnd.NormFloat64())
			genomes = append(genomes, a)
		}
		fmt.Println(float64(time.Now().Sub(start)) / float64(time.Second))
	}

	return getLoops(genomes[0]), genomes[0].Genome
}

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

// Plane is a plane
type Plane struct {
	Ox, Oy, Oz    float64
	P1x, P1y, P1z float64
	P2x, P2y, P2z float64
	P3x, P3y, P3z float64
}

func (p Plane) I(a, T float64, loop Worldline, x float64) float64 {
	r1x, r1y, r1z := p.P1x-p.P2x, p.P1y-p.P2y, p.P1z-p.P2z
	r2x, r2y, r2z := p.P3x-p.P2x, p.P3y-p.P2y, p.P3z-p.P2z

	nx := r1y*r2z - r1z*r2y
	ny := r1z*r2x - r1x*r2z
	nz := r1x*r2y - r1y*r2x

	intersections := 0.0
	t := math.Sqrt(T)
	N := len(loop.Line)
	for i := 0; i < N+1; i++ {
		v1, v2 := loop.Line[(i+N-1)%N], loop.Line[i%N]
		l1x, l1y, l1z := x+t*v1[0]-p.Ox, t*v1[1]-p.Oy, t*v1[2]-p.Oz
		l2x, l2y, l2z := x+t*v2[0]-p.Ox, t*v2[1]-p.Oy, t*v2[2]-p.Oz
		rx, ry, rz := l2x-l1x, l2y-l1y, l2z-l1z
		xi := -rx*(nx*l1x+ny*l1y+nz*l1z)/(nx*rx+ny*ry+nz*rz) + l1x
		yi := -ry*(nx*l1x+ny*l1y+nz*l1z)/(nx*rx+ny*ry+nz*rz) + l1y
		zi := -rz*(nx*l1x+ny*l1y+nz*l1z)/(nx*rx+ny*ry+nz*rz) + l1z
		if (xi > l1x && xi < l2x) || (xi > l2x && xi < l1x) ||
			(yi > l1y && yi < l2y) || (yi > l2y && yi < l1y) ||
			(zi > l1z && zi < l2z) || (zi > l2z && zi < l1z) {
			if xi == 0 &&
				yi > p.P2y && yi < p.P1y &&
				zi > p.P3z && zi < p.P2z {
				intersections++
			}
		}
	}

	return Lambda * T * intersections
}

// W integrate over wilson loops
func W(a, T float64, loop Worldline, x float64) float64 {
	intersections := 0.0
	a /= 2
	t := math.Sqrt(T)
	N := len(loop.Line)
	for i := 0; i < N+1; i++ {
		v1, v2 := loop.Line[(i+N-1)%N], loop.Line[i%N]
		if x1, x2 := x+t*v1[0], x+t*v2[0]; (x1 < -a && x2 > a) ||
			(x1 > a && x2 < -a) {
			intersections += 2
		} else if (x1 < a && x2 > a) ||
			(x1 > a && x2 < a) ||
			(x1 < -a && x2 > -a) ||
			(x1 > -a && x2 < -a) {
			intersections++
		}
	}

	return Lambda * T * intersections
}

// Result is a wilson loop integration
type Result struct {
	Intersections float64
	Length        float64
}

// V compute energy for plate separation a
func V(loops []Worldline, a, T float64) float64 {
	w := func(x float64) float64 {
		done := make(chan Result, 8)
		process := func(a, T float64, i int) {
			intersections := W(a, T, loops[i], x)
			done <- Result{
				Intersections: intersections,
				Length:        loops[i].Length,
			}
		}
		length, sum, denominator, flight, i := len(loops), 0.0, 0.0, 0, 0
		for i < length && flight < CPUs {
			go process(a, T, i)
			flight++
			i++
		}
		for i < length {
			result := <-done
			flight--
			sum += math.Exp(-result.Intersections) * result.Length
			denominator += result.Length

			go process(a, T, i)
			flight++
			i++
		}
		for i := 0; i < flight; i++ {
			result := <-done
			sum += math.Exp(-result.Intersections) * result.Length
			denominator += result.Length
		}
		return (sum/denominator - 1)
	}
	b := a
	if T > 1 {
		b *= T
	}
	return quad.Fixed(w, -b, b, 1000, nil, 0)
}

func main() {
	flag.Parse()

	fmt.Println("CPUs=", CPUs)

	if *FlagGA {
		loops, freq := MakeLoopsGA(1024, 1)
		N := len(loops[0].Line)

		points := make(plotter.XYs, 0, N)
		for n, value := range loops[0].Line {
			points = append(points, plotter.XY{X: float64(n), Y: value[0]})
		}

		p := plot.New()

		p.Title.Text = "x vs time"
		p.X.Label.Text = "x"
		p.Y.Label.Text = "time"

		scatter, err := plotter.NewScatter(points)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		p.Add(scatter)

		err = p.Save(8*vg.Inch, 8*vg.Inch, "loop.png")
		if err != nil {
			panic(err)
		}

		points = make(plotter.XYs, 0, N)
		for n, value := range freq[0] {
			points = append(points, plotter.XY{X: float64(n), Y: real(value)})
		}

		p = plot.New()

		p.Title.Text = "freq vs time"
		p.X.Label.Text = "freq"
		p.Y.Label.Text = "time"

		scatter, err = plotter.NewScatter(points)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		p.Add(scatter)

		err = p.Save(8*vg.Inch, 8*vg.Inch, "freq.png")
		if err != nil {
			panic(err)
		}

		points = make(plotter.XYs, 0, N)
		for _, value := range freq[0] {
			points = append(points, plotter.XY{X: real(value), Y: imag(value)})
		}

		p = plot.New()

		p.Title.Text = "real vs imag"
		p.X.Label.Text = "real"
		p.Y.Label.Text = "image"

		scatter, err = plotter.NewScatter(points)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		p.Add(scatter)

		err = p.Save(8*vg.Inch, 8*vg.Inch, "complex.png")
		if err != nil {
			panic(err)
		}

		return
	}

	if *FlagCompare {
		factor := -1 / (2 * math.Pow(4*math.Pi, D/2))
		m := Lambda / 100

		compare := func(t func(N, Loops int) []Worldline) plotter.XYs {
			points := make(plotter.XYs, 0, 8)
			size := 2
			for i := 0; i < 13; i++ {
				loops := t(size, size)
				for i := range loops {
					loops[i].ComputeLength()
				}

				a := 1.0
				f := func(T float64) float64 {
					return math.Exp(-square(m)*T) / math.Pow(T, 1+D/2) * V(loops, a, T)
				}
				e := factor * quad.Fixed(f, 1/square(1e15), math.Inf(1), 200, nil, 0)
				if math.IsNaN(e) {
					e = 0
				}
				fmt.Println(size, e)
				points = append(points, plotter.XY{X: float64(size), Y: e})

				size *= 2
			}
			return points
		}

		p := plot.New()

		p.Title.Text = "size vs energy"
		p.X.Label.Text = "size"
		p.Y.Label.Text = "energy"

		points := compare(MakeLoopsFFT)
		scatter, err := plotter.NewScatter(points)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.Color = color.RGBA{0xFF, 0, 0, 0xFF}
		p.Add(scatter)

		points = compare(MakeLoopsFFT2)
		scatter, err = plotter.NewScatter(points)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.Color = color.RGBA{0, 0, 0xFF, 0xFF}
		p.Add(scatter)

		err = p.Save(8*vg.Inch, 8*vg.Inch, "size.png")
		if err != nil {
			panic(err)
		}

		return
	}

	if *FlagInner {
		loops := MakeLoopsFFT2(1024, 1024)
		for i := range loops {
			loops[i].ComputeLength()
		}

		/*
			line := func(a, T float64, loop Worldline, x float64) float64 {
				intersections := 0.0
				a /= 2
				t := math.Sqrt(T)
				N := len(loop.Line)
				for i := 0; i < N+1; i++ {
					v1, v2 := loop.Line[(i+N-1)%N], loop.Line[i%N]
					if x1, x2 := x+t*v1[0], x+t*v2[0]; (x1 < a && x2 > a) ||
						(x1 > a && x2 < a) {
						intersections++
					}
				}

				return Lambda * T * intersections
			}*/

		p1 := Plane{
			0.5, 0, 0,
			0.0, .5, .5,
			0.0, -.5, .5,
			0.0, -.5, -.5,
		}

		circle := func(a, T float64, loop Worldline, x float64) float64 {
			intersections := 0.0
			t := math.Sqrt(T)
			N := len(loop.Line)
			r := a / 2
			rSquared := square(r)
			for i := 0; i < N+1; i++ {
				v1, v2 := loop.Line[(i+N-1)%N], loop.Line[i%N]
				r1 := square(x+t*v1[0]-a) + square(v1[1])
				r2 := square(x+t*v2[0]-a) + square(v2[1])
				if (r1 < rSquared && r2 > rSquared) ||
					(r1 > rSquared && r2 < rSquared) {
					intersections++
				} else {
					// https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
					Ax, Ay := t*v1[0]+x, v1[1]
					Bx, By := t*v2[0]+x, v2[1]
					Cx := a
					Cy := 0.0

					// compute the euclidean distance between A and B
					LAB := math.Sqrt(square(Bx-Ax) + square(By-Ay))

					// compute the direction vector D from A to B
					Dx := (Bx - Ax) / LAB
					Dy := (By - Ay) / LAB

					// the equation of the line AB is x = Dx*t + Ax, y = Dy*t + Ay with 0 <= t <= LAB.

					// compute the distance between the points A and E, where
					// E is the point of AB closest the circle center (Cx, Cy)
					t := Dx*(Cx-Ax) + Dy*(Cy-Ay)

					// compute the coordinates of the point E
					Ex := t*Dx + Ax
					Ey := t*Dy + Ay

					// compute the euclidean distance between E and C
					LEC := math.Sqrt(square(Ex-Cx) + square(Ey-Cy))

					if LEC < r {
						// compute distance from t to circle intersection point
						dt := math.Sqrt(rSquared - square(LEC))

						// compute first intersection point
						Fx := (t-dt)*Dx + Ax
						Fy := (t-dt)*Dy + Ay

						// compute second intersection point
						Gx := (t+dt)*Dx + Ax
						Gy := (t+dt)*Dy + Ay

						if ((Fx > Ax && Fx < Bx && Fy > Ay && Fy < By) ||
							(Fx > Bx && Fx < Ax && Fy > By && Fy < Ay)) &&
							((Gx > Ax && Gx < Bx && Gy > Ay && Gy < By) ||
								(Gx > Bx && Gx < Ax && Gy > By && Gy < Ay)) {
							intersections += 2
						}
					}
				}
			}
			return intersections
		}

		pointsLine := make(plotter.XYs, 0, 10)
		for x := 0; x < 200; x++ {
			x, intersections := float64(x)*.01, 0.0
			for _, loop := range loops {
				intersections += p1.I(1, 1, loop, x)
			}
			pointsLine = append(pointsLine, plotter.XY{X: x, Y: intersections})
		}

		pointsCircle := make(plotter.XYs, 0, 10)
		for x := 0; x < 200; x++ {
			x, intersections := float64(x)*.01, 0.0
			for _, loop := range loops {
				intersections += circle(1, 1, loop, x)
			}
			pointsCircle = append(pointsCircle, plotter.XY{X: x, Y: intersections})
		}

		pointsLineT := make(plotter.XYs, 0, 10)
		for x := -100; x < 300; x++ {
			x, intersections := float64(x)*.01, 0.0
			for _, loop := range loops {
				intersections += p1.I(1, 4, loop, x)
			}
			pointsLineT = append(pointsLineT, plotter.XY{X: x, Y: intersections})
		}

		pointsCircleT := make(plotter.XYs, 0, 10)
		for x := -100; x < 300; x++ {
			x, intersections := float64(x)*.01, 0.0
			for _, loop := range loops {
				intersections += circle(1, 4, loop, x)
			}
			pointsCircleT = append(pointsCircleT, plotter.XY{X: x, Y: intersections})
		}

		p := plot.New()

		p.Title.Text = "x vs intersections"
		p.X.Label.Text = "x"
		p.Y.Label.Text = "intersections"

		scatter, err := plotter.NewScatter(pointsLine)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.Color = color.RGBA{0xFF, 0, 0, 0xFF}
		p.Add(scatter)

		scatter, err = plotter.NewScatter(pointsCircle)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.Color = color.RGBA{0, 0, 0xFF, 0xFF}
		p.Add(scatter)

		scatter, err = plotter.NewScatter(pointsLineT)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.Color = color.RGBA{0, 0, 0, 0xFF}
		p.Add(scatter)

		scatter, err = plotter.NewScatter(pointsCircleT)
		if err != nil {
			panic(err)
		}
		scatter.GlyphStyle.Radius = vg.Length(1)
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.Color = color.RGBA{0xFF, 0, 0xFF, 0xFF}
		p.Add(scatter)

		err = p.Save(8*vg.Inch, 8*vg.Inch, "inner.png")
		if err != nil {
			panic(err)
		}
		return
	}

	fmt.Println("making loops...")
	//loops := MakeLoopsFFT(1024, 1024)
	loops := MakeLoopsFFT2(1024, 1024)
	//loops := MakeLoopsZeta(1024, 1024)

	for i := range loops {
		loops[i].ComputeLength()
	}

	fmt.Println("simulating...")
	factor := -1 / (2 * math.Pow(4*math.Pi, D/2))
	m := Lambda / 100

	points := make(plotter.XYs, 0, 10)
	for i := 1; i < 100; i++ {
		a := float64(i)
		f := func(T float64) float64 {
			return math.Exp(-square(m)*T) / math.Pow(T, 1+D/2) * V(loops, a, T)
		}
		e := factor * quad.Fixed(f, 1/square(1e15), math.Inf(1), 200, nil, 0)
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
