// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"
	"fmt"
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
	d = 2
	// D is the number of space time dimensions
	D = d + 1
	// N is the number of points in the Worldline
	N = 1024
	// Loops is the number of loops
	Loops = 1024
	// Lambda is a plate factor
	Lambda = 1e-3
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
)

// MakeLoopsZeta makes a loop using the zeta function
func MakeLoopsZeta() [][][d]float64 {
	loops := make([][][d]float64, 1)
	x := 2.0
	for i := 0; i < N; i++ {
		var point [d]float64
		for j := 0; j < d; j++ {
			point[j] = mathext.Zeta(float64(x), float64(x))
		}
		x += .01
		loops[0] = append(loops[0], point)
	}

	var min [d]float64
	var max [d]float64
	for i := 0; i < d; i++ {
		min[i] = math.MaxFloat64
		max[i] = -math.MaxFloat64
	}
	for i := 0; i < N; i++ {
		for j := 0; j < d; j++ {
			r := loops[0][i][j]
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
			loops[0][i][j] /= norm[j]
		}
		points = append(points, plotter.XY{X: loops[0][i][0], Y: loops[0][i][1]})
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
func MakeLoopsGA() ([][][d]float64, [][]complex128) {
	target := -8.9518852409623

	factor := -1 / (2 * math.Pow(4*math.Pi, D/2))
	m := 1.0
	f := func(T float64) float64 {
		return math.Exp(-math.Pow(m, 2)*T) / math.Pow(T, 1+D/2)
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

	getLoops := func(genome Genome) [][][d]float64 {
		loops := make([][][d]float64, 1)

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
			loops[0] = append(loops[0], point)
		}

		return loops
	}

	done := make(chan bool, 8)
	fitness := func(i int) {
		loops := getLoops(genomes[i])

		diff := factor*V(loops, 1) - target
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
func MakeLoopsMulti() [][][d]float64 {
	loops := make([][][d]float64, Loops)
	rnd := rand.New(rand.NewSource(int64(1)))
	y := make([][][]complex128, 0, d)
	for i := 0; i < d; i++ {
		y = append(y, make([][]complex128, 0, Loops))
		for j := 0; j < Loops; j++ {
			y[i] = append(y[i], make([]complex128, 0, N))
		}
	}
	factor := math.Sqrt(2)
	for i := 0; i < d; i++ {
		for j := 0; j < Loops; j++ {
			for k := 0; k < N; k++ {
				if j == 0 && k == 0 {
					y[i][j] = append(y[i][j], 0)
					continue
				}
				y[i][j] = append(y[i][j], complex(rnd.Float64()*factor, rnd.Float64()*factor))
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
			loops[loop] = append(loops[loop], point)
		}
	}

	if *FlagGraph {
		for loop := 0; loop < Loops; loop++ {
			points := make(plotter.XYs, 0, N)
			for i := 0; i < N; i++ {
				points = append(points, plotter.XY{X: loops[loop][i][0], Y: loops[loop][i][1]})
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

// W integrate over wilson loops
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

// Result is a wilson loop integration
type Result struct {
	Intersections float64
	Length        float64
}

// V compute energy for plate separation a
func V(loops [][][d]float64, a float64) float64 {
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
		for i < len(loops) && flight < CPUs {
			go process(a, i)
			flight++
			i++
		}
		for i < len(loops) {
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
		return (sum/denominator - 1)
	}
	return quad.Fixed(w, -1, 1, 1000, nil, 0)
}

func main() {
	flag.Parse()

	if *FlagGA {
		loops, freq := MakeLoopsGA()

		points := make(plotter.XYs, 0, N)
		for n, value := range loops[0] {
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

	fmt.Println("making loops...")
	loops := MakeLoopsMulti()
	//loops := MakeLoopsZeta()

	fmt.Println("simulating...")
	factor := -1 / (2 * math.Pow(4*math.Pi, D/2))
	m := 1.0
	f := func(T float64) float64 {
		return math.Exp(-math.Pow(m, 2)*T) / math.Pow(T, 1+D/2)
	}
	factor *= quad.Fixed(f, 0, math.Inf(1), 1000, nil, 0)

	points := make(plotter.XYs, 0, 10)
	for a := 1.0; a > 0; a -= .01 {
		e := factor * V(loops, a)
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
