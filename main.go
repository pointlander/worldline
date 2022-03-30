// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"
	"fmt"
	"image/color"
	"math"
	"runtime"

	"gonum.org/v1/gonum/integrate/quad"
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

// https://www.mathsisfun.com/algebra/vectors-cross-product.html
func crossProduct(ax, ay, az, bx, by, bz float64) (cx, cy, cz float64) {
	cx = ay*bz - az*by
	cy = az*bx - ax*bz
	cz = ax*by - ay*bx
	return
}

func min(a ...float64) float64 {
	min := math.MaxFloat64
	for _, value := range a {
		if value < min {
			min = value
		}
	}
	return min
}

func max(a ...float64) float64 {
	max := -math.MaxFloat64
	for _, value := range a {
		if value > max {
			max = value
		}
	}
	return max
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

// Plane is a plane
type Plane struct {
	Ox, Oy, Oz    float64
	P1x, P1y, P1z float64
	P2x, P2y, P2z float64
}

// I calculates worldline intersections with a plane
// https://johannesbuchner.github.io/intersection/intersection_line_plane.html
// https://www.math.usm.edu/lambers/mat169/fall09/lecture25.pdf
// https://tutorial.math.lamar.edu/classes/calciii/eqnsofplanes.aspx
func (p Plane) I(a, T float64, loop Worldline, x, y, z float64) float64 {
	nx, ny, nz := crossProduct(p.P1x, p.P1y, p.P1z, p.P2x, p.P2y, p.P2z)
	intersections := 0.0
	t := math.Sqrt(T)
	N := len(loop.Line)
	for i := 0; i < N+1; i++ {
		v1, v2 := loop.Line[(i+N-1)%N], loop.Line[i%N]
		l1x, l1y, l1z := x+t*v1[0]-p.Ox, y+t*v1[1]-p.Oy, z+t*v1[2]-p.Oz
		l2x, l2y, l2z := x+t*v2[0]-p.Ox, y+t*v2[1]-p.Oy, z+t*v2[2]-p.Oz
		rx, ry, rz := l2x-l1x, l2y-l1y, l2z-l1z
		numerator, denominator := (nx*l1x + ny*l1y + nz*l1z), (nx*rx + ny*ry + nz*rz)
		if denominator == 0 {
			continue
		}
		xi := l1x - rx*numerator/denominator
		yi := l1y - ry*numerator/denominator
		zi := l1z - rz*numerator/denominator
		if ((l1x == l2x) || (xi > l1x && xi < l2x) || (xi > l2x && xi < l1x)) &&
			((l1y == l2y) || (yi > l1y && yi < l2y) || (yi > l2y && yi < l1y)) &&
			((l1z == l2z) || (zi > l1z && zi < l2z) || (zi > l2z && zi < l1z)) {
			p3x, p3y, p3z := -p.P1x, -p.P1y, -p.P1z
			p4x, p4y, p4z := -p.P2x, -p.P2y, -p.P2z
			minx, maxx := min(p.P1x, p.P2x, p3x, p4x), max(p.P1x, p.P2x, p3x, p4x)
			miny, maxy := min(p.P1y, p.P2y, p3y, p4y), max(p.P1y, p.P2y, p3y, p4y)
			minz, maxz := min(p.P1z, p.P2z, p3z, p4z), max(p.P1z, p.P2z, p3z, p4z)
			if xi == 0 &&
				yi > miny && yi < maxy &&
				zi > minz && zi < maxz {
				intersections++
			} else if yi == 0 &&
				xi > minx && xi < maxx &&
				zi > minz && zi < maxz {
				intersections++
			} else if zi == 0 &&
				xi > minx && xi < maxx &&
				yi > miny && yi < maxy {
				intersections++
			} else if xi > minx && xi < maxx &&
				yi > miny && yi < maxy &&
				zi > minz && zi < maxz {
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

		p1 := Plane{
			0.5, 0, 0,
			0.0, .5, .5,
			0.0, -.5, .5,
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
				intersections += p1.I(1, 1, loop, x, 0, 0)
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
				intersections += p1.I(1, 4, loop, x, 0, 0)
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
