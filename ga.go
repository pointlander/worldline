// Copyright 2022 The Worldline Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"time"

	"github.com/mjibson/go-dsp/fft"
	"gonum.org/v1/gonum/integrate/quad"
)

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
