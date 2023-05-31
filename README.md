<p align="center">
	<h1 align="center">💨 Bloe</h1>
</p>
<p align="center">Lattice Boltzmann Method Simulator in Rust</p>

<h3 align="center">
	<a href="https://ndbaker1.github.io/bloe/">
		Online Demo ☁
	</a>
</h2>

## Usage ⚡

To run the code locally, clone the project and run:

```bash
cargo run
```

This will execute the [visualization demo](./src/bin/macroquad-demo.rs) locally

## Configuration ⚙

Several parameters/constants of the simulation can be changed including:
 
* Lattice Grid Dimensions
* Lattice Density Distribution (to influence directional Flow)
* Boundary Objects
* Simulation Step Size

## Motivation 🤔

A dead-simple implementation of the lattice Boltzmann method with as few dependencies as possible in my favorite programming language 🦀.

> The lattice Boltzmann method is an application of Boltzmann particle (microscopic) principles in a lattice grid to simulate computation fluid dynamics (CFD) without directly solving the Navier-Stokes equations (macroscopic)
