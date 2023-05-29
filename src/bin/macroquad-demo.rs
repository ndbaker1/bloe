use macroquad::prelude::*;

use bloe;

#[allow(dead_code)]
enum Mode {
    OneCircle,
    Scene,
}

#[allow(dead_code)]
enum Plot {
    Density,
    Velocity,
    Vorticity,
}

const MODE: Mode = Mode::OneCircle;
const PLOT: Plot = Plot::Vorticity;

const XDIM: usize = 200;
const YDIM: usize = 100;

// initial external forces
const FORCE_APPLIED: f32 = 30.0;

#[macroquad::main("LBM Simulator")]
async fn main() {
    let width = screen_width() / XDIM as f32;
    let height = screen_height() / YDIM as f32;

    let mut sim = bloe::LBM::<XDIM, YDIM>::new();

    // random natural force over the entire field
    for i in 0..XDIM {
        for j in 0..YDIM {
            for k in 0..bloe::NDIR {
                sim.f[i][j][k] += 0.001 * rand::gen_range(0.0, 1.0);
            }
        }
    }

    const B: usize = 30;
    for j in B..YDIM - B {
        // increase force in the right (->) direction
        // along the 1st column
        sim.f[0][j][1] += FORCE_APPLIED * (1.0 + 0.2 * rand::gen_range(0.0, 1.0));
    }

    // initialize the values of the lattice field
    // using a given base average density
    for i in 0..XDIM {
        for j in 0..YDIM {
            // recompute a new average density in the field
            let rho: f32 = sim.f[i][j].iter().sum();
            for k in 0..bloe::NDIR {
                sim.f[i][j][k] *= 1.0 / rho;
            }
        }
    }

    // simple circular boundary to use for testing
    struct Circle<const R: isize, const X: isize, const Y: isize>;
    impl<const RADIUS: isize, const X: isize, const Y: isize> bloe::Boundary for Circle<RADIUS, X, Y> {
        fn contains(&self, x: usize, y: usize) -> bool {
            let dx = X - x as isize;
            let dy = Y - y as isize;
            dx * dx + dy * dy < RADIUS.pow(2)
        }
    }

    struct Square<const S: isize, const X: isize, const Y: isize>;
    impl<const SIDE: isize, const X: isize, const Y: isize> bloe::Boundary for Square<SIDE, X, Y> {
        fn contains(&self, x: usize, y: usize) -> bool {
            let dx = X - x as isize;
            let dy = Y - y as isize;
            dx * dx < SIDE.pow(2) && dy * dy < SIDE.pow(2)
        }
    }

    match MODE {
        Mode::Scene => {
            // place a boundary to block the flow
            sim.boundaries.push(&Circle::<20, 140, 30>);
            sim.boundaries.push(&Square::<25, 90, 60>);
            sim.boundaries.push(&Circle::<10, 40, 50>);
        }
        Mode::OneCircle => {
            const XDIM_M: isize = XDIM as isize / 2;
            const YDIM_M: isize = YDIM as isize / 2;
            sim.boundaries.push(&Circle::<10, XDIM_M, YDIM_M>);
        }
    };

    loop {
        clear_background(WHITE);

        sim.run(1);

        for i in 0..XDIM {
            for j in 0..YDIM {
                let w = match PLOT {
                    Plot::Density => sim.rho[i][j],
                    Plot::Velocity => (sim.u[i][j].0.powi(2) + sim.u[i][j].1.powi(2)).sqrt() * 10.0,
                    Plot::Vorticity => vorticity(&sim.u, i, j) * 10.0,
                };

                let color = match sim.boundaries.iter().any(|b| b.contains(i, j)) {
                    true => BROWN,
                    false => Color {
                        r: w,
                        g: w,
                        b: w,
                        a: 1.0,
                    },
                };

                draw_rectangle(i as f32 * width, j as f32 * height, width, height, color);
            }
        }

        next_frame().await;
    }
}

fn vorticity(u_mat: &[[(f32, f32); YDIM]], i: usize, j: usize) -> f32 {
    let x_f = u_mat[(i + 1).rem_euclid(XDIM)][j].0;
    let x_i = u_mat[(i as isize - 1).rem_euclid(XDIM as isize) as usize][j].0;

    let y_f = u_mat[i][(j + 1).rem_euclid(YDIM)].1;
    let y_i = u_mat[i][(j as isize - 1).rem_euclid(YDIM as isize) as usize].1;

    let vx = x_f - x_i;
    let vy = y_f - y_i;

    vy - vx
}
