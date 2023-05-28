use bloe::{Boundary, LBM};
use macroquad::prelude::*;

const XDIM: usize = 200;
const YDIM: usize = 200;

#[macroquad::main("LBM Sim")]
async fn main() {
    let width = screen_width() / XDIM as f32;
    let height = screen_height() / YDIM as f32;

    // initial average density
    let density_0 = 500.0;
    let max_compression_factor = 2.0;

    let mut sim = LBM::<XDIM, YDIM>::new(density_0);

    // add initial external forces
    for j in 0..YDIM {
        sim.f[0][j][1] *= 15.0;
    }

    // update the macro parameters
    sim.update_macro();

    // simple circular boundary to use for testing
    struct Circle<const X: isize>;
    impl<const RADIUS: isize> Boundary for Circle<RADIUS> {
        fn contains(&self, x: usize, y: usize) -> bool {
            let dx = XDIM as isize / 2 - x as isize;
            let dy = YDIM as isize / 2 - y as isize;
            dx * dx + dy * dy < RADIUS.pow(2)
        }
    }

    // place a boundary to block the flow
    sim.boundaries.push(&Circle::<30>);

    loop {
        clear_background(WHITE);

        sim.run(1);

        for i in 0..XDIM {
            for j in 0..YDIM {
                let density = sim.rho[i][j];
                let color = if density == 0.0 {
                    BLACK
                } else {
                    Color {
                        r: 0.0,
                        b: density / (density_0 * max_compression_factor),
                        g: 0.0,
                        a: 1.0,
                    }
                };

                draw_rectangle(i as f32 * width, j as f32 * height, width, height, color);
            }
        }

        next_frame().await;
    }
}
