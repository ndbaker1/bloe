/// Lattice Boltzmann Method Simulator for Computational Fluid Dynamics (CFD)
///
/// * applies particle bound-back to boundaries
pub struct LBM<'b, const XDIM: usize, const YDIM: usize> {
    /// collision timescale parameter
    pub tau: f64,
    /// density at each lattice
    pub rho: [[f64; YDIM]; XDIM],
    /// velocity at each lattice
    pub u: [[(f64, f64); YDIM]; XDIM],
    /// boundaries where fluid cannot pass or exist
    pub boundaries: Vec<&'b dyn Boundary>,
    /// lattice field
    pub f: [[[f64; NDIR]; YDIM]; XDIM],
}

trait Dot<T> {
    /// Compute the dot product of two vector-like objects
    fn dot(self, other: Self) -> T;
}

impl<T> Dot<T> for (T, T)
where
    T: std::ops::Mul<Output = T> + std::ops::Add<Output = T>,
{
    fn dot(self, other: Self) -> T {
        self.0 * other.0 + self.1 * other.1
    }
}

/// Figures which can be used to represent a boundary within the
/// Lattice Boltzmann Method simulation
pub trait Boundary {
    /// determine whether a point in within the boundary of the figure
    fn contains(&self, x: usize, y: usize) -> bool;
}

/// Number of points that exist in the 2D lattice based on the D2Q9 scheme
const NDIR: usize = 9;

impl<'b, const X: usize, const Y: usize> LBM<'b, X, Y> {
    pub fn new(rho_0: f64) -> Self {
        let mut lbm = Self {
            tau: 0.53,
            rho: Self::create_field(rho_0),
            u: Self::create_field((0.0, 0.0)),
            f: Self::create_field([1.0; NDIR]),
            boundaries: Vec::default(),
        };

        // initialize the values of the lattice field
        // using a given base average density
        for i in 0..X {
            for j in 0..Y {
                let rho: f64 = lbm.f[i][j].iter().sum();
                for k in 0..NDIR {
                    lbm.f[i][j][k] *= rho_0 / rho;
                }
            }
        }

        lbm.update_macro();

        lbm
    }

    #[inline]
    fn create_field<T: Copy>(fill: T) -> [[T; Y]; X] {
        [[fill; Y]; X]
    }

    /// execute a given number of steps of the simulation
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.clean_boundaries();

            self.update_macro();
            self.collide();
            self.stream();
        }
    }

    /// zero-out values within boundaries before performing any additional computations
    fn clean_boundaries(&mut self) {
        for i in 0..X {
            for j in 0..Y {
                if self.boundaries.iter().any(|b| b.contains(i, j)) {
                    self.f[i][j] = [0.0; NDIR];
                    self.u[i][j] = (0.0, 0.0);
                    self.rho[i][j] = 0.0;
                }
            }
        }
    }

    /// stream step of the Lattice Boltzmann Method,
    /// where the distributions of the outer surrounding lattice sites
    /// are propagated to neighboring latticies.
    fn stream(&mut self) {
        let mut f_new = Self::create_field([0.0; NDIR]);

        for k in 0..NDIR {
            let (vx, vy) = SITE_VECS[k];

            for i in 0..X {
                for j in 0..Y {
                    // ignore processing points inside of boundaries,
                    // where fluid is unable to move (or doesn't exist)
                    if self.boundaries.iter().any(|b| b.contains(i, j)) {
                        continue;
                    }

                    let i_new = (i as f64 + vx).rem_euclid(X as _);
                    let j_new = (j as f64 + vy).rem_euclid(Y as _);

                    let i_new = i_new as usize;
                    let j_new = j_new as usize;

                    // use particle bounce-back when computing streaming operations
                    // that collide with a boundary.
                    // this involves reverse the vector within the current lattice (i,j),
                    // rather than propagating to the point within the boundary (i_new, j_new).
                    if self.boundaries.iter().any(|b| b.contains(i_new, j_new)) {
                        f_new[i][j][SITE_REV[k]] = self.f[i][j][k];
                    } else {
                        f_new[i_new][j_new][k] = self.f[i][j][k];
                    }
                }
            }
        }

        self.f = f_new;
    }

    /// collision step of the Lattice Boltzmann Method,
    /// in which the lattice is updated using the equilibirum state formula:
    ///
    /// ```math
    /// F_eq = w * ρ * ( 1 + (3)(c⋅u) + (3/2)(c⋅u)^2 - (9/2)(u⋅u) )
    /// ```
    fn collide(&mut self) {
        for i in 0..X {
            for j in 0..Y {
                for k in 0..NDIR {
                    let c = SITE_VECS[k];
                    let u = self.u[i][j];

                    let p1 = 3.0 * c.dot(u);
                    // p2 = 9 / 2 * c.dot(u) ^ 2
                    let p2 = p1.powi(2) / 2.0;
                    let p3 = (3.0 / 2.0) * u.dot(u);

                    // equilibrium state lattice point
                    let f_eq_ijk = self.rho[i][j] * WEIGHTS[k] * (1.0 + p1 + p2 - p3);

                    // update the lattice site using the LBM step
                    self.f[i][j][k] += -(self.f[i][j][k] - f_eq_ijk) / self.tau;
                }
            }
        }
    }

    /// recompute the values of `ρ` (density) and `u` (velocity)
    /// based on the current state of the lattice field.
    ///
    /// ## Note
    /// if any value or parameters of the simulator are manually changed
    /// outside the jurisdiction of the [`LBM::run()`], then this function
    /// should be explicitly called.
    fn update_macro(&mut self) {
        for i in 0..X {
            for j in 0..Y {
                // density is the sum of the particle distribution within the lattice
                let rho: f64 = self.f[i][j].iter().sum();

                // velocity of the lattice is the weighted sum of each 9 lattice sites
                // incorporating their respective velocities
                let x: f64 = self.f[i][j]
                    .iter()
                    .enumerate()
                    .map(|(k, s)| SITE_VECS[k].0 * s)
                    .sum();

                let y: f64 = self.f[i][j]
                    .iter()
                    .enumerate()
                    .map(|(k, s)| SITE_VECS[k].1 * s)
                    .sum();

                self.rho[i][j] = rho;
                self.u[i][j] = if rho == 0.0 {
                    (0.0, 0.0)
                } else {
                    (x / rho, y / rho)
                };
            }
        }
    }
}

/// coordinates describing positions of sites in the D2Q9 lattice scheme
const SITE_VECS: [(f64, f64); 9] = [
    (0.0, 0.0),   // center
    (1.0, 0.0),   // right
    (1.0, 1.0),   // bottom-right
    (0.0, 1.0),   // bottom
    (-1.0, 1.0),  // bottom-left
    (-1.0, 0.0),  // left
    (-1.0, -1.0), // top-left
    (0.0, -1.0),  // top
    (1.0, -1.0),  // top-right
];

/// weights of sites in the D2Q9 lattice scheme
const WEIGHTS: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 9.0,
    1.0 / 36.0,
];

/// reversed lattice directions used for reflective boundaries
const SITE_REV: [usize; 9] = [0, 5, 6, 7, 8, 1, 2, 3, 4];

#[cfg(test)]
mod test {
    use super::*;

    /// This is an example of the workflow for a simple single-phase flow
    /// with an obstacle.
    #[test]
    fn sample() {
        const XDIM: usize = 9;
        const YDIM: usize = 9;

        let mut sim = LBM::<XDIM, YDIM>::new(100.0);

        // add initial external forces
        for i in 0..XDIM {
            for j in 0..YDIM {
                sim.f[i][j][1] *= 2.0;
            }
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
        // (circle with a radius of 1)
        sim.boundaries.push(&Circle::<1>);

        let average_density = sim.rho[0][0];
        let expected_mass = (XDIM * YDIM - 1) as f64 * average_density;

        // after running the simulation with 0 steps nothing should change,
        // which includes the fact that no boundaries will be check and the
        // distributions should be completely smooth
        sim.run(0);

        // the center should be zero'd out in terms of density and velocity
        // once we step the simulation containing a placed boundary
        sim.run(1);

        assert_eq!(sim.rho[XDIM / 2][YDIM / 2], 0.0);
        assert_eq!(sim.u[XDIM / 2][YDIM / 2], (0.0, 0.0));
        for k in sim.f[XDIM / 2][YDIM / 2] {
            assert_eq!(k, 0.0);
        }

        // continue running the simulation and checking for stable behavior
        sim.run(200);

        // display the rho values of the field
        let mut mass = 0.0;

        println!("density grid:");
        for row in sim.rho {
            for density in row {
                print!("{:8.2?}", density);
                mass += density;
            }
            println!("");
        }

        println!("mass: {}", mass);
        // division can drift. allow rounding
        assert_eq!(mass.round(), expected_mass.round());

        fn round_to(v: f64, p: u32) -> f64 {
            (v * 10u32.pow(p) as f64).round()
        }

        // verify that the system has reached a stable state with no motion
        let mut final_it = -1;
        const PRECISION: u32 = 3;
        for it in 0..10000 {
            sim.run(1);
            if sim.u.iter().all(|r| {
                r.iter()
                    .all(|u| round_to(u.0, PRECISION) == 0.0 && round_to(u.1, PRECISION) == 0.0)
            }) {
                final_it = it;
                break;
            }
        }

        println!("convergence after: {}", final_it);
        assert_ne!(final_it, -1);
    }
}
