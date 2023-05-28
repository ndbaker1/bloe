/// Lattice Boltzmann Method Simulator
pub struct LBM<'b, const XDIM: usize, const YDIM: usize> {
    /// collision timescale parameter
    pub tau: f64,
    /// particle bound-back applied to boundaries
    pub pbb: bool,
    /// density at each lattice
    pub rho: [[f64; YDIM]; XDIM],
    /// velocity at each lattice
    pub u: [[(f64, f64); YDIM]; XDIM],
    /// boundaries
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
            pbb: true,
        };

        // initialize the values of the lattice field
        // using a given base average density
        for i in 0..X {
            for j in 0..Y {
                let rho: f64 = lbm.f[i][j].iter().sum();
                for k in 0..NDIR {
                    lbm.f[i][j][k] *= rho_0 / rho;
                }

                lbm.f[i][j][1] *= 5.0;
            }
        }

        lbm.calc_macro();

        lbm
    }

    #[inline]
    fn create_field<T: Copy>(fill: T) -> [[T; Y]; X] {
        [[fill; Y]; X]
    }

    pub fn sim(&mut self, steps: usize) {
        for _ in 0..steps {
            self.clean_boundaries();
            self.calc_macro();
            self.collide();
            self.stream();

            if self.pbb {
                self.reflect_boundaries();
            }
        }
    }

    fn clean_boundaries(&mut self) {
        for i in 0..X {
            for j in 0..Y {
                if self.boundaries.iter().any(|b| b.contains(i, j)) {
                    for k in 0..NDIR {
                        self.f[i][j][k] = 0.0;
                    }

                    self.rho[i][j] = 0.0;
                    self.u[i][j] = (0.0, 0.0);
                }
            }
        }
    }

    /// stream step of the Lattice Boltzmann Method,
    /// where the distributions of the outer surrounding lattice sites
    /// are propagated to neighboring latticies
    fn stream(&mut self) {
        let mut f_new = Self::create_field([0.0; NDIR]);

        for k in 0..NDIR {
            let (vx, vy) = SITE_VECS[k];

            for i in 0..X {
                for j in 0..Y {
                    let i_new = (i as f64 + vx).rem_euclid(X as _) as usize;
                    let j_new = (j as f64 + vy).rem_euclid(Y as _) as usize;

                    if self.boundaries.iter().any(|b| b.contains(i_new, j_new)) {
                        f_new[i][j][k] = self.f[i][j][k];
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

                    // uodate the lattice site using the LBM step
                    self.f[i][j][k] -= (self.f[i][j][k] - f_eq_ijk) / self.tau;
                }
            }
        }
    }

    /// recompute the values of `ρ` (density) and `u` (velocity)
    /// based on the current state of the lattice field
    fn calc_macro(&mut self) {
        for i in 0..X {
            for j in 0..Y {
                // density is the sum of the particle distribution within the lattice
                let rho = self.f[i][j].iter().sum();

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
                if rho > 0.0 {
                    self.u[i][j] = (x / rho, y / rho);
                }
            }
        }
    }

    fn reflect_boundaries(&mut self) {
        for boundary in &self.boundaries {
            for i in 0..X {
                for j in 0..Y {
                    if boundary.contains(i, j) {
                        // reflect this point
                    }
                }
            }
        }
    }
}

/// weights of sites in the D2Q9 lattice
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

/// coordinates describing positions of sites in the D2Q9 lattice
const SITE_VECS: [(f64, f64); 9] = [
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.0),
    (0.0, 1.0),
    (-1.0, 1.0),
    (-1.0, 0.0),
    (-1.0, -1.0),
    (0.0, -1.0),
    (1.0, -1.0),
];

/// reversed lattice directions used for reflective boundaries
const SITE_REV: [usize; 9] = [0, 5, 6, 7, 8, 1, 2, 3, 4];

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_all() {
        const XDIM: usize = 10;
        const YDIM: usize = 10;

        // This is an example of the workflow for a simple single-phase flow
        // with an obstacle.
        let mut s = LBM::<XDIM, YDIM>::new(100.0);

        struct Circle<const X: isize>;
        impl<const RADIUS: isize> Boundary for Circle<RADIUS> {
            fn contains(&self, x: usize, y: usize) -> bool {
                let dx = XDIM as isize / 2 - x as isize;
                let dy = YDIM as isize / 2 - y as isize;
                dx * dx + dy * dy < RADIUS.pow(2)
            }
        }

        s.boundaries.push(&Circle::<3>);

        s.sim(10);

        for (i, row) in s.rho.iter().enumerate() {
            for (j, _) in row.iter().enumerate() {
                print!("{:.2?} ", s.rho[i][j]);
            }
            println!("");
        }
    }
}
