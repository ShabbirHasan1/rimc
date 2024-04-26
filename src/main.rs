use clap::Parser;
use ndarray::Array2;
use rand::prelude::*;
use std::time::Instant;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Dimension for the NxN ensemble
    #[clap(short, long, default_value_t = 10)]
    dim: usize,
}

#[derive(Debug, Clone, Copy)]
struct Params {
    mag_field: f32,
    coupling_const: f32,
    beta: f32,
}

#[derive(Debug)]
struct Ensemble {
    arr: Array2<i8>,
    params: Params,
    rng: ThreadRng,
}

trait System {
    fn new(dim: usize, params: Params) -> Self;
    fn init_ensemble(&mut self);
    fn calc_site_energy(&self, i: usize, j: usize) -> Option<f32>;
    fn calc_site_energy_diff(&self, i: usize, j: usize) -> Option<f32>;
    fn calc_state_energy(&self) -> Option<f32>;
    fn step(&mut self) -> Result<(), &'static str>;
    fn update_site(&mut self, i: usize, j: usize) -> Result<(), &'static str>;
}

impl System for Ensemble {
    fn new(dim: usize, params: Params) -> Self {
        Ensemble {
            arr: Array2::<i8>::zeros((dim, dim)),
            params,
            rng: rand::thread_rng(),
        }
    }

    /// Initials the ensemble with a normal distribution of [-1, 1] states.
    fn init_ensemble(&mut self) {
        let options = [-1, 1];

        for el in self.arr.iter_mut() {
            *el = *options.choose(&mut self.rng).unwrap();
        }
    }

    /// Calculates a given sites energy for the Ising model.
    fn calc_site_energy(&self, i: usize, j: usize) -> Option<f32> {
        let coupling_const = &self.params.coupling_const;

        let s_k = self.arr[[i, j]] as f32;
        let mut site_energy: f32 = 0f32;

        let nrows = self.arr.shape()[0] as i128;
        let ncols = self.arr.shape()[1] as i128;

        let offsets = [-1, 1];
        let mut x: usize;
        let mut y: usize;

        for &n in &offsets {
            x = (i as i128 + n).rem_euclid(nrows) as usize;
            y = (j as i128 + n).rem_euclid(ncols) as usize;

            site_energy += self.arr[[x, j]] as f32;
            site_energy += self.arr[[i, y]] as f32;
        }
        site_energy *= -coupling_const * s_k;

        Some(site_energy)
    }

    /// Calculates the site energy difference if thie site at [i, j] were to flip.
    fn calc_site_energy_diff(&self, i: usize, j: usize) -> Option<f32> {
        let mut site_energy = self.calc_site_energy(i, j).unwrap();
        site_energy *= -2.;

        Some(site_energy)
    }

    /// Calculate the current states energy for all sites.
    fn calc_state_energy(&self) -> Option<f32> {
        let nrows = self.arr.shape()[0];
        let ncols = self.arr.shape()[1];

        let mut state_energy: f32 = 0f32;

        for i in 0..nrows {
            for j in 0..ncols {
                state_energy += self
                    .calc_site_energy(i, j)
                    .expect("Expected value, but got None");
            }
        }
        Some(state_energy)
    }

    /// Iterates through one step and decides to flip a spin.
    fn step(&mut self) -> Result<(), &'static str> {
        //const COORD_NUM: f32 = 4.0;

        //let J = self.params.coupling_const;
        let beta = self.params.beta;

        let nrows = self.arr.shape()[0];
        let ncols = self.arr.shape()[1];

        let i: usize = self.rng.gen_range(0..nrows);
        let j: usize = self.rng.gen_range(0..ncols);
        let r: f32 = self.rng.gen();

        let accept_ratio: f32;
        let mut exponent = self.calc_site_energy_diff(i, j).unwrap();
        match exponent {
            x if x <= 0.0 => accept_ratio = 1.0,
            x if x > 0.0 => {
                exponent *= -beta;
                accept_ratio = exponent.exp();
            }
            _ => return Err("Invalid number to calculate Acceptance Ratio"),
        }

        if r < accept_ratio {
            self.update_site(i, j).unwrap();
        }

        Ok(())
    }

    /// Toggle the current state from -1 to 1 or from 1 to -1
    fn update_site(&mut self, i: usize, j: usize) -> Result<(), &'static str> {
        let spin = &self.arr[[i, j]];

        match spin {
            -1 => self.arr[[i, j]] = 1,
            1 => self.arr[[i, j]] = -1,
            _ => return Err("Invalid state to flip"),
        }

        Ok(())
    }
}

fn write_state(ensemble: &Ensemble, filename: &str) -> std::io::Result<()> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    for row in ensemble.arr.rows() {
        let line = row
            .iter()
            .map(|&x| x.to_string())
            .collect::<Vec<String>>()
            .join(" ");
        writeln!(writer, "{}", line)?;
    }

    Ok(())
}

fn main() {
    let cli = Cli::parse();

    let dim = cli.dim;

    let my_params = Params {
        beta: 1.0,
        mag_field: 0.0,
        coupling_const: 1.0,
    };

    let mut my_ensemble = Ensemble::new(dim, my_params.clone());

    println!("Ensemble {} x {}\n{:#?}", dim, dim, my_ensemble);

    my_ensemble.init_ensemble();
    println!("Init\n{:#?}", my_ensemble.arr);
    write_state(&my_ensemble, "init.txt").unwrap();

    let t_start = Instant::now();
    for _ in 0..my_ensemble.arr.len() {
        my_ensemble.step().unwrap();
    }
    let duration = t_start.elapsed();
    println!("Finished\n{:#?}", my_ensemble.arr);
    println!(
        "Duration: {:?} for {} steps",
        duration,
        my_ensemble.arr.len()
    );

    write_state(&my_ensemble, "final.txt").unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn setup() -> Ensemble {
        let params = Params {
            coupling_const: 1.0,
            mag_field: 0.0,
            beta: 1.0,
        };
        let mut ensemble = Ensemble::new(4, params.clone());

        // Force system state, i.e., set up array with alternating -1 and 1
        for i in 0..4 {
            for j in 0..4 {
                if (i + j) % 2 == 0 {
                    ensemble.arr[[i, j]] = 1;
                } else {
                    ensemble.arr[[i, j]] = -1;
                }
            }
        }

        ensemble
    }

    #[test]
    fn test_site_energy() {
        // Sum over +1 spin with 4 neighbouring spins that are -1
        let ensemble = setup();
        let s_k = ensemble.arr[[1, 1]] as f32;
        let value = -ensemble.params.coupling_const * s_k * -4.0;

        assert_relative_eq!(
            ensemble.calc_site_energy(1, 1).unwrap(),
            value,
            epsilon = 1E-6
        );
    }

    #[test]
    fn test_site_energy_diff() {
        let ensemble = setup();
        let s_k = ensemble.arr[[1, 1]] as f32;
        let value = -ensemble.params.coupling_const * s_k * -4.0 * -2.0;
        assert_relative_eq!(
            ensemble.calc_site_energy_diff(1, 1).unwrap(),
            value,
            epsilon = 1E-6
        );
    }

    #[test]
    fn test_state_energy() {
        let ensemble = setup();
        // Each site is -J * s_k * sum_i=0..4 s_i, which is
        // -1.0 * -/+1.0 * +/-4.0 * 16.0 = 64
        let value: f32 = 64.0;
        assert_relative_eq!(ensemble.calc_state_energy().unwrap(), value, epsilon = 1E-6);
    }
}
