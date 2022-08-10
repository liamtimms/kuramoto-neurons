extern crate ndarray;

use ndarray_npy::write_npy;
use std::fmt::Write;

use clap::Parser;
// use ndarray::Zip;
// use rayon::prelude::*;

mod onedimensional;
mod twodimensional;
mod utils;

/// CLI arguments
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Number of dimensions of the system
    #[clap(short, long, value_parser, default_value_t = 2)]
    dimension: i32,

    /// Number of oscillators
    #[clap(short, long, value_parser, default_value_t = 20)]
    number: usize,

    /// Number of steps to run the simulation
    #[clap(short, long, value_parser, default_value_t = 2000)]
    tmax: usize,
}

fn main() {
    // all parameters
    let cli = Args::parse();

    let dimension = cli.dimension;
    let n: usize = cli.number; //Either the number of oscillators (1D) or the height/width of the square array
    let tmax: usize = cli.tmax; //the number of total time steps we are simulating (counts the initial conditions)

    // let dimension = 2;
    // let n: usize = 20; //Either the number of oscillators (1D) or the height/width of the square array
    // let tmax: usize = 2000; //the number of total time steps we are simulating (counts the initial conditions)

    let dt = 0.01; //the time step we are using
                   //
    let spreadinomega: f64 = 0.25; //the standard deviation in the natural frequencies
                                   //
    let g: f64 = 1.0; //The initial coupling constant
    let epsilon: f64 = 0.1; //Speed of coupling change
    let timemetric: f64 = 2.0; //Zannette's "T" used only for initial conditions if timemetric2 != it.
    let tinitial = 200; // the time step we start the simulation at
                        // tinitial must be less than tmax and longer than longest time delay in the system

    let _mode: i32 = 0; //Sets initial conditions:
                        // 0->random phases, only one currently implemented
                        // 1->same phase cluster, 2->plane wave cluster, 3->circles cluster
    let clustersize = 0; //cluster=N for whole array to be driven
    let drivingfrequency: f64 = 1.0; //frequency that the oscillaotrs in the cluster will be forced to move at

    if dimension == 1 {
        let phi = onedimensional::run(
            n,
            tmax,
            dt,
            spreadinomega,
            g,
            epsilon,
            timemetric,
            tinitial,
            clustersize,
            drivingfrequency,
        );

        let mut phi_save_name = "phi_".to_string();
        write!(phi_save_name, "N{}-tmax{}-dim{}.npy", n, tmax, dimension).unwrap();
        write_npy(&phi_save_name, &phi);
    } else if dimension == 2 {
        let phi = twodimensional::run(
            n,
            tmax,
            dt,
            spreadinomega,
            g,
            epsilon,
            timemetric,
            tinitial,
            clustersize,
            drivingfrequency,
        );
        let mut phi_save_name = "phi_".to_string();
        write!(phi_save_name, "N{}-tmax{}-dim{}.npy", n, tmax, dimension).unwrap();
        write_npy(&phi_save_name, &phi);
    } else {
        println!("Dimension {} not implemented", dimension);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_finds_max() {
        // my first test!
        assert_eq!(4.0, utils::max(&3.5, &4.0));
    }
}
