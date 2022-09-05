extern crate ndarray;

use core::panic;
use ndarray::Array;
use ndarray_npy::write_npy;
use std::fmt::Write;
use std::path::PathBuf;

use clap::Parser;

mod onedimensional;
mod twodimensional;
mod utils;

/// CLI arguments
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// The directory to save the results to
    #[clap(short, long, value_parser, default_value_t = String::from("."))]
    output_dir: String,

    /// Number of dimensions of the system
    #[clap(short, long, value_parser, default_value_t = 1)]
    dimension: i32,

    /// Number of oscillators
    #[clap(short, long, value_parser, default_value_t = 20)]
    number: usize,

    /// Number of steps to run the simulation
    #[clap(short, long, value_parser, default_value_t = 500)]
    timesim: usize,
}

fn main() {
    // all parameters
    let cli = Args::parse();

    let dimension = cli.dimension;
    if dimension != 1 && dimension != 2 {
        panic!("Dimension must be 1 or 2, but was {dimension}");
    }

    let output_dir: PathBuf = PathBuf::from(cli.output_dir);
    if !output_dir.exists() {
        panic!("Output directory not found. Create it first.");
    }

    let n: usize = cli.number; //Either the number of oscillators (1D) or the height/width of the square array
    let timesim: usize = cli.timesim; //the number of total time steps to run the simulation for

    let dt = 0.01; //the time step we are using
                   //
    let spreadinomega: f64 = 0.10; //the standard deviation in the natural frequencies
                                   //
    let g: f64 = 1.0; //The initial coupling constant
    let epsilon: f64 = 0.1; //Speed of coupling change
    let timemetric: f64 = 2.0; //Zannette's "T" used only for initial conditions if timemetric2 != it.

    // let tinitial = 200; // the time step we start the simulation at
    // tinitial must be less than tmax and longer than longest time delay in the system so we
    // calculate it in the code now

    let _mode: i32 = 0; //Sets initial conditions:
                        // 0 -> random phases, ONLY ONE CURRENTLY IMPLEMENTED
                        // 1 -> same phase cluster,
                        // 2 -> plane wave cluster,
                        // 3 -> circles cluster
    let clustersize = n / 2 as usize; //cluster=N for whole array to be driven
    let drivingfrequency: f64 = 1.0; //frequency that the oscillaotrs in the cluster will be forced to move at

    // setting up the output file name for later
    let mut phi_save_name = "phi_".to_string();
    write!(phi_save_name, "N{}-time{}-D{}.npy", n, timesim, dimension).unwrap();
    let phi_save_name = output_dir.join(phi_save_name);

    let phi = match dimension {
        1 => onedimensional::run(
            n,
            timesim,
            dt,
            spreadinomega,
            g,
            epsilon,
            timemetric,
            clustersize,
            drivingfrequency,
        ),
        2 => twodimensional::run(
            n,
            timesim,
            dt,
            spreadinomega,
            g,
            epsilon,
            timemetric,
            clustersize,
            drivingfrequency,
        ),
        _ => panic!("Dimension must be 1 or 2, but was {}", dimension),
    };

    match write_npy(&phi_save_name, &phi) {
        Ok(_) => println!("Saved to {}", phi_save_name.display()),
        Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
    }

    // if dimension == 1 {
    //     // calling the 1D simulation
    //     let phi = onedimensional::run(
    //         n,
    //         timesim,
    //         dt,
    //         spreadinomega,
    //         g,
    //         epsilon,
    //         timemetric,
    //         clustersize,
    //         drivingfrequency,
    //     );
    //
    //     // save data with error handling
    //     match write_npy(&phi_save_name, &phi) {
    //         Ok(_) => println!("Saved to {}", phi_save_name.display()),
    //         Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
    //     }
    // } else if dimension == 2 {
    //     // calling the 2D simulation
    //     let phi = twodimensional::run(
    //         n,
    //         timesim,
    //         dt,
    //         spreadinomega,
    //         g,
    //         epsilon,
    //         timemetric,
    //         clustersize,
    //         drivingfrequency,
    //     );
    //     match write_npy(&phi_save_name, &phi) {
    //         Ok(_) => println!("Saved to {}", phi_save_name.display()),
    //         Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
    //     }
    // } else {
    //     // we haven't implemented anything for higher dimensions yet
    //     println!("Dimension {} not implemented", dimension);
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_finds_max() {
        /// my first test!
        assert_eq!(4.0, utils::max(&3.5, &4.0));
    }
}
