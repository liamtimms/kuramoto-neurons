extern crate ndarray;

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
    #[clap(short, long, value_parser, default_value_t = 2000)]
    timesim: usize,

    /// Size of time step using runge-kutta
    #[clap(long, value_parser, default_value_t = 0.01)]
    dt: f64,

    /// Time metric. Zannette's "T"
    #[clap(long, value_parser, default_value_t = 2.0)]
    timemetric: f64,

    /// Standard deviation in natural frequency distribution
    #[clap(short, long, value_parser, default_value_t = 0.25)]
    spreadinomega: f64,

    /// Initial Coupling strength
    #[clap(short, long, value_parser, default_value_t = 1.0)]
    g: f64,

    /// Speed of coupling evolution
    #[clap(short, long, value_parser, default_value_t = 0.1)]
    epsilon: f64,

    /// driving frequency
    #[clap(short, long, value_parser, default_value_t = 1.0)]
    drivingfrequency: f64,

    /// Size of cluster
    #[clap(short, long, value_parser, default_value_t = 0)]
    clustersize: usize,
}

fn main() {
    // all parameters
    let cli = Args::parse();

    let dimension = cli.dimension;
    if dimension != 1 && dimension != 2 {
        panic!("Dimension must be 1 or 2");
    }

    let output_dir: PathBuf = PathBuf::from(cli.output_dir);
    if !output_dir.exists() {
        panic!("Output directory does not exist");
    }

    let n: usize = cli.number; //Either the number of oscillators (1D) or the height/width of the square array
    let timesim: usize = cli.timesim; //the number of total time steps we are simulating (counts the initial conditions)

    // let dimension = 2;
    // let n: usize = 20; //Either the number of oscillators (1D) or the height/width of the square array
    // let tmax: usize = 2000; //the number of total time steps we are simulating (counts the initial conditions)

    let dt = cli.dt; //the time step we are using

    let spreadinomega: f64 = cli.spreadinomega; //the standard deviation in the natural frequencies
                                                //
    let g: f64 = cli.g; //The initial coupling constant
    let epsilon: f64 = cli.epsilon; //Speed of coupling change
    let timemetric: f64 = cli.timemetric; //Zannette's "T" used only for initial conditions if timemetric2 != it.

    // let tinitial = 200; // the time step we start the simulation at
    // tinitial must be less than tmax and longer than longest time delay in the system so we
    // calculate it in the code now

    let _mode: i32 = 0; //Sets initial conditions:
                        // 0 -> random phases, ONLY ONE CURRENTLY IMPLEMENTED
                        // 1 -> same phase cluster,
                        // 2 -> plane wave cluster,
                        // 3 -> circles cluster
    let clustersize: usize = cli.clustersize; //cluster=N for whole array to be driven
    if clustersize > n {
        panic!("Cluster size must be less than or equal to the number of oscillators");
    }
    let drivingfrequency: f64 = cli.drivingfrequency; //frequency that the oscillaotrs in the cluster will be forced to move at

    // setting up the output file name for later
    let mut phi_save_name = "phi_".to_string();
    write!(phi_save_name, "N{}-tmax{}-D{}.npy", n, timesim, dimension).unwrap();
    let phi_save_name = output_dir.join(phi_save_name);

    if dimension == 1 {
        // calling the 1D simulation
        let phi = onedimensional::run(
            n,
            timesim,
            dt,
            spreadinomega,
            g,
            epsilon,
            timemetric,
            clustersize,
            drivingfrequency,
        );

        // save data with error handling
        match write_npy(&phi_save_name, &phi) {
            Ok(_) => println!("Saved to {}", phi_save_name.display()),
            Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
        }
    } else if dimension == 2 {
        // calling the 2D simulation
        let phi = twodimensional::run(
            n,
            timesim,
            dt,
            spreadinomega,
            g,
            epsilon,
            timemetric,
            clustersize,
            drivingfrequency,
        );
        match write_npy(&phi_save_name, &phi) {
            Ok(_) => println!("Saved to {}", phi_save_name.display()),
            Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
        }
    } else {
        // we haven't implemented anything for higher dimensions yet
        println!("Dimension {} not implemented", dimension);
    }
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
