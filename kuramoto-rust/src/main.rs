//! Coupled Kuramoto oscillator simulation with learning and delay.
//!
//! Provides a way to simulate 1D and 2D Kuramoto oscillators with their coupling
//! governed by Hebbian learning rules and a time delay to simulate axons.
//! It is:
//! - Very fast for 1D case
//! - Faster than alternatives for 2D case
//! - Parallelized for large number of oscillators in 2D
//!
//! # Examples
//!
//! ```bash
//! $ mkdir simulation_data
//! $ ./kuramoto_rust -d 1 -n 50 -o simulation_data/
//! ```
//!
//! Multiple simulations can be run in parallel or series
//! using a tool like GNU parallel and some basic bash scripting.
//!
//! # Arguments
//!
//! Takes a set of command line arguments to specify the simulation parameters.
//! The basic parameters are:
//! - `-h` - help message, information about the parameters
//! - `-d` - dimension of the oscillator space
//! - `-n` - number of oscillators (1D case) or number of oscillators per dimension (2D case)
//! - `-t` - number of time steps to simulate
//! - `-o` - output directory
//!
//! Parameters to further specify the model and initial conditions:
//! - `-s` - spread (standard deviation) in natural fequencies (omega) of each oscillator
//! - `-e` - epsilon for Hebbian learning, i.e. learning rate
//! - `-c` - seperate cluster size
//! - `-g` - initial coupling strength
//! - `--timemetric` - Zannette's time metric
//! - `--drivingfrequency` - frequency to drive the oscillators in the cluster during initial
//! conditions
//!
//! # Output
//!
//! The output is a two `.npy` files, one for the phases and one for the coupling strengths.
//! These store the values of the phases and coupling strengths at each time step. These can then
//! be loaded and investigated further in Python with `numpy` and `matplotlib` (or other python plotting
//! libraries like `seaborn`.) See github README for more information.

#![warn(missing_docs)]
extern crate ndarray;

use ndarray_npy::write_npy;
use std::path::PathBuf;

use clap::Parser;

/// One dimensional model
mod onedimensional;
/// Two dimensional model
mod twodimensional;
/// Utility functions useful for all models
mod utils;

/// CLI arguments for the program
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The directory to save the results to
    #[arg(short, long, value_parser, default_value_t = String::from("."))]
    output_dir: String,

    /// Number of dimensions of the system
    #[arg(short, long, value_parser, default_value_t = 1)]
    dimension: i32,

    /// Number of oscillators
    #[arg(short, long, value_parser, default_value_t = 20)]
    number: usize,

    /// Number of steps to run the simulation
    #[arg(short, long, value_parser, default_value_t = 2000)]
    timesim: usize,

    /// Size of time step using runge-kutta
    #[arg(long, value_parser, default_value_t = 0.01)]
    dt: f64,

    /// Time metric. Zannette's "T"
    #[arg(long, value_parser, default_value_t = 2.0)]
    timemetric: f64,

    /// Standard deviation in natural frequency distribution
    #[arg(short, long, value_parser, default_value_t = 0.25)]
    spreadinomega: f64,

    /// Initial Coupling strength
    #[arg(long, value_parser, default_value_t = 1.0)]
    g: f64,

    /// Speed of coupling evolution
    #[arg(short, long, value_parser, default_value_t = 0.1)]
    epsilon: f64,

    /// driving frequency
    #[arg(long, value_parser, default_value_t = 1.0)]
    drivingfrequency: f64,

    /// Size of cluster
    #[arg(short, long, value_parser, default_value_t = 0)]
    clustersize: usize,
}

/// Main function parses arguments, calls the simulation modules and saves the results
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
    let drivingfrequency: f64 = cli.drivingfrequency; //frequency that the oscillaotrs in the cluster will be forced to move at during initial conditions

    // setting up the output file name for later
    let run_params = utils::RunParams {
        n,
        timesim,
        dt,
        spreadinomega,
        g,
        epsilon,
        timemetric,
        clustersize,
        drivingfrequency,
    };

    let identifier = format!(
        "N{n}-D{dimension}-t{timesim}-T{timemetric}-ostd{spreadinomega}-g{g}-e{epsilon}-C{clustersize}-do{drivingfrequency}.npy"
    );

    let phi_base = "phi_".to_string();
    let phi_save_name = utils::update_name(&output_dir, &phi_base, &identifier);

    let k_base = "K_".to_string();
    let k_save_name = utils::update_name(&output_dir, &k_base, &identifier);

    if dimension == 1 {
        // calling the 1D simulation
        let (phi, kcoupling) = onedimensional::run(&run_params);

        // save data with error handling
        match write_npy(&phi_save_name, &phi) {
            Ok(_) => println!("Saved to {}", phi_save_name.display()),
            Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
        }
        match write_npy(&k_save_name, &kcoupling) {
            Ok(_) => println!("Saved to {}", k_save_name.display()),
            Err(e) => println!("Error saving to {}: {}", k_save_name.display(), e),
        }
    } else if dimension == 2 {
        // calling the 2D simulation
        let (phi, kcoupling) = twodimensional::run(&run_params);

        match write_npy(&phi_save_name, &phi) {
            Ok(_) => println!("Saved to {}", phi_save_name.display()),
            Err(e) => println!("Error saving to {}: {}", phi_save_name.display(), e),
        }
        match write_npy(&k_save_name, &kcoupling) {
            Ok(_) => println!("Saved to {}", k_save_name.display()),
            Err(e) => println!("Error saving to {}: {}", k_save_name.display(), e),
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
        // my first test!
        assert_eq!(4.0, utils::max(&3.5, &4.0));
    }
}
