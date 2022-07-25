#[macro_use]
extern crate ndarray;

use ndarray::prelude::*;

use std::cmp;

fn kuramoto_controller() {}

fn calculate_delays(
    timemetric: &f64,
    N: &usize,
    dt: &f64,
    dimension: &i32,
    tau: &mut Array<f64, Ix2>,
) {
    let Z = timemetric / dt;
    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut d: f64 = 0.0;

    if *dimension == 1 {
        for i in 0..*N {
            for j in 0..*N {
                // tau[i, j]=trunc(Z/N*min(abs(i-j),N-abs(i-j)));
                x = (i as f64 - j as f64).abs();
                y = (*N as f64 - x).abs();
                if x < y {
                    d = x;
                } else {
                    d = y;
                }
                tau[[i, j]] = Z / d;
                // tau[[i, j]] = Z / N as f64 * cmp::min(x, N as f64 - x);
            }
        }
    }
}

fn main() {
    // all parameters
    let dimension = 1;
    let g = 1.0; //The initial coupling constant
    let N: usize = 20; //Either the number of oscillators (1D) or the height/width of the square array
    let epsilon = 0.1;
    let dt = 0.01;
    let timemetric: f64 = 2.0; //Zannette's "T" used only for initial conditions if timemetric2 != it.
    let tinitial = 1000;

    let mode = 0; //Sets initial conditions: 0->random phases, 1->same phase cluster, 2->plane wave cluster, 3->circles cluster
    let clustersize = N; //cluster=N for whole array to be driven
    let drivingfrequency = 1; //frequency that the oscillaotrs in the cluster will be forced to move at
    let spreadinomega = 0; //the standard deviation in the natural frequencies
    let tmax = 40000; //the number of total time steps we are simulating (counts the initial conditions)

    if dimension == 1 {
        let epsilon_counter = 0;

        let mut phi: Array<f64, Ix2> = Array::zeros((N, tmax)); // 2D array, must be big enough to hold largest delay and total run time
        let mut phicurrent: Array<f64, _> = Array::zeros(N);
        let mut omega: Array<f64, _> = Array::zeros(N);
        let mut Kcoupling: Array<f64, Ix2> = Array::zeros((N, N));
        let mut tau: Array<f64, Ix2> = Array::zeros((N, N));
        let mut FinalFrequencies: Array<f64, _> = Array::zeros(N);
        let mut OrderParameter: Array<f64, _> = Array::zeros((4, tmax));
        let mut OrderParameter2Cluster: Array<f64, _> = Array::zeros((4, tmax));
        let mut OrderParameter2ClusterConnected: Array<f64, _> = Array::zeros((4, tmax));
        let mut tt: Array<f64, _> = Array::zeros(N * tmax);
        let mut freqs: Array<f64, _> = Array::zeros(N * tmax);

        let mut alpha: Array<f64, _> = Array::ones((N, N));
        let mut connectionmatrix: Array<f64, _> = Array::ones((N, N));
        Kcoupling.fill(g);
        calculate_delays(&timemetric, &N, &dt, &dimension, &mut tau);
    }

    println!("Hello, world!");
}
