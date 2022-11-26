use indicatif::ProgressBar;
use ndarray::parallel::prelude::*;
use ndarray::prelude::*;
use ndarray::Array;
use ndarray_rand::rand::distributions::{Distribution, Uniform};
use ndarray_rand::rand_distr;
use ndarray_rand::RandomExt;
use ndarray_stats::QuantileExt;

use crate::utils;

/// Two dimensional model
struct TwoDimensional {
    phi: Array<f64, Ix3>,
    kcoupling: Array<f64, Ix5>,
    omega: Array<f64, Ix2>,
    tau: Array<usize, Ix4>,
    connectionmatrix: Array<f64, Ix4>,
    alpha: Array<f64, Ix4>,
}

impl TwoDimensional {
    fn tmax(&self) -> usize {
        self.phi.shape()[2]
    }
}

fn get_tmax(phi: &Array<f64, Ix3>) -> usize {
    // set the maximum time step
    phi.shape()[2]
}

fn find_delay(i: &usize, j: &usize, k: &usize, l: &usize, n: &usize, z: &f64) -> usize {
    // find the delay between oscillator at i, j and it's neighbor at k, l
    // tau[i][j][k][l]=  trunc(Z/N*sqrt( (min(abs(i-k), N-abs(i-k)))^2 +(min(abs(j-l), N-abs(j-l)))^2 ))
    let x = (*i as f64 - *k as f64).abs();
    let x_prime = *n as f64 - x; // wrap around
    let y = (*j as f64 - *l as f64).abs();
    let y_prime = *n as f64 - y; // wrap around
    let d = (utils::min(x, x_prime).powi(2) + utils::min(y, y_prime).powi(2)).sqrt();

    ((z / *n as f64) * d).round() as usize
}

fn initialize_phi(n: &usize, clustersize: &usize, tmax: &usize) -> Array<f64, Ix3> {
    // replaces "function initialconditions1(N,clustersize, dimension)" in igor
    let mut phi: Array<f64, Ix3> = Array::zeros((*n, *n, *tmax));

    let mut rng = ndarray_rand::rand::thread_rng();
    let circle = Uniform::new(0.0, 2.0 * std::f64::consts::PI);

    for i in 0..*clustersize {
        for j in 0..*clustersize {
            phi[[i, j, 0]] = std::f64::consts::PI;
        }
    }
    for i in *clustersize..*n {
        for j in *clustersize..*n {
            phi[[i, j, 0]] = circle.sample(&mut rng);
        }
    }

    phi
}

fn calculate_delays(timemetric: &f64, n: &usize, dt: &f64) -> Array<usize, Ix4> {
    // calculate the delays for each oscillator
    let z: f64 = timemetric / dt;
    let mut tau: Array<usize, Ix4> = Array::zeros((*n, *n, *n, *n));

    for i in 0..*n {
        for j in 0..*n {
            for k in 0..*n {
                for l in 0..*n {
                    if i != j {
                        tau[[i, j, k, l]] = find_delay(&i, &j, &k, &l, &*n, &z);
                    } else {
                        tau[(i, j, k, l)] = 0;
                    }
                }
            }
        }
    }
    tau
}

fn set_tinitial(tau: &Array<usize, Ix4>) -> usize {
    // set the initial time step
    match tau.argmax() {
        Ok(max_tau) => {
            let max_tau = tau[max_tau];
            max_tau as usize * 3 // fill out thrice the max delay
        }
        Err(_) => {
            println!("Error: no maximum tau found");
            0 // in principle this should never happen
        }
    }
}

fn driver(
    mut phi: Array<f64, Ix3>,
    dt: &f64,
    drivingfrequency: &f64,
    n: &usize,
    clustersize: &usize,
    tinitial: &usize,
    omega: &Array<f64, Ix2>,
) -> Array<f64, Ix3> {
    for t in 0..*tinitial {
        for i in 0..*n {
            for j in 0..*n {
                let timewithstep = t + 1;
                if i < *clustersize {
                    // phi[[i, timewithstep]] = phi[[i, t]] + drivingfrequency * dt;
                    phi[[i, j, timewithstep]] = phi[[i, j, t]] + drivingfrequency * dt;
                } else if i >= *clustersize {
                    // 0<t<tinitial phase assignments
                    phi[[i, j, timewithstep]] = phi[[i, j, t]] + omega[[i, j]] * dt;
                }
            }
        }
    }
    phi
}

fn update_oscillator(
    t: &usize,
    i: &usize,
    j: &usize,
    phi: &Array<f64, Ix3>,
    tau: &Array<usize, Ix4>,
    omega: &Array<f64, Ix2>,
    kcoupling: &Array<f64, Ix5>,
    connectionmatrix: &Array<f64, Ix4>,
    dt: &f64,
) -> f64 {
    let mut summation = 0.0;
    let phi_current = phi[[*i, *j, *t]];
    let n = phi.shape()[0];
    for k in 0..n {
        for l in 0..n {
            let timewithdelay = t - tau[[*i, *j, k, l]] as usize;
            let diff = phi[[k, l, timewithdelay]] - phi_current;
            summation +=
                kcoupling[[*i, *j, k, l, *t]] * diff.sin() * connectionmatrix[[*i, *j, k, l]];
        }
    }
    // returning a "phi_new"
    phi_current + (omega[[*i, *j]] + (summation / n as f64)) * dt
}

fn update_oscillator_parallel(
    t: &usize,
    i: &usize,
    j: &usize,
    phi: &Array<f64, Ix3>,
    tau: &Array<usize, Ix4>,
    omega: &Array<f64, Ix2>,
    kcoupling: &Array<f64, Ix5>,
    connectionmatrix: &Array<f64, Ix4>,
    dt: &f64,
) -> f64 {
    // SLOWER than update_oscillator in all tests
    let phi_current = phi[[*i, *j, *t]];
    let n = phi.shape()[0];
    // based on https://www.anycodings.com/1questions/5406370/how-to-use-rayon-for-parallel-calculation-of-pi
    let summation: f64 = (0..n)
        .into_par_iter()
        .map(|k| {
            let mut inner_sum = 0.0;
            for l in 0..n {
                let timewithdelay = t - tau[[*i, *j, k, l]] as usize;
                let diff = phi[[k, l, timewithdelay]] - phi_current;
                inner_sum +=
                    kcoupling[[*i, *j, k, l, *t]] * diff.sin() * connectionmatrix[[*i, *j, k, l]];
            }
            inner_sum
        })
        .reduce(|| 0.0, |a, b| a + b);

    // returning a "phi_new"
    phi_current + (omega[[*i, *j]] + (summation / n as f64)) * dt
}

fn model(
    // parameters for oscillators
    mut phi: Array<f64, Ix3>,
    mut kcoupling: Array<f64, Ix5>,
    omega: &Array<f64, Ix2>,
    // parameters of the simulation itself
    dt: &f64,
    tinitial: &usize,
    n: &usize,
    epsilon: &f64,
    // parameters for connections
    tau: &Array<usize, Ix4>,
    connectionmatrix: &Array<f64, Ix4>,
    alpha: &Array<f64, Ix4>,
) -> (Array<f64, Ix3>, Array<f64, Ix5>) {
    // replaces "function model(dt, tinitial, tmax, N, epsilon, dimension)" in igor
    // let mut phicurrent: Array<f64, Ix2> = Array::zeros((*n, *n));

    let tmax = get_tmax(&phi);
    let pb = ProgressBar::new(((tmax - 1) - *tinitial) as u64);

    // BEGINING OF TIME (2D)
    println!("starting model");
    for t in *tinitial..(tmax - 1) {
        pb.inc(1); // update progressbar

        // Evolution of the Oscillators
        if *n < 100 {
            for i in 0..*n {
                for j in 0..*n {
                    phi[[i, j, t + 1]] = update_oscillator(
                        &t,
                        &i,
                        &j,
                        &phi,
                        &tau,
                        &omega,
                        &kcoupling,
                        &connectionmatrix,
                        dt,
                    );
                }
            }
        } else {
            for i in 0..*n {
                for j in 0..*n {
                    phi[[i, j, t + 1]] = update_oscillator_parallel(
                        &t,
                        &i,
                        &j,
                        &phi,
                        &tau,
                        &omega,
                        &kcoupling,
                        &connectionmatrix,
                        dt,
                    );
                }
            }
        }

        //Evolution of Coupling
        for i in 0..*n {
            for j in 0..*n {
                for k in 0..*n {
                    for l in 0..*n {
                        let timewithdelay = t - tau[[i, j, k, l]] as usize;
                        let factor = phi[[i, j, t]] - phi[[k, l, timewithdelay]];
                        kcoupling[[i, j, k, l, t + 1]] = kcoupling[[i, j, k, l, t]]
                            + (epsilon
                                * (alpha[[i, j, k, l]] * factor.cos()
                                    - kcoupling[[i, j, k, l, t]])
                                * dt);
                    }
                }
            }
        }
    }
    // END OF TIME (1D)
    pb.finish_with_message("Done with modeling");
    (phi, kcoupling)
}

fn setup(
    n: usize,
    spreadinomega: f64,
    timemetric: &f64,
    dt: &f64,
    timesim: &usize,
    g: f64,
    clustersize: &usize,
) -> (
    TwoDimensional, // parameters for oscillators
    usize,
) {
    // contains immutable parameters describing the initial state of the system
    let alpha: Array<f64, _> = Array::ones((n, n, n, n)); // coupling speed eviolution
    let connectionmatrix: Array<f64, _> = Array::ones((n, n, n, n)); // whether coupled or not

    let omega: Array<f64, Ix2> =
        // Array::random(n, rand_distr::Normal::new(1.0, spreadinomega).unwrap());
    Array::random((n, n), rand_distr::Normal::new(1.0, spreadinomega).unwrap());

    let tau: Array<usize, Ix4> = calculate_delays(&timemetric, &n, &dt);

    // initial conditions
    let tinitial: usize = set_tinitial(&tau);
    println!("tinitial = {}", tinitial);
    let tmax: usize = tinitial + timesim;

    let phi: Array<f64, Ix3> = initialize_phi(&n, &clustersize, &tmax);

    let mut kcoupling: Array<f64, Ix5> = Array::zeros((n, n, n, n, tmax)); // coupling strength
    kcoupling.fill(g);

    let two_dimensional_plane = TwoDimensional {
        phi,
        kcoupling,
        omega,
        tau,
        connectionmatrix,
        alpha,
    };

    (two_dimensional_plane, tinitial)
}

pub fn run(run_params: &utils::RunParams) -> (Array<f64, Ix3>, Array<f64, Ix5>) {
    let (mut two_dimensional_plane, tinitial) = setup(
        run_params.n,
        run_params.spreadinomega,
        &run_params.timemetric,
        &run_params.dt,
        &run_params.timesim,
        run_params.g,
        &run_params.clustersize,
    );

    // driving the cluster for the initial time
    two_dimensional_plane.phi = driver(
        two_dimensional_plane.phi,
        &run_params.dt,
        &run_params.drivingfrequency,
        &run_params.n,
        &run_params.clustersize,
        &tinitial,
        &two_dimensional_plane.omega,
    );

    // run the actual model
    (two_dimensional_plane.phi, two_dimensional_plane.kcoupling) = model(
        two_dimensional_plane.phi,
        two_dimensional_plane.kcoupling,
        &two_dimensional_plane.omega,
        &run_params.dt,
        &tinitial,
        &run_params.n,
        &run_params.epsilon,
        &two_dimensional_plane.tau,
        &two_dimensional_plane.connectionmatrix,
        &two_dimensional_plane.alpha,
    );

    (two_dimensional_plane.phi, two_dimensional_plane.kcoupling)
}
