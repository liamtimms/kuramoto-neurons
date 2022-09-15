use indicatif::ProgressBar;
use ndarray::prelude::*;
use ndarray::Array;
use ndarray_rand::rand::distributions::{Distribution, Uniform};
use ndarray_rand::rand_distr;
use ndarray_rand::RandomExt;
use ndarray_stats::QuantileExt;

use crate::utils;

fn find_delay(i: &usize, j: &usize, n: &usize, z: &f64) -> usize {
    // find the delay between oscillator i and it's neighbor j
    let x = (*i as f64 - *j as f64).abs();
    let y = *n as f64 - x; // wrap around otherside of the circle
    let d = utils::min(x, y);
    ((z / *n as f64) * d).round() as usize
}

fn calculate_delays(timemetric: &f64, n: &usize, dt: &f64) -> Array<usize, Ix2> {
    // calculate the delays for each oscillator
    let z: f64 = timemetric / dt;
    let mut tau: Array<usize, Ix2> = Array::zeros((*n, *n));
    for i in 0..*n {
        for j in 0..*n {
            if i != j {
                tau[[i, j]] = find_delay(&i, &j, &*n, &z);
            } else {
                tau[(i, j)] = 0 as usize;
            }
        }
    }
    tau
}

fn set_tinitial(tau: &Array<usize, Ix2>) -> usize {
    // set the initial time step
    match tau.argmax() {
        Ok(max_tau) => {
            let max_tau = tau[max_tau];
            max_tau as usize * 3 // fill out thrice the delay
        }
        Err(_) => {
            println!("Error: no maximum tau found");
            0 // in principle this should never happen
        }
    }
}

fn calc_order_params(
    phi: &Array<f64, Ix2>,
    n: &usize,
    tmax: &usize,
) -> (Array<f64, Ix2>, Array<f64, Ix2>, Array<f64, Ix2>) {
    // calculate order parameters
    // This is reproduced from the original code
    // and CLEARLY NOT OPTIMIZED
    let mut sinsum = 0.0;
    let mut cossum = 0.0;

    let mut sinsum2 = 0.0;
    let mut cossum2 = 0.0;

    let mut sinsum3 = 0.0;
    let mut cossum3 = 0.0;

    let mut phicurrent: Array<f64, Ix1> = Array::zeros(*n);

    let mut order_parameter: Array<f64, Ix2> = Array::zeros((5, *tmax));
    let mut order_parameter2_cluster: Array<f64, Ix2> = Array::zeros((5, *tmax));
    let mut order_parameter2_cluster_connected: Array<f64, Ix2> = Array::zeros((5, *tmax));

    for m in 0..5 {
        for t in 0..*tmax {
            for i in 0..*n {
                phicurrent[i] =
                    phi[[i, t]] - (2.0 * std::f64::consts::PI * m as f64) * i as f64 / *n as f64;
                sinsum += (phicurrent[i]).sin();
                cossum += (phicurrent[i]).sin();
                sinsum2 += (2.0 * phicurrent[i]).sin();
                cossum2 += (2.0 * phicurrent[i]).cos();
                sinsum3 += (3.0 * phicurrent[i]).sin();
                cossum3 += (3.0 * phicurrent[i]).cos();
            }

            let op_minus = (sinsum.powi(2) + cossum.powi(2)).sqrt() / *n as f64;
            let op_minus2 = ((((sinsum2.powi(2) + cossum2.powi(2)).sqrt()).abs() / *n as f64
                - op_minus)
                .powi(2))
            .sqrt();
            let op_minus3 = ((((sinsum3.powi(2) + cossum3.powi(2)).sqrt()).abs() / *n as f64
                - op_minus)
                .powi(2))
            .sqrt();

            sinsum = 0.0;
            cossum = 0.0;
            sinsum2 = 0.0;
            cossum2 = 0.0;
            sinsum3 = 0.0;
            cossum3 = 0.0;

            for i in 0..*n {
                phicurrent[i] =
                    phi[[i, t]] + (2.0 * std::f64::consts::PI * m as f64) * i as f64 / *n as f64;
                sinsum += (phicurrent[i]).sin();
                cossum += (phicurrent[i]).sin();
                sinsum2 += (2.0 * phicurrent[i]).sin();
                cossum2 += (2.0 * phicurrent[i]).cos();
                sinsum3 += (3.0 * phicurrent[i]).sin();
                cossum3 += (3.0 * phicurrent[i]).cos();
            }

            let op_plus = (sinsum.powi(2) + cossum.powi(2)).sqrt() / *n as f64;
            let op_plus2 = ((((sinsum2.powi(2) + cossum2.powi(2)).sqrt()).abs() / *n as f64
                - op_plus)
                .powi(2))
            .sqrt();
            let op_plus3 = ((((sinsum3.powi(2) + cossum3.powi(2)).sqrt()).abs() / *n as f64
                - op_plus)
                .powi(2))
            .sqrt();

            sinsum = 0.0;
            cossum = 0.0;
            sinsum2 = 0.0;
            cossum2 = 0.0;
            sinsum3 = 0.0;
            cossum3 = 0.0;

            order_parameter[[m, t]] = utils::max(&op_plus, &op_minus);
            order_parameter2_cluster[[m, t]] = utils::max(&op_plus2, &op_minus2);
            order_parameter2_cluster_connected[[m, t]] = utils::max(&op_plus3, &op_minus3);
        }
    }
    (
        order_parameter,
        order_parameter2_cluster,
        order_parameter2_cluster_connected,
    )
}

fn initialize_phi(n: &usize, clustersize: &usize, tmax: &usize) -> Array<f64, Ix2> {
    // replaces "function initialconditions1(N,clustersize, dimension)" in igor
    let mut phi: Array<f64, Ix2> = Array::zeros((*n, *tmax));

    let mut rng = ndarray_rand::rand::thread_rng();
    let circle = Uniform::new(0.0, 2.0 * std::f64::consts::PI);

    for i in 0..*clustersize {
        phi[[i, 0]] = std::f64::consts::PI;
    }
    for i in *clustersize..*n {
        phi[[i, 0]] = circle.sample(&mut rng);
    }
    phi
}

fn driver(
    mut phi: Array<f64, Ix2>,
    dt: &f64,
    drivingfrequency: &f64,
    n: &usize,
    clustersize: &usize,
    tinitial: &usize,
    omega: &Array<f64, Ix1>,
) -> Array<f64, Ix2> {
    // replaces "function driver(dt, drivingfrequency, clustersize, dimension)" in igor
    for t in 0..*tinitial {
        for i in 0..*n {
            let timewithstep = t + 1;
            if i < *clustersize {
                phi[[i, timewithstep]] = phi[[i, t]] + drivingfrequency * dt;
            } else if i >= *clustersize {
                //0<t<tinitial phase assignments
                phi[[i, timewithstep]] = phi[[i, t]] + omega[i] * dt;
            }
        }
    }
    phi
}

fn update_oscillator(
    t: &usize,
    i: &usize,
    phi: &Array<f64, Ix2>,
    tau: &Array<usize, Ix2>,
    omega: &Array<f64, Ix1>,
    kcoupling: &Array<f64, Ix3>,
    connectionmatrix: &Array<f64, Ix2>,
    dt: &f64,
) -> f64 {
    let mut summation = 0.0;
    let phi_current = phi[[*i, *t]];
    let n = phi.shape()[0];
    for j in 0..n {
        let diff = phi[[j, (t - tau[[*i, j]])]] - phi_current;
        summation += kcoupling[[*i, j, *t]] * diff.sin() * connectionmatrix[[*i, j]];
    }
    let phi_new = phi_current + (omega[*i] + (summation / n as f64)) * dt;
    phi_new
}

// fn update_oscillator_parallel(
//     t: &usize,
//     i: &usize,
//     phi: &Array<f64, Ix2>,
//     tau: &Array<usize, Ix2>,
//     omega: &Array<f64, Ix1>,
//     kcoupling: &Array<f64, Ix2>,
//     connectionmatrix: &Array<f64, Ix2>,
//     dt: &f64,
// ) -> f64 {
//     // SLOWER than update_oscillator in all tests
//     let phi_current = phi[[*i, *t]];
//     let n = phi.shape()[0];
//     // based on https://www.anycodings.com/1questions/5406370/how-to-use-rayon-for-parallel-calculation-of-pi
//     let summation: f64 = (0..n)
//         .into_par_iter()
//         .map(|j| {
//             let timewithdelay = t - tau[[*i, j]] as usize;
//             let factor = phi[[j, timewithdelay]] - phi_current;
//             ((kcoupling[[*i, j]] * factor.sin() * connectionmatrix[[*i, j]]) as f64) as f64
//         })
//         .reduce(|| 0.0, |a, b| a + b);
//
//     let phi_new = phi_current + (omega[*i] + (summation / n as f64)) * dt; /* error here */
//     phi_new
// }

fn model(
    // parameters for oscillators
    mut phi: Array<f64, Ix2>,
    mut kcoupling: Array<f64, Ix3>,
    omega: &Array<f64, Ix1>,
    // parameters of the simulation itself
    dt: &f64,
    tinitial: &usize,
    tmax: &usize,
    n: &usize,
    epsilon: &f64,
    // parameters for connections
    tau: &Array<usize, Ix2>,
    connectionmatrix: &Array<f64, Ix2>,
    alpha: &Array<f64, Ix2>,
) -> (Array<f64, Ix2>, Array<f64, Ix3>) {
    // replaces "function model(dt, tinitial, tmax, N, epsilon, dimension)" in igor
    // let mut phicurrent: Array<f64, Ix1> = Array::zeros(*n);

    let pb = ProgressBar::new(((*tmax - 1) - *tinitial) as u64);

    // BEGINING OF TIME (1D)
    for t in *tinitial..(*tmax - 1) {
        pb.inc(1); // update progressbar

        // Evolution of the Oscillators
        for i in 0..*n {
            phi[[i, t + 1]] = update_oscillator(
                &t,
                &i,
                &phi,
                &tau,
                &omega,
                &kcoupling,
                &connectionmatrix,
                &dt,
            );
        }

        //Evolution of Coupling
        for i in 0..*n {
            for j in 0..*n {
                let timewithdelay = t - tau[[i, j]] as usize;
                let factor = phi[[i, t]] - phi[[j, timewithdelay]];
                kcoupling[[i, j, t + 1]] = kcoupling[[i, j, t]]
                    + (epsilon * (alpha[[i, j]] * factor.cos() - kcoupling[[i, j, t]]) * dt);
            }
        }
    }
    // END OF TIME (1D)
    pb.finish_with_message("Done with modeling");
    return (phi, kcoupling);
}

fn op_compare(
    tmax: &usize,
    order_parameter: &Array<f64, Ix2>,
    order_parameter2_cluster: &Array<f64, Ix2>,
    order_parameter2_cluster_connected: &Array<f64, Ix2>,
) {
    let mut current_max = 0.0;

    let mut current_string = String::new();
    for m in 0..5 {
        if order_parameter[[m, *tmax - 1]] > current_max {
            current_max = order_parameter[[m, *tmax - 1]];
            current_string =
                format!("It is single cluster, mode={m}, with OrderParameter={current_max}");
        } else if order_parameter2_cluster[[m, *tmax - 1]] > current_max {
            current_max = order_parameter2_cluster[[m, *tmax - 1]];
            current_string =
                format!("It is unconnected 2 cluster, mode={m}, with OrderParameter={current_max}");
        } else if order_parameter2_cluster_connected[[m, *tmax - 1]] > current_max {
            current_max = order_parameter2_cluster_connected[[m, *tmax - 1]];
            current_string =
                format!("It is connected 2 cluster, mode={m}, with OrderParameter={current_max}");
        }
    }
    println!("{}", current_string);
}

pub fn run(
    n: usize,
    timesim: usize,
    dt: f64,
    spreadinomega: f64,

    g: f64,
    epsilon: f64,
    timemetric: f64,

    clustersize: usize,
    drivingfrequency: f64,
) -> (Array<f64, Ix2>, Array<f64, Ix3>) {
    let alpha: Array<f64, _> = Array::ones((n, n)); // coupling speed eviolution
    let connectionmatrix: Array<f64, _> = Array::ones((n, n)); // whether coupled or not

    let omega: Array<f64, Ix1> =
        Array::random(n, rand_distr::Normal::new(1.0, spreadinomega).unwrap()); // random frequencies of oscillators

    // initial conditions
    let tau: Array<usize, Ix2> = calculate_delays(&timemetric, &n, &dt);
    let tinitial: usize = set_tinitial(&tau);
    println!("tinitial = {}", tinitial);
    let tmax: usize = &tinitial + &timesim;

    let mut kcoupling: Array<f64, Ix3> = Array::zeros((n, n, tmax)); // coupling strength
    kcoupling.fill(g);

    let mut phi: Array<f64, Ix2> = initialize_phi(&n, &clustersize, &tmax);
    // println!("phi shape: {:?}", phi.shape()[0]);

    // driving the cluster for the initial time
    phi = driver(
        phi,
        &dt,
        &drivingfrequency,
        &n,
        &clustersize,
        &tinitial,
        &omega,
    );

    // run the actual model
    (phi, kcoupling) = model(
        phi,
        kcoupling,
        &omega,
        &dt,
        &tinitial,
        &tmax,
        &n,
        &epsilon,
        &tau,
        &connectionmatrix,
        &alpha,
    );

    // calculating the order parameter
    let (order_parameter, order_parameter2_cluster, order_parameter2_cluster_connected) =
        calc_order_params(&phi, &n, &tmax);

    op_compare(
        &tmax,
        &order_parameter,
        &order_parameter2_cluster,
        &order_parameter2_cluster_connected,
    );

    // TODO: figure out saving the data
    // possibly consider moving model into a struct
    (phi, kcoupling)
}
