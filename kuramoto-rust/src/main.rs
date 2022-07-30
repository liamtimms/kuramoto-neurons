extern crate ndarray;

use ndarray::prelude::*;
use ndarray::Array;
use ndarray_rand::rand_distr;
use ndarray_rand::RandomExt;

fn max(a: &f64, b: &f64) -> f64 {
    if *a >= *b {
        *a
    } else {
        *b
    }
}

fn calculate_delays(timemetric: &f64, n: &usize, dt: &f64, dimension: &i32) -> Array<usize, Ix2> {
    let z: f64 = timemetric / dt;
    let mut d: f64 = 0.0;

    if *dimension == 1 {
        let mut tau: Array<usize, Ix2> = Array::zeros((*n, *n));
        for i in 0..*n {
            for j in 0..*n {
                if i != j {
                    let x = (i as f64 - j as f64).abs();
                    let y = *n as f64 - x;
                    // pick minimum value
                    if x <= y {
                        d = x;
                    } else if x > y {
                        d = y;
                    } else {
                        panic!("Error in calculating delays");
                    }
                    // input calculated delay into tau
                    tau[[i, j]] = (z / (*n as f64) * d).round() as usize;
                } else {
                    tau[(i, j)] = 0 as usize;
                }
                // tau[i, j]=trunc(Z/N*min(abs(i-j),N-abs(i-j)));
                // tau[[i, j]] = Z / N as f64 * cmp::min(x, N as f64 - x);
                // print!(
                //     "for i={}, j={}, x={}, y={}, d={}, tau = {}\n",
                //     i,
                //     j,
                //     x,
                //     y,
                //     d,
                //     tau[[i, j]]
                // );
            }
        }
        tau
    } else {
        let tau: Array<usize, Ix2> = Array::zeros((*n, *n));
        tau
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
            let mut op_minus = 0.0;
            let mut op_plus = 0.0;
            let mut op_minus2 = 0.0;
            let mut op_plus2 = 0.0;
            let mut op_minus3 = 0.0;
            let mut op_plus3 = 0.0;

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

            op_minus = (sinsum.powi(2) + cossum.powi(2)).sqrt() / *n as f64;
            op_minus2 = ((((sinsum2.powi(2) + cossum2.powi(2)).sqrt()).abs() / *n as f64
                - op_minus)
                .powi(2))
            .sqrt();
            op_minus3 = ((((sinsum3.powi(2) + cossum3.powi(2)).sqrt()).abs() / *n as f64
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
            op_plus = (sinsum.powi(2) + cossum.powi(2)).sqrt() / *n as f64;
            op_plus2 = ((((sinsum2.powi(2) + cossum2.powi(2)).sqrt()).abs() / *n as f64 - op_plus)
                .powi(2))
            .sqrt();
            op_plus3 = ((((sinsum3.powi(2) + cossum3.powi(2)).sqrt()).abs() / *n as f64 - op_plus)
                .powi(2))
            .sqrt();

            sinsum = 0.0;
            cossum = 0.0;
            sinsum2 = 0.0;
            cossum2 = 0.0;
            sinsum3 = 0.0;
            cossum3 = 0.0;

            order_parameter[[m, t]] = max(&op_plus, &op_minus);
            order_parameter2_cluster[[m, t]] = max(&op_plus2, &op_minus2);
            order_parameter2_cluster_connected[[m, t]] = max(&op_plus3, &op_minus3);
        }
    }
    (
        order_parameter,
        order_parameter2_cluster,
        order_parameter2_cluster_connected,
    )
}

fn final_freqs(phi: &Array<f64, Ix2>, n: &usize, dt: &f64, tmax: &usize) -> Array<f64, Ix1> {
    let delta_t = 20;
    let mut finalfrequencies: Array<f64, Ix1> = Array::zeros(*n);
    for i in 0..*n {
        let comparisontime = *tmax - (delta_t + 1);
        finalfrequencies[i] =
            (phi[[i, *tmax - 1]] - phi[[i, comparisontime]]) / (delta_t as f64 * dt);
        finalfrequencies[i] = finalfrequencies[i] * dt;
    }
    finalfrequencies
}

fn initialize_phi(n: &usize, clustersize: &usize, tmax: &usize) -> Array<f64, Ix2> {
    // replaces "function initialconditions1(N,clustersize, dimension)" in igor
    let mut phi: Array<f64, Ix2> = Array::zeros((*n, *tmax));
    for i in 0..*clustersize {
        phi[[i, 0]] = 1.0;
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

fn model(
    // parameters for oscillators
    mut phi: Array<f64, Ix2>,
    mut kcoupling: Array<f64, Ix2>,
    omega: &Array<f64, Ix1>,

    dt: &f64,
    tinitial: &usize,
    tmax: &usize,
    n: &usize,
    epsilon: &f64,

    // parameters for connections
    tau: &Array<usize, Ix2>,
    connectionmatrix: &Array<f64, Ix2>,
    alpha: &Array<f64, Ix2>,
) -> (Array<f64, Ix2>, Array<f64, Ix2>) {
    let mut summation = 0.0;
    let mut phicurrent: Array<f64, Ix1> = Array::zeros(*n);

    //BEGINING OF TIME (1D)
    for t in *tinitial..(*tmax - 1) {
        //Evolution of Phase
        for i in 0..*n {
            phicurrent[i] = phi[[i, t]];
            for j in 0..*n {
                let timewithdelay = t - tau[[i, j]] as usize;
                let factor = phi[[j, timewithdelay]] - phicurrent[i];
                summation = summation + kcoupling[[i, j]] * factor.sin() * connectionmatrix[[i, j]];
            }
            let timewithstep = t + 1;
            phi[[i, timewithstep]] = phicurrent[i] + (omega[i] + (summation / *n as f64)) * dt; /* error here */
            summation = 0.0;
        }

        //Evolution of Coupling
        for i in 0..*n {
            for j in 0..*n {
                let timewithdelay = t - tau[[i, j]] as usize;
                let factor = phicurrent[i] - phi[[j, timewithdelay]];
                kcoupling[[i, j]] = kcoupling[[i, j]]
                    + (epsilon * (alpha[[i, j]] * factor.cos() - kcoupling[[i, j]]) * dt);
            }
        }
    } //END OF TIME (1D)
    return (phi, kcoupling);
}

fn op_compare(
    tmax: &usize,
    order_parameter: &Array<f64, Ix2>,
    order_parameter2_cluster: &Array<f64, Ix2>,
    order_parameter2_cluster_connected: &Array<f64, Ix2>,
) {
    let mut current_max = 0.0;

    for m in 0..5 {
        if order_parameter[[m, *tmax - 1]] > current_max {
            current_max = order_parameter[[m, *tmax - 1]];
        }
    }
}

fn main() {
    // all parameters
    let dimension = 1;
    let g = 1.0; //The initial coupling constant
    let n: usize = 20; //Either the number of oscillators (1D) or the height/width of the square array
    let epsilon = 0.1;
    let dt = 0.01;
    let timemetric: f64 = 2.0; //Zannette's "T" used only for initial conditions if timemetric2 != it.
    let tinitial = 1000;

    let mode = 0; //Sets initial conditions: 0->random phases, 1->same phase cluster, 2->plane wave cluster, 3->circles cluster
    let clustersize = n; //cluster=N for whole array to be driven
    let drivingfrequency: f64 = 1.0; //frequency that the oscillaotrs in the cluster will be forced to move at
    let spreadinomega = 0; //the standard deviation in the natural frequencies
    let tmax = 2000; //the number of total time steps we are simulating (counts the initial conditions)

    if dimension == 1 {
        let mut kcoupling: Array<f64, Ix2> = Array::zeros((n, n));
        kcoupling.fill(g);

        // TODO: figure out the use of these "waves" and variables from old code
        let mut _tt: Array<f64, _> = Array::zeros(n * tmax);
        let mut _freqs: Array<f64, _> = Array::zeros(n * tmax);
        let _epsilon_counter = 0;

        // TODO: better documenation on the differences between these
        let alpha: Array<f64, _> = Array::ones((n, n));
        let connectionmatrix: Array<f64, _> = Array::ones((n, n));

        let tau: Array<usize, Ix2> = calculate_delays(&timemetric, &n, &dt, &dimension);

        // TODO: check implmentation of this normal distribution
        let omega: Array<f64, Ix1> = Array::random(n, rand_distr::StandardNormal);

        // initial conditions
        let mut phi: Array<f64, Ix2> = initialize_phi(&n, &clustersize, &tmax);

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

        let finalfrequencies = final_freqs(&phi, &n, &dt, &tmax);

        // TODO: figure out saving the data
        // possibly consider moving model into a struct
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_finds_max() {
        // my first test!
        assert_eq!(4.0, max(&3.5, &4.0));
    }
}
