#[macro_use]
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

fn calculate_delays(timemetric: &f64, N: &usize, dt: &f64, dimension: &i32) -> Array<f64, Ix2> {
    let Z: f64 = timemetric / dt;
    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut d: f64 = 0.0;

    if *dimension == 1 {
        let mut tau: Array<f64, Ix2> = Array::zeros((*N, *N));
        for i in 0..*N {
            for j in 0..*N {
                if i != j {
                    x = (i as f64 - j as f64).abs();
                    y = *N as f64 - x;
                    // pick minimum value
                    if x <= y {
                        d = x;
                    } else if x > y {
                        d = y;
                    }
                    // input calculated delay into tau
                    tau[[i, j]] = Z / (*N as f64) * d;
                } else {
                    tau[(i, j)] = 0.0;
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
        let mut tau: Array<f64, Ix2> = Array::zeros((*N, *N));
        tau
    }
}

fn OPs(
    phi: Array<f64, Ix2>,
    N: &usize,
    tmax: &usize,
) -> (Array<f64, Ix2>, Array<f64, Ix2>, Array<f64, Ix2>) {
    // calculate order parameters
    // This is reproduced from the original code
    // and clearly not optimal.
    let mut sinsum = 0.0;
    let mut cossum = 0.0;

    let mut sinsum2 = 0.0;
    let mut cossum2 = 0.0;

    let mut sinsum3 = 0.0;
    let mut cossum3 = 0.0;

    let mut OPminus = 0.0;
    let mut OPplus = 0.0;

    let mut OPminus2 = 0.0;
    let mut OPplus2 = 0.0;

    let mut OPminus3 = 0.0;
    let mut OPplus3 = 0.0;

    let mut phicurrent: Array<f64, Ix1> = Array::zeros(*N);

    let mut OrderParameter: Array<f64, Ix2> = Array::zeros((5, *tmax));
    let mut OrderParameter2Cluster: Array<f64, Ix2> = Array::zeros((5, *tmax));
    let mut OrderParameter2ClusterConnected: Array<f64, Ix2> = Array::zeros((5, *tmax));

    for m in 0..5 {
        for t in 0..*tmax {
            for i in 0..*N {
                phicurrent[i] =
                    phi[[i, t]] - (2.0 * std::f64::consts::PI * m as f64) * i as f64 / *N as f64;
                sinsum += (phicurrent[i]).sin();
                cossum += (phicurrent[i]).sin();
                sinsum2 += (2.0 * phicurrent[i]).sin();
                cossum2 += (2.0 * phicurrent[i]).cos();
                sinsum3 += (3.0 * phicurrent[i]).sin();
                cossum3 += (3.0 * phicurrent[i]).cos();
            }

            OPminus = (sinsum.powi(2) + cossum.powi(2)).sqrt() / *N as f64;
            OPminus2 = ((((sinsum2.powi(2) + cossum2.powi(2)).sqrt()).abs() / *N as f64 - OPminus)
                .powi(2))
            .sqrt();
            OPminus3 = ((((sinsum3.powi(2) + cossum3.powi(2)).sqrt()).abs() / *N as f64 - OPminus)
                .powi(2))
            .sqrt();

            sinsum = 0.0;
            cossum = 0.0;
            sinsum2 = 0.0;
            cossum2 = 0.0;
            sinsum3 = 0.0;
            cossum3 = 0.0;

            for i in 0..*N {
                phicurrent[i] =
                    phi[[i, t]] + (2.0 * std::f64::consts::PI * m as f64) * i as f64 / *N as f64;
                sinsum += (phicurrent[i]).sin();
                cossum += (phicurrent[i]).sin();
                sinsum2 += (2.0 * phicurrent[i]).sin();
                cossum2 += (2.0 * phicurrent[i]).cos();
                sinsum3 += (3.0 * phicurrent[i]).sin();
                cossum3 += (3.0 * phicurrent[i]).cos();
            }
            OPplus = (sinsum.powi(2) + cossum.powi(2)).sqrt() / *N as f64;
            OPplus2 = ((((sinsum2.powi(2) + cossum2.powi(2)).sqrt()).abs() / *N as f64 - OPplus)
                .powi(2))
            .sqrt();
            OPplus3 = ((((sinsum3.powi(2) + cossum3.powi(2)).sqrt()).abs() / *N as f64 - OPplus)
                .powi(2))
            .sqrt();

            sinsum = 0.0;
            cossum = 0.0;
            sinsum2 = 0.0;
            cossum2 = 0.0;
            sinsum3 = 0.0;
            cossum3 = 0.0;

            OrderParameter[[m, t]] = max(&OPplus, &OPminus);
            OrderParameter2Cluster[[m, t]] = max(&OPplus2, &OPminus2);
            OrderParameter2ClusterConnected[[m, t]] = max(&OPplus3, &OPminus3);

            OPminus = 0.0;
            OPplus = 0.0;
            OPminus2 = 0.0;
            OPplus2 = 0.0;
            OPminus3 = 0.0;
            OPplus3 = 0.0;
        }
    }
    (
        OrderParameter,
        OrderParameter2Cluster,
        OrderParameter2ClusterConnected,
    )
}

fn initialize_phi(N: &usize, clustersize: &usize, tmax: &usize) -> Array<f64, Ix2> {
    // replaces "function initialconditions1(N,clustersize, dimension)" in igor
    let mut phi: Array<f64, Ix2> = Array::zeros((*N, *tmax));
    for i in 0..*clustersize {
        phi[[i, 0]] = 1.0;
    }
    phi
}

fn driver(
    mut phi: Array<f64, Ix2>,
    dt: &f64,
    drivingfrequency: &f64,
    N: &usize,
    clustersize: &usize,
    tinitial: &usize,
    omega: &Array<f64, Ix1>,
) -> Array<f64, Ix2> {
    // replaces "function driver(dt, drivingfrequency, clustersize, dimension)" in igor
    for t in 0..*tinitial {
        for i in 0..*N {
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
    mut phi: Array<f64, Ix2>,
    mut Kcoupling: Array<f64, Ix2>,
    dt: &f64,
    tinitial: &usize,
    tmax: &usize,
    N: &usize,
    epsilon: &f64,
    tau: &Array<f64, Ix2>,
    omega: &Array<f64, Ix1>,
    connectionmatrix: &Array<f64, Ix2>,
    alpha: &Array<f64, Ix2>,
) -> (Array<f64, Ix2>, Array<f64, Ix2>) {
    let mut summation = 0.0;
    let mut phicurrent: Array<f64, Ix1> = Array::zeros(*N);
    let mut factor = 0.0;
    let mut timewithdelay = 0;
    let mut timewithstep = 0;

    //BEGINING OF TIME (1D)
    for t in *tinitial..(*tmax - 1) {
        //Evolution of Phase
        for i in 0..*N {
            phicurrent[i] = phi[[i, t]];
            for j in 0..*N {
                timewithdelay = t - tau[[i, j]] as usize;
                factor = phi[[j, timewithdelay]] - phicurrent[i];
                summation = summation + Kcoupling[[i, j]] * factor.sin() * connectionmatrix[[i, j]];
            }
            timewithstep = t + 1;
            phi[[i, timewithstep]] = phicurrent[i] + (omega[i] + (summation / *N as f64)) * dt; /* error here */
            summation = 0.0;
        }

        //Evolution of Coupling
        for i in 0..*N {
            for j in 0..*N {
                timewithdelay = t - tau[[i, j]] as usize;
                factor = phicurrent[i] - phi[[j, timewithdelay]];
                Kcoupling[[i, j]] = Kcoupling[[i, j]]
                    + (epsilon * (alpha[[i, j]] * factor.cos() - Kcoupling[[i, j]]) * dt);
            }
        }
    } //END OF TIME (1D)
    return (phi, Kcoupling);
}

fn OPcompare(
    tmax: &usize,
    OrderParameter: &Array<f64, Ix2>,
    OrderParameter2Cluster: &Array<f64, Ix2>,
    OrderParameter2ClusterConnected: &Array<f64, Ix2>,
) {
    let mut currentMax = 0.0;

    for m in 0..5 {
        if OrderParameter[[m, *tmax - 1]] > currentMax {
            currentMax = OrderParameter[[m, *tmax - 1]];
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
    let drivingfrequency: f64 = 1.0; //frequency that the oscillaotrs in the cluster will be forced to move at
    let spreadinomega = 0; //the standard deviation in the natural frequencies
    let tmax = 2000; //the number of total time steps we are simulating (counts the initial conditions)

    if dimension == 1 {
        let epsilon_counter = 0;

        let mut Kcoupling: Array<f64, Ix2> = Array::zeros((N, N));
        let mut FinalFrequencies: Array<f64, _> = Array::zeros(N);
        let mut tt: Array<f64, _> = Array::zeros(N * tmax);
        let mut freqs: Array<f64, _> = Array::zeros(N * tmax);

        let mut alpha: Array<f64, _> = Array::ones((N, N));
        let mut connectionmatrix: Array<f64, _> = Array::ones((N, N));
        Kcoupling.fill(g);

        let tau: Array<f64, Ix2> = calculate_delays(&timemetric, &N, &dt, &dimension);
        let omega: Array<f64, Ix1> = Array::random(N, rand_distr::StandardNormal);
        // initial conditions
        let mut phi: Array<f64, Ix2> = initialize_phi(&N, &clustersize, &tmax);

        // driving the cluster for the initial time
        phi = driver(
            phi,
            &dt,
            &drivingfrequency,
            &N,
            &clustersize,
            &tinitial,
            &omega,
        );

        // run the actual model
        (phi, Kcoupling) = model(
            phi,
            Kcoupling,
            &dt,
            &tinitial,
            &tmax,
            &N,
            &epsilon,
            &tau,
            &omega,
            &connectionmatrix,
            &alpha,
        );

        // calculating the order parameter
        let (OrderParameter, OrderParameter2Cluster, OrderParameter2ClusterConnected) =
            OPs(phi, &N, &tmax);
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
