use ndarray::prelude::*;
use ndarray::Array;

// Currently unused struct
struct KuramotoCircle {
    // a 1D circle of oscillators (2nd dimension is the time)
    phi: Array<f64, Ix2>,       // the phase of each oscillator in rads
    omega: Array<f64, Ix2>,     // the angular velocity of each oscillator
    kcoupling: Array<f64, Ix2>, // the coupling constant of each oscillator
    alpha: f64,                 // a constant that controls the strength of the coupling
    epsilon: f64,               // constant to adjust dynamic level of coupling
}

pub fn max(a: &f64, b: &f64) -> f64 {
    // quick function to return the max of two f64 values
    if *a >= *b {
        *a
    } else {
        *b
    }
}

pub fn min(x: f64, x_prime: f64) -> f64 {
    if x <= x_prime {
        x
    } else if x > x_prime {
        x_prime
    } else {
        0.0
    }
}
