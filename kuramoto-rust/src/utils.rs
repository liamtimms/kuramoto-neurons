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

pub fn find_delay(i: &usize, j: &usize, n: &usize, z: &f64) -> usize {
    // find the delay between oscillator i and it's neighbor j
    let x = (*i as f64 - *j as f64).abs();
    let y = *n as f64 - x; // wrap around otherside of the circle
    let d = {
        // find minimum distance
        if x <= y {
            x
        } else if x > y {
            y
        } else {
            0.0
        }
    };
    ((z / *n as f64) * d).round() as usize
}
