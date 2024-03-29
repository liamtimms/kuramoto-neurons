use std::path::Path;
use std::path::PathBuf;

pub struct RunParams {
    pub n: usize,
    pub timesim: usize,
    pub dt: f64,
    pub spreadinomega: f64,
    pub g: f64,
    pub epsilon: f64,
    pub timemetric: f64,
    pub clustersize: usize,
    pub drivingfrequency: f64,
}

pub fn update_name(output_dir: &Path, base_name: &str, identifier: &str) -> PathBuf {
    let mut new_name = base_name.to_owned();
    new_name.push_str(identifier);
    output_dir.join(new_name)
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
