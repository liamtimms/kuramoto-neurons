use ndarray::prelude::*;
use ndarray::{Array, ErrorKind, ShapeError};

// Help from friend on how to do this
pub enum MyArray<A> {
    Array2(Array2<A>),
    Array3(Array3<A>),
}

// ArrayD is a dynamic-dimensional array from ndarray
impl<A> From<MyArray<A>> for ArrayD<A> {
    fn from(arr: MyArray<A>) -> Self {
        match arr {
            MyArray::Array2(array2) => array2.into_dyn(),
            MyArray::Array3(array3) => array3.into_dyn(),
        }
    }
}

// need to understand array traits better
impl<A> TryFrom<ArrayD<A>> for MyArray<A> {
    type Error = ShapeError;

    fn try_from(arr: ArrayD<A>) -> Result<MyArray<A>, Self::Error> {
        match arr.ndim() {
            2 => arr.into_dimensionality::<Ix2>().map(MyArray::Array2),
            3 => arr.into_dimensionality::<Ix3>().map(MyArray::Array3),
            _ => Err(ShapeError::from_kind(ErrorKind::IncompatibleShape)),
        }
    }
}

// really high level stuff here
impl<A: Default> MyArray<A> {
    fn apply<R>(&mut self, func: impl FnOnce(&mut ArrayD<A>) -> R) -> R {
        // Calling 'Into::into' consumes "self"
        // we make a copy of "self" to avoid consuming it
        let my_array = std::mem::replace(self, MyArray::Array2(Array2::default([0, 0])));
        let mut dyn_array: ArrayD<A> = my_array.into();
        let result = func(&mut dyn_array);
        *self = dyn_array.try_into().unwrap();
        result
    }
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
