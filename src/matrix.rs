use num_traits::Zero;
use std::ops::BitXor;

pub trait Number: Copy + Eq + Ord + BitXor<Self> + Zero {}

impl Number for u8 {}

pub trait MatrixTrait<T: Number>: HasElements<T> {
    fn rank(&self) -> usize;
    fn kernel(&self) -> Vec<Vec<T>>;
    fn echelon_form(&self) -> (Self, Vec<(usize, usize)>)
    where
        Self: Sized;
    fn image(&self) -> Vec<Vec<T>>;
    fn is_reduced_echelon(&self) -> bool;
    fn nrows(&self) -> usize {
        self.elements().len()
    }
    fn ncols(&self) -> usize {
        self.elements().get(0).map_or(0, |row| row.len())
    }
    fn get_pivot(vec: &Vec<T>) -> Option<usize> {
        vec.iter().position(|&x| !x.is_zero())
    }
}

pub trait HasElements<T: Number> {
    fn elements(&self) -> &Vec<Vec<T>>;
}

#[derive(Clone, Debug)]
pub struct Matrix<T: Number> {
    pub elements: Vec<Vec<T>>,
}

impl<T: Number> Matrix<T> {
    pub fn new(elements: Vec<Vec<T>>) -> Self {
        Self { elements }
    }
}

impl<T: Number> HasElements<T> for Matrix<T> {
    fn elements(&self) -> &Vec<Vec<T>> {
        &self.elements
    }
}
