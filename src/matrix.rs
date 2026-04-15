use num_traits::{One, Zero};
use std::ops::{BitAnd, BitXor, Shr};

pub trait Number:
    Copy
    + Eq
    + Ord
    + BitXor<Output = Self>
    + BitAnd<Output = Self>
    + Shr<usize, Output = Self>
    + Zero
    + One
{
    fn into_usize(self) -> usize;
}

impl Number for u8 {
    fn into_usize(self) -> usize {
        self as usize
    }
}
impl Number for u16 {
    fn into_usize(self) -> usize {
        self as usize
    }
}
impl Number for u32 {
    fn into_usize(self) -> usize {
        self as usize
    }
}
impl Number for u64 {
    fn into_usize(self) -> usize {
        self as usize
    }
}

pub trait MatrixTrait<T: Number>: MatrixCommon<T> {
    fn rank(&self) -> usize;
    fn kernel(&self) -> Vec<Vec<T>>;
    fn echelon_form(&self) -> (Self, Vec<(usize, usize)>)
    where
        Self: Sized;
    fn image(&self) -> Vec<Vec<T>>;
    fn is_reduced_echelon(&self) -> bool;
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

impl<T: Number> MatrixCommon<T> for Matrix<T> {
    fn nrows(&self) -> usize {
        self.elements.len()
    }
    fn ncols(&self) -> usize {
        self.elements.first().map_or(0, |r| r.len())
    }
    fn row(&self, r: usize) -> &[T] {
        &self.elements[r]
    }
}

pub trait MatrixCommon<T: Number> {
    fn nrows(&self) -> usize;
    fn ncols(&self) -> usize;
    fn row(&self, r: usize) -> &[T];
    fn get_pivot(row: &[T]) -> Option<usize> {
        row.iter().position(|&x| !x.is_zero())
    }
}
