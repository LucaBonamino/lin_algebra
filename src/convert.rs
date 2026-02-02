use crate::gf2_matrix::GF2Matrix;
use crate::int_gf2_matrix::{BitOrder, InterGF2Matrix};
use crate::matrix::{Matrix, Number};

impl<T: Number> From<&InterGF2Matrix<T>> for GF2Matrix {
    fn from(int_matrix: &InterGF2Matrix<T>) -> Self {
        int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::MSB)
    }
}

impl<T: Number> From<InterGF2Matrix<T>> for GF2Matrix {
    fn from(int_matrix: InterGF2Matrix<T>) -> Self {
        int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::MSB)
    }
}

impl<T: Number> Matrix<T>
where
    T: Number + From<usize>,
{
    pub fn from_vec_to_int_msb(v: &[u8]) -> T {
        let mut x = 0usize;
        let n = v.len();

        for (i, &bit) in v.iter().enumerate() {
            assert!(bit == 0 || bit == 1);
            if bit == 1 {
                x |= 1 << (n - 1 - i);
            }
        }

        T::from(x)
    }
}

impl<T: Number> From<GF2Matrix> for InterGF2Matrix<T>
where
    T: Number + From<usize>,
{
    fn from(gf2_matrix: GF2Matrix) -> Self {
        let ve: Vec<T> = gf2_matrix
            .elements
            .iter()
            .map(|row| Matrix::<T>::from_vec_to_int_msb(&row))
            .collect();
        Self::from_vec(ve)
    }
}

impl<T: Number> From<&GF2Matrix> for InterGF2Matrix<T>
where
    T: Number + From<usize>,
{
    fn from(gf2_matrix: &GF2Matrix) -> Self {
        let ve: Vec<T> = gf2_matrix
            .elements
            .iter()
            .map(|row| Matrix::<T>::from_vec_to_int_msb(&row))
            .collect();
        Self::from_vec(ve)
    }
}

impl<T: Number> From<Vec<T>> for InterGF2Matrix<T> {
    fn from(value: Vec<T>) -> Self {
        Self::from_vec(value)
    }
}

impl<T: Number> From<&Vec<T>> for InterGF2Matrix<T> {
    fn from(value: &Vec<T>) -> Self {
        Self::from_vec_referenced(value)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_int_mtrix_to_gf2_matrix_convertion() {
        let elements = vec![0, 1, 2, 4, 8];
        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = GF2Matrix::from(int_matrix);
        let expected = vec![
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 1],
            vec![0, 0, 1, 0],
            vec![0, 1, 0, 0],
            vec![1, 0, 0, 0],
        ];
        assert_eq!(gf2_matrix.elements, expected);

        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix: GF2Matrix = int_matrix.into();
        assert_eq!(gf2_matrix.elements, expected);
    }

    #[test]
    fn test_int_mtrix_to_gf2_matrix_convertion_by_ref() {
        let elements = vec![0, 1, 2, 4, 8];
        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = GF2Matrix::from(&int_matrix);
        let expected = vec![
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 1],
            vec![0, 0, 1, 0],
            vec![0, 1, 0, 0],
            vec![1, 0, 0, 0],
        ];
        assert_eq!(gf2_matrix.elements, expected);

        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix: GF2Matrix = (&int_matrix).into();
        assert_eq!(gf2_matrix.elements, expected);
    }
}
