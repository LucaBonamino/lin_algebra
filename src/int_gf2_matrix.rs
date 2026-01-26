use crate::{GF2Matrix, matrix::Number};

#[derive(Clone, Copy, Debug)]
pub enum BitOrder {
    LSB,
    MSB,
}

/// An intermediate matrix representation where each row is encoded
/// as an integer value.
///
/// Each element of `elements` represents one row of the matrix,
/// with its bits encoding the entries of that row.
///
/// This type is typically used as a compact or efficient
/// representation before expanding into an explicit GF(2) matrix.
pub struct InterGF2Matrix<T: Number> {
    elements: Vec<T>,
    n: usize
}

impl<T: Number> InterGF2Matrix<T>{

    /// Creates a new integer-encoded matrix.
    ///
    /// # Arguments
    ///
    /// * `elements` - A vector where each element encodes one row as bits.
    /// * `n` - The number of columns (bits) per row.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::int_gf2_matrix::InterGF2Matrix;
    /// 
    /// let m = InterGF2Matrix::new(vec![0b1011u8, 0b0101u8], 4);
    /// ```
    pub fn new(elements: Vec<T>, n: usize) -> Self{
        Self {elements: elements, n: n}
    }

    /// Returns the number of rows in the matrix.
    ///
    /// # Returns
    /// 
    /// number of rows: this is equal to the number of elements stored internally.
    pub fn nrows(&self) -> usize{
        self.elements.len()
    }

    /// Returns the number of columns in the matrix.
    ///
    /// # Returns
    /// 
    /// Number of columns: this corresponds to the number of bits extracted from each row.
    pub fn ncols(&self) -> usize{
        self.n
    }

    /// Returns the integer-encoded value of a specific row.
    ///
    /// # Arguments
    /// 
    /// * `row_index`- usize, index of the row.
    /// 
    /// # Returns
    /// 
    /// Row corresponding to index `row_index`
    /// 
    /// # Panics
    ///
    /// Panics if `row_index` is out of bounds.
    pub fn row(&self, row_index: usize) -> T {
        self.elements[row_index]
    }


    /// Converts the integer-encoded matrix into an explicit GF(2) matrix.
    ///
    /// Each row is expanded into a vector of bits (`0` or `1`),
    /// according to the specified bit order.
    ///
    /// # Bit Order
    ///
    /// - `BitOrder::LSB`: the least-significant bit is placed first
    ///   (column 0 corresponds to bit 0).
    /// - `BitOrder::MSB`: the most-significant bit is placed first
    ///   (column 0 corresponds to bit `n - 1`).
    ///
    /// # Returns
    ///
    /// A `GF2Matrix` whose entries are elements of GF(2),
    /// represented as `u8` values (`0` or `1`).
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::int_gf2_matrix::{InterGF2Matrix, BitOrder};
    /// 
    /// let m = InterGF2Matrix::new(vec![0b0010u8], 4);
    /// let gf2 = m.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
    ///
    /// // Result: [[0, 1, 0, 0]]
    /// ```
    pub fn from_int_matrix_to_gf2_matrix(&self, bit_order: BitOrder) -> GF2Matrix {

        let mut matrix_elements = Vec::with_capacity(self.nrows());
        match bit_order {
            BitOrder::LSB => {
                for &row in &self.elements{
                    let mut r = Vec::with_capacity(self.n);
                    for i in 0..self.n {
                        r.push((((row >> i) & T::one()) != T::zero()) as u8);
                    }
                    matrix_elements.push(r);
                    }
                },
            BitOrder::MSB => {
                for &row in &self.elements{
                    let mut r = Vec::with_capacity(self.n);
                    for i in (0..self.n).rev() {
                        r.push((((row >> i) & T::one()) != T::zero()) as u8);
                    }
                    matrix_elements.push(r);
                }
            }
        };
        GF2Matrix::new(matrix_elements)

    }
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_from_int_matrix_to_gf2_matrix_u8_lsb(){
        let elements = vec![0,1,2,4,8];
        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
        let expected = vec![
            vec![0,0,0,0],
            vec![1,0,0,0],
            vec![0,1,0,0],
            vec![0,0,1,0],
            vec![0,0,0,1],
        ];
        assert_eq!(gf2_matrix.elements, expected);

        let elements = vec![0,1,2,3,5];
        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
        let expected = vec![
            vec![0,0,0,0],
            vec![1,0,0,0],
            vec![0,1,0,0],
            vec![1,1,0,0],
            vec![1,0,1,0],
        ];
        assert_eq!(gf2_matrix.elements, expected);
    }

    #[test]
    fn test_from_int_matrix_to_gf2_matrix_u8_msb(){
        let elements = vec![0,1,2,4, 8];
        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::MSB);
        let expected = vec![
            vec![0,0,0,0],
            vec![0,0,0,1],
            vec![0,0,1,0],
            vec![0,1,0,0],
            vec![1,0,0,0],
        ];
        assert_eq!(gf2_matrix.elements, expected);

        let elements = vec![0,1,2,3,5];
        let int_matrix = InterGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::MSB);
        let expected = vec![
            vec![0,0,0,0],
            vec![0,0,0,1],
            vec![0,0,1,0],
            vec![0,0,1,1],
            vec![0,1,0,1],
        ];
        assert_eq!(gf2_matrix.elements, expected);
    }
}