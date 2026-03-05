use crate::{matrix::Number, GF2Matrix};

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
#[derive(Clone)]
pub struct PackedGF2Matrix<T: Number> {
    elements: Vec<T>,
    n: usize,
}


impl<T: Number> PackedGF2Matrix<T> {
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
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    ///
    /// let m = PackedGF2Matrix::new(vec![0b1011u8, 0b0101u8], 4);
    /// ```
    pub fn new(elements: Vec<T>, n: usize) -> Self {
        Self {
            elements: elements,
            n: n,
        }
    }

    /// Returns the number of rows in the matrix.
    ///
    /// # Returns
    ///
    /// number of rows: this is equal to the number of elements stored internally.
    pub fn nrows(&self) -> usize {
        self.elements.len()
    }

    /// Returns the number of columns in the matrix.
    ///
    /// # Returns
    ///
    /// Number of columns: this corresponds to the number of bits extracted from each row.
    pub fn ncols(&self) -> usize {
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
    /// # use lin_algebra::packed_gf2_matrix::{PackedGF2Matrix, BitOrder};
    ///
    /// let m = PackedGF2Matrix::new(vec![0b0010u8], 4);
    /// let gf2 = m.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
    ///
    /// // Result: [[0, 1, 0, 0]]
    /// ```
    pub fn from_int_matrix_to_gf2_matrix(&self, bit_order: BitOrder) -> GF2Matrix {
        let mut matrix_elements = Vec::with_capacity(self.nrows());
        match bit_order {
            BitOrder::LSB => {
                for &row in &self.elements {
                    let mut r = Vec::with_capacity(self.n);
                    for i in 0..self.n {
                        r.push((((row >> i) & T::one()) != T::zero()) as u8);
                    }
                    matrix_elements.push(r);
                }
            }
            BitOrder::MSB => {
                for &row in &self.elements {
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

    pub fn from_vec(vect: Vec<T>) -> Self {
        let n = vect.len();
        Self {
            elements: vect,
            n: n,
        }
    }

    pub fn from_vec_referenced(vect: &Vec<T>) -> Self {
        let n = vect.len();
        Self {
            elements: vect.to_vec(),
            n: n,
        }
    }

    fn get_element(&self, row: usize, col: usize) -> u8 {
        ((self.elements[row].into_usize() >> (self.ncols() - 1 - col)) & 1) as u8
    }

    fn swap_rows(&mut self, idx1: usize, idx2: usize) {
        self.elements.swap(idx1, idx2);
    }

    pub fn echelon_form(&self) -> (Self, Vec<(usize, usize)>) {
        let mut m_copy = self.clone();
        let mut lead = 0;
        let mut operations: Vec<(usize, usize)> = Vec::new();
        for r in 0..self.nrows() {
            if lead >= self.ncols() {
                break;
            }
            let mut i = r;
            while m_copy.get_element(i, lead) == 0 {
                i += 1;
                if i == self.nrows() {
                    i = r;
                    lead += 1;
                    if lead == self.ncols() {
                        return (m_copy, operations);
                    }
                }
            }
            m_copy.swap_rows(r, i);
            if r != i {
                operations.push((r, i));
                operations.push((i, r));
                operations.push((r, i));
            }
            for i in 0..self.nrows() {
                if i != r && m_copy.get_element(i, lead) == 1 {
                    m_copy.elements[i] = m_copy.elements[i] ^ m_copy.elements[r];
                    operations.push((i, r));
                }
            }
            lead += 1
        }
        (m_copy, operations)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_from_int_matrix_to_gf2_matrix_u8_lsb() {
        let elements = vec![0, 1, 2, 4, 8];
        let int_matrix = PackedGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
        let expected = vec![
            vec![0, 0, 0, 0],
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 0, 0, 1],
        ];
        assert_eq!(gf2_matrix.elements, expected);

        let elements = vec![0, 1, 2, 3, 5];
        let int_matrix = PackedGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
        let expected = vec![
            vec![0, 0, 0, 0],
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![1, 1, 0, 0],
            vec![1, 0, 1, 0],
        ];
        assert_eq!(gf2_matrix.elements, expected);
    }

    #[test]
    fn test_from_int_matrix_to_gf2_matrix_u8_msb() {
        let elements = vec![0, 1, 2, 4, 8];
        let int_matrix = PackedGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::MSB);
        let expected = vec![
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 1],
            vec![0, 0, 1, 0],
            vec![0, 1, 0, 0],
            vec![1, 0, 0, 0],
        ];
        assert_eq!(gf2_matrix.elements, expected);

        let elements = vec![0, 1, 2, 3, 5];
        let int_matrix = PackedGF2Matrix::<u8>::new(elements.clone(), 4);
        let gf2_matrix = int_matrix.from_int_matrix_to_gf2_matrix(BitOrder::MSB);
        let expected = vec![
            vec![0, 0, 0, 0],
            vec![0, 0, 0, 1],
            vec![0, 0, 1, 0],
            vec![0, 0, 1, 1],
            vec![0, 1, 0, 1],
        ];
        assert_eq!(gf2_matrix.elements, expected);
    }
}
