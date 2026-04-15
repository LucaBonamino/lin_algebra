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

    /// Creates a one-row bit-packed matrix from an owned vector of row data.
    ///
    /// Each element of `vect` is interpreted as one packed row of the matrix.
    /// The number of columns is inferred from the length of the vector.
    ///
    /// # Parameters
    ///
    /// - `vect`: the vector containing the packed row elements.
    ///
    /// # Returns
    ///
    /// A `PackedGF2Matrix` whose internal storage is initialized from `vect`.
    ///
    /// # Notes
    ///
    /// - The input vector is moved into the matrix without cloning.
    /// - The resulting matrix has `vect.len()` columns in its packed representation.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// # use lin_algebra::matrix::Number;
    /// let m = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    /// let (echelon, ops) = m.echelon_form();
    ///
    /// assert_eq!(m.ncols(), 3);
    /// ```
    pub fn from_vec(vect: Vec<T>) -> Self {
        let n = vect.len();
        Self {
            elements: vect,
            n: n,
        }
    }

    /// Creates a one-row bit-packed matrix from a referenced vector of row data.
    ///
    /// Each element of `vect` is interpreted as one packed row of the matrix.
    /// The number of columns is inferred from the length of the vector.
    ///
    /// # Parameters
    ///
    /// - `vect`: a reference to the vector containing the packed row elements.
    ///
    /// # Returns
    ///
    /// A `PackedGF2Matrix` whose internal storage is initialized from a clone
    /// of `vect`.
    ///
    /// # Notes
    ///
    /// - The input vector is cloned using `to_vec()`.
    /// - The resulting matrix owns its internal storage independently of the
    ///   original vector.
    /// - The resulting matrix has `vect.len()` columns in its packed representation.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let rows = vec![0b1010u8, 0b0110u8];
    /// let m = PackedGF2Matrix::from_vec_referenced(&rows);
    ///
    /// assert_eq!(m.ncols(), 2);
    /// ```
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

    /// Computes the row echelon form of the bit-packed matrix over GF(2).
    ///
    ///
    /// In addition to the transformed matrix, it records the sequence of
    /// elementary row operations applied during the elimination process.   
    ///
    /// # Row Operations
    ///
    /// The following row operations may be recorded:
    ///
    /// - `(i, j)` with `i != j` represents adding row `j` to row `i`
    ///   (that is, `row_i <- row_i + row_j` over GF(2)).
    /// - A row swap between `r` and `i` is encoded as three additions:
    ///   `(r, i)`, `(i, r)`, `(r, i)`.
    ///
    /// This encoding is valid over GF(2), where swapping two rows can be
    /// expressed as a sequence of XOR-based row additions.
    ///
    /// # Returns
    ///
    /// A pair `(echelon_matrix, operations)` where:
    ///
    /// - `echelon_matrix` is the row echelon form of the matrix.
    /// - `operations` is the list of row operations used to obtain it.
    ///
    /// # Notes
    ///
    /// - The computation is performed on a clone of the matrix; the original
    ///   matrix is not modified.
    /// - Since the matrix is over GF(2), pivots are always `1`.
    /// - The elimination clears all other `1`s in each pivot column, so the
    ///   result is closer to a reduced row echelon form than to a strictly
    ///   upper-triangular echelon form.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// # use lin_algebra::matrix::Number;
    /// let m = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    /// let (echelon, ops) = m.echelon_form();
    /// // `echelon` contains the transformed matrix over GF(2),
    /// // and `ops` contains the recorded row operations.
    /// ```
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
