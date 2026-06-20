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

    fn get_packed_bit(value: T, len: usize, idx: usize) -> u8 {
        ((value.into_usize() >> (len - 1 - idx)) & 1) as u8
    }

    fn toggle_packed_bit(value: &mut T, len: usize, idx: usize) {
        *value = *value ^ (T::one() << (len - 1 - idx));
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

    /// Reduces this matrix in place to row-reduced echelon form over GF(2).
    ///
    /// This method modifies the matrix directly and returns the ordered list of row
    /// operations used during elimination.
    ///
    /// # Row Operations
    ///
    /// The returned operations are ordered and can be applied to the original matrix
    /// to reproduce the final state of `self`.
    ///
    /// The following row operations may be recorded:
    ///
    ///  - `(i, j)` with `i != j` represents adding row `j` to row `i`: `row_i <- row_i + row_j`.
    /// - A row swap between rows `r` and `i` is encoded as three row additions:
    ///  `(r, i)`, `(i, r)`, `(r, i)`.
    ///
    /// # Notes
    /// - This method does not clone the matrix.
    /// - The original contents of the matrix are overwritten.
    /// - This method allocates a new `Vec` for the returned operations.
    /// - If you want to reuse the operations vector allocation across repeated
    /// calls, use [`Self::echelon_form_in_place_with_ops`].
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let mut m = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    ///
    /// let ops = m.echelon_form_in_place();
    ///
    /// // `m` is now in row-reduced echelon form.
    /// // `ops` contains the recorded row operations.
    /// ```
    pub fn echelon_form_in_place(&mut self) -> Vec<(usize, usize)> {
        let mut operations = Vec::new();
        self.echelon_form_in_place_with_ops(&mut operations);
        operations
    }

    /// Reduces this matrix in place to row-reduced echelon form over GF(2),
    /// reusing the provided operations buffer.
    /// This method is intended for performance-sensitive code that repeatedly
    /// computes row-reduced echelon forms and wants to avoid allocating a new
    /// operations vector on every call.
    ///
    /// The provided `operations` vector is cleared before new operations are
    /// recorded. Its existing capacity is retained.
    ///
    /// # Row operations
    ///
    /// The recorded operations are ordered and can be applied to the original matrix
    /// to reproduce the final state of `self`.
    ///
    /// The following row operations may be recorded:
    ///
    ///  - `(i, j)` with `i != j` represents adding row `j` to row `i`: `row_i <- row_i + row_j`.
    ///  - A row swap between rows `r` and `i` is encoded as three row additions:
    /// `(r, i)`, `(i, r)`, `(r, i)`.
    ///
    ///  # Notes
    ///
    ///  - This method does not clone the matrix.
    /// - The original contents of the matrix are overwritten.
    /// - `operations.clear()` is called before recording new operations.
    /// - Clearing the vector removes old operations but keeps the allocated capacity for reuse.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let original = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    /// let mut work = original.clone();
    /// let mut ops = Vec::new();
    ///
    /// work.clone_from(&original);
    /// work.echelon_form_in_place_with_ops(&mut ops);
    /// // `work` is now in row-reduced echelon form.
    /// // `ops` contains the recorded row operations.
    /// ```
    pub fn echelon_form_in_place_with_ops(&mut self, operations: &mut Vec<(usize, usize)>) {
        operations.clear();

        let mut lead = 0;

        for r in 0..self.nrows() {
            if lead >= self.ncols() {
                break;
            }

            let mut i = r;

            while self.get_element(i, lead) == 0 {
                i += 1;

                if i == self.nrows() {
                    i = r;
                    lead += 1;

                    if lead == self.ncols() {
                        return;
                    }
                }
            }

            self.swap_rows(r, i);

            if r != i {
                operations.push((r, i));
                operations.push((i, r));
                operations.push((r, i));
            }

            for i in 0..self.nrows() {
                if i != r && self.get_element(i, lead) == 1 {
                    self.elements[i] = self.elements[i] ^ self.elements[r];
                    operations.push((i, r));
                }
            }

            lead += 1;
        }
    }

    /// Converts this matrix into row-reduced echelon form over GF(2).
    ///
    /// This method consumes the matrix, performs Gauss-Jordan-style elimination
    /// in place, and returns the transformed matrix together with the recorded row operations.
    ///
    /// Because this method takes ownership of `self`, it does not clone the matrix.
    /// If the original matrix must be preserved, clone it before calling this method or use [`Self::echelon_form`].
    ///
    /// # Returns
    ///
    /// A pair `(reduced_matrix, operations)` where:
    ///
    /// - `reduced_matrix` is the row-reduced echelon form of the matrix.
    /// - `operations` is the ordered list of row operations used to obtain it.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let m = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    /// let (rref, ops) = m.into_echelon_form();
    ///
    /// // `m` has been consumed.
    /// // `rref` contains the transformed matrix over GF(2).
    /// ```
    pub fn into_echelon_form(mut self) -> (Self, Vec<(usize, usize)>) {
        let operations = self.echelon_form_in_place();
        (self, operations)
    }

    /// Computes a row-reduced echelon form of the bit-packed matrix over GF(2).
    ///
    /// This method does not modify the original matrix. It clones `self`, performs
    /// elimination on the clone, and returns the transformed matrix together with
    /// the sequence of elementary row operations used during elimination.
    ///
    /// Over GF(2), row addition is implemented as XOR, and every nonzero pivot is
    /// equal to `1`.
    ///
    /// # Row Operations
    ///
    /// The returned operations are ordered and can be applied to the original matrix
    /// to reproduce the returned echelon matrix.
    ///
    /// The following row operations may be recorded:
    ///
    /// - `(i, j)` with `i != j` represents adding row `j` to row `i`:
    ///   `row_i <- row_i + row_j`.
    /// - A row swap between rows `r` and `i` is encoded as three row additions:
    ///   `(r, i)`, `(i, r)`, `(r, i)`.
    ///
    /// This swap encoding is valid over GF(2), where row addition is XOR.
    ///
    /// # Returns
    ///
    /// A pair `(echelon_matrix, operations)` where:
    ///
    /// - `echelon_matrix` is the transformed matrix.
    /// - `operations` is the list of row operations used to obtain it.
    ///
    /// # Notes
    /// - The original matrix is not modified.
    /// - This method allocates a cloned matrix.
    /// - If you want to avoid cloning, use [`Self::into_echelon_form`] or [`Self::echelon_form_in_place`].
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let m = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    ///
    /// let (echelon, ops) = m.echelon_form();
    ///
    /// // `m` is unchanged.
    /// // `echelon` contains the transformed matrix over GF(2).
    /// // `ops` contains the recorded row operations.
    /// ```
    pub fn echelon_form(&self) -> (Self, Vec<(usize, usize)>) {
        self.clone().into_echelon_form()
    }

    fn get_pivot(&self, row: usize) -> Option<usize> {
        for col in 0..self.ncols() {
            if self.get_element(row, col) == 1 {
                return Some(col);
            }
        }

        None
    }

    /// Checks if a GF(2) matrix is in reduced row echelon form (RREF).
    ///
    /// # Returns
    /// `true` if the matrix is in reduced row echelon form; otherwise, `false`.
    pub fn is_reduced_echelon(&self) -> bool {
        let mut previous_pivot: Option<usize> = None;
        let mut seen_zero_row = false;

        for row in 0..self.nrows() {
            match self.get_pivot(row) {
                None => {
                    seen_zero_row = true;
                }
                Some(pivot) => {
                    if seen_zero_row {
                        return false;
                    }

                    if let Some(previous) = previous_pivot {
                        if pivot <= previous {
                            return false;
                        }
                    }

                    previous_pivot = Some(pivot);
                }
            }
        }

        true
    }

    /// Computes the rank of the linear applocation represented by a packed GF(2) matrix which is already in echelon form.
    ///
    /// # Returns
    /// An integer representing the rank of the matrix.
    ///
    pub fn rank_echelon(&self) -> usize {
        self.elements
            .iter()
            .filter(|&&row| row != T::zero())
            .count()
    }

    /// Computes the rank of the linear applocation represented by a GF(2) Bitpacked matrix.
    ///
    /// It first converts the matrix to its RREF before computing the rank.
    ///
    /// # Returns
    /// An integer representing the rank of the matrix.
    ///
    pub fn rank(&self) -> usize {
        let echelon = self.echelon_form();
        echelon.0.rank_echelon()
    }

    /// Computes the rank of the linear applocation represented by a packed GF(2) matrix if the matrix is already in echelon form.
    /// If the matrix is not in reduces echelon form, it returns None.
    /// # Returns
    /// An integer representing the rank of the matrix.
    ///
    pub fn rank_if_echelon(&self) -> Option<usize> {
        if self.is_reduced_echelon() {
            Some(self.rank_echelon())
        } else {
            None
        }
    }

    fn get_value_element(&self, value: T, col: usize) -> u8 {
        Self::get_packed_bit(value, self.ncols(), col)
    }

    fn toggle_value_bit(&self, value: &mut T, col: usize) {
        Self::toggle_packed_bit(value, self.ncols(), col);
    }

    /// Computes a basis of the kernel of a bit-packed GF(2) matrix already in reduced echelon form.
    ///
    /// The returned vectors are also bit-packed using MSB order.
    pub fn kernel_echelon_form(&self) -> Vec<T> {
        assert!(
            self.is_reduced_echelon(),
            "kernel_echelon_form expects the matrix to be in echelon form"
        );

        let cols = self.ncols();

        let mut pivots: Vec<(usize, usize)> = Vec::new(); // (pivot_col, pivot_row)
        let mut is_pivot_col = vec![false; cols];

        for row in 0..self.nrows() {
            if let Some(pivot_col) = self.get_pivot(row) {
                pivots.push((pivot_col, row));
                is_pivot_col[pivot_col] = true;
            }
        }

        let mut kernel_basis: Vec<T> = Vec::new();

        for free_col in 0..cols {
            if is_pivot_col[free_col] {
                continue;
            }

            let mut kernel_vector = T::zero();

            // Set the free variable to 1.
            self.toggle_value_bit(&mut kernel_vector, free_col);

            // Solve pivot variables by back-substitution.
            for &(pivot_col, pivot_row) in pivots.iter().rev() {
                let mut sum = 0u8;

                for col in (pivot_col + 1)..cols {
                    let a = self.get_element(pivot_row, col);
                    let x = self.get_value_element(kernel_vector, col);
                    sum ^= a & x;
                }

                if sum == 1 {
                    self.toggle_value_bit(&mut kernel_vector, pivot_col);
                }
            }

            kernel_basis.push(kernel_vector);
        }

        kernel_basis
    }

    /// Computes a basis of the kernel of any packed GF(2) matrix.
    pub fn kernel(&self) -> Vec<T> {
        let (echelon, _) = self.echelon_form();
        echelon.kernel_echelon_form()
    }

    /// Computes a basis for the image of a bit-packed GF(2) matrix already in
    /// reduced row echelon form.
    ///
    /// The returned basis vectors are the nonzero rows of the matrix, stored in the
    /// same bit-packed row representation used internally.
    ///
    /// # Returns
    ///
    /// A vector of bit-packed rows, each representing a basis vector of the image.
    ///
    /// # Panics
    ///
    /// Panics if the matrix is not in reduced row echelon form.
    ///
    /// # Notes
    ///
    /// - This method does not clone the matrix.
    /// - This method assumes that nonzero rows of the reduced matrix form the image
    ///   basis.
    /// - Use [`Self::image`] if the matrix may not already be reduced.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let m = PackedGF2Matrix::new(vec![0b100u8, 0b010u8, 0b000u8], 3);
    ///
    /// let image_basis = m.image_echelon_form();
    ///
    /// assert_eq!(image_basis, vec![0b100u8, 0b010u8]);
    /// ```
    pub fn image_echelon_form(&self) -> Vec<T> {
        assert!(
            self.is_reduced_echelon(),
            "image_echelon_form expects the matrix to be in reduced echelon form"
        );

        self.elements
            .iter()
            .copied()
            .filter(|&row| row != T::zero())
            .collect()
    }

    /// Computes a basis for the image of the linear map represented by this
    /// bit-packed GF(2) matrix.
    ///
    /// If the matrix is not already in reduced row echelon form, this method first
    /// computes its echelon form. Despite the historical name `echelon_form`, that
    /// method performs Gauss-Jordan-style elimination and returns a row-reduced
    /// echelon form.
    ///
    /// The returned basis vectors are stored in the same bit-packed row
    /// representation used internally by the matrix.
    ///
    /// # Returns
    ///
    /// A vector of bit-packed rows, each representing a basis vector of the image.
    ///
    /// # Notes
    ///
    /// - The original matrix is not modified.
    /// - If the matrix is already in reduced echelon form, no matrix clone is
    ///   needed to extract the image basis.
    /// - If the matrix is not reduced, this method allocates a transformed matrix
    ///   by calling [`Self::echelon_form`].
    /// - Zero rows are ignored.
    /// - Nonzero rows of the reduced matrix form the returned basis.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use lin_algebra::packed_gf2_matrix::PackedGF2Matrix;
    /// let m = PackedGF2Matrix::new(vec![0b110u8, 0b101u8, 0b011u8], 3);
    ///
    /// let image_basis = m.image();
    ///
    /// // Each entry of `image_basis` is a bit-packed row.
    /// ```
    pub fn image(&self) -> Vec<T> {
        if self.is_reduced_echelon() {
            self.image_echelon_form()
        } else {
            let (echelon, _) = self.echelon_form();
            echelon.image_echelon_form()
        }
    }

    /// Applies a sequence of GF(2) row operations to a packed vector.
    ///
    /// The packed vector is interpreted as having length `len`.
    ///
    /// Each operation `(i, j)` represents:
    ///
    /// `v_i <- v_i + v_j`
    ///
    /// Since the entries are in GF(2), addition is XOR.
    fn apply_operations_packed(operations: &[(usize, usize)], mut value: T, len: usize) -> T {
        for &(i, j) in operations {
            if Self::get_packed_bit(value, len, j) == 1 {
                Self::toggle_packed_bit(&mut value, len, i);
            }
        }

        value
    }

    /// Returns `true` if all bits from `start` to `len - 1` are zero.
    fn packed_suffix_is_zero(value: T, len: usize, start: usize) -> bool {
        for idx in start..len {
            if Self::get_packed_bit(value, len, idx) == 1 {
                return false;
            }
        }

        true
    }

    /// Keeps the first `new_len` bits of a packed vector of length `old_len`.
    ///
    /// This is useful after elimination, where the first `self.ncols()` entries
    /// contain the solution and the remaining entries correspond to zero rows.
    fn truncate_packed_prefix(value: T, old_len: usize, new_len: usize) -> T {
        let mut result = T::zero();

        for idx in 0..new_len {
            if Self::get_packed_bit(value, old_len, idx) == 1 {
                Self::toggle_packed_bit(&mut result, new_len, idx);
            }
        }

        result
    }

    /// Returns the column of the bit-packed matrix at index `idx`, packed as a
    /// value of type `T`.
    ///
    /// The returned packed vector has length `self.nrows()`.
    pub fn column_packed(&self, idx: usize) -> T {
        assert!(
            idx < self.ncols(),
            "column index out of bounds: index is {}, but matrix has {} columns",
            idx,
            self.ncols()
        );

        let mut column = T::zero();

        for row in 0..self.nrows() {
            if self.get_element(row, idx) == 1 {
                Self::toggle_packed_bit(&mut column, self.nrows(), row);
            }
        }

        column
    }

    /// Solves the linear system `A * x = b` over GF(2), where `A` is this
    /// bit-packed matrix and `b` is a packed right-hand side vector.
    ///
    /// The right-hand side vector `b` is interpreted as a packed vector of length
    /// `self.nrows()`. The returned solution vector is packed as a value of type
    /// `T` with length `self.ncols()`.
    ///
    /// # Returns
    ///
    /// A packed integer representing the solution vector `x`.
    ///
    /// # Panics
    ///
    /// Panics if:
    ///
    /// - the matrix does not have full column rank.
    /// - the system is inconsistent.
    pub fn solve(&self, b: T) -> T {
        let (echelon, operations) = self.echelon_form();
        let rank = echelon.rank_echelon();

        if rank < self.ncols() {
            panic!("Matrix must have full rank");
        }

        let solved_b = Self::apply_operations_packed(&operations, b, self.nrows());

        if !Self::packed_suffix_is_zero(solved_b, self.nrows(), rank) {
            panic!("Linear system is inconsistent");
        }

        Self::truncate_packed_prefix(solved_b, self.nrows(), self.ncols())
    }

    /// Solves the matrix equation `A * X = Y` over GF(2), where `A` is this
    /// bit-packed matrix and `Y` is a bit-packed right-hand side matrix.
    ///
    /// The returned solution matrix `X` is also bit-packed.
    ///
    /// If `A` has shape `m x n` and `Y` has shape `m x k`, then the returned
    /// matrix has shape `n x k`.
    ///
    /// # Returns
    ///
    /// A bit-packed matrix `X` such that `self * X = y`.
    ///
    /// # Panics
    ///
    /// Panics if:
    ///
    /// - `self.nrows() != y.nrows()`.
    /// - the matrix does not have full column rank.
    /// - the system is inconsistent for at least one column of `y`.
    pub fn solve_matrix_system(&self, y: &PackedGF2Matrix<T>) -> PackedGF2Matrix<T> {
        assert_eq!(
            self.nrows(),
            y.nrows(),
            "left-hand side and right-hand side must have the same number of rows"
        );

        let (echelon, operations) = self.echelon_form();
        let rank = echelon.rank_echelon();

        if rank < self.ncols() {
            panic!("Matrix must have full rank");
        }

        let n_rows = self.ncols(); // rows of X
        let n_cols = y.ncols(); // columns of X

        let mut solution_rows = vec![T::zero(); n_rows];

        for col in 0..n_cols {
            let rhs_col = y.column_packed(col);

            let solved_col = Self::apply_operations_packed(&operations, rhs_col, self.nrows());

            if !Self::packed_suffix_is_zero(solved_col, self.nrows(), rank) {
                panic!("Linear system is inconsistent");
            }

            let solution_col = Self::truncate_packed_prefix(solved_col, self.nrows(), self.ncols());

            for row in 0..n_rows {
                if Self::get_packed_bit(solution_col, n_rows, row) == 1 {
                    Self::toggle_packed_bit(&mut solution_rows[row], n_cols, col);
                }
            }
        }

        PackedGF2Matrix::new(solution_rows, n_cols)
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
