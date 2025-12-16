use crate::matrix::{MatrixTrait, Matrix};

/// GF(2) matrix.
/// 
/// Implements the trait MatrixTrait: needs to implement
/// - rank
/// - kernel
/// - echelon_form
/// - is_reduced_echelon
/// - image
pub type GF2Matrix = Matrix<u8>;

impl MatrixTrait<u8> for GF2Matrix {
    
    /// Checks if a GF(2) matrix is in reduced row echelon form (RREF).
    ///
    /// # Returns
    /// `true` if the matrix is in reduced row echelon form; otherwise, `false`.
    fn is_reduced_echelon(&self) -> bool{
        let nrows = self.nrows();
        let mut old_piv = 0;
        let mut found_zero_row = false;
        for i in 0..nrows{
            let piv = GF2Matrix::get_pivot(&self.elements[i]);
            match piv {
                None => {
                    found_zero_row = true;
                    continue;
                }
                Some(piv) => {
                    if piv < old_piv || found_zero_row == true {
                        return false;
                    }
                    for j in 0..nrows {
                        if j != i && self.elements[j][piv] != 0 {
                            return false;
                        }
                    }
                    old_piv = piv;

                }
            }
        }
        return true;
    }

    /// Computes the rank of the linear applocation represented by a GF(2) matrix.
    /// 
    /// If the matrix is not in echelon form (RREF),
    /// it first converts the matrix to its RREF before computing the rank.
    /// 
    /// # Returns
    /// An integer representing the rank of the matrix.
    ///
    fn rank(&self) -> usize{
        if self.is_reduced_echelon(){
            return self.rank_echelon_form()
        }
        else {
            let (ech_form, _) = self.echelon_form();
            return ech_form.rank_echelon_form();   
        }

    }

    /// Computes the base of the kernel of the linear applocation represented by a GF(2) matrix.
    ///
    /// If the matrix is not in reduced row echelon form (RREF),
    /// it first converts the matrix to its RREF before computing the kernel.
    ///
    /// # Returns
    /// A vector of row vectors, each representing a basis vector of the kernel.
    fn kernel(&self)-> Vec<Vec<u8>>{
        if self.is_reduced_echelon(){
            println!("{:?}", self.elements);
            return self.kernel_echelon_form();
        }
        else {
            let (ech_form, _) = self.echelon_form();
            println!("{:?}", ech_form);
            return ech_form.kernel_echelon_form();
        }
    }

    /// Computes the reduced echelon form (RREF) of a GF(2) matrix along with the history of all row operations applied.
    /// 
    /// # Returns
    /// A tuple where the first element is the RREF form of the matrix and 
    /// the second is a Vec containing the history of the row operations applied to the matrix to compute the RREF.
    /// Each element of the row operations vector, is a tuple heving the modified row as first element and the row to which it has been summed as second element:
    ///     R1 -> R1 + R2 is represented the entry (R1, R2)
    /// The swap of two rows is represented as 3 entries:
    ///     swap(R1, R2) is represented as (R1, R2), (R2,R1), (R1,R2) 
    fn echelon_form(&self) -> (Self, Vec<(usize, usize)>) {
        let mut m_copy = self.clone();
        let rows = m_copy.nrows();
        let cols = m_copy.ncols();
        let mut operations: Vec<(usize, usize)> = Vec::new();
        let mut lead = 0;

        for r in 0..rows {
            if lead >= cols {
                break;
            }
            let mut i = r;
            while m_copy.elements[i][lead] == 0 {
                i += 1;
                if i == rows {
                    i = r;
                    lead += 1;
                    if lead == cols {
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
            for i in 0..rows {
                if i != r && m_copy.elements[i][lead] == 1 {
                    for j in 0..cols {
                        m_copy.elements[i][j] = (m_copy.elements[i][j] + m_copy.elements[r][j]) % 2;
                    }
                    operations.push((i, r));
                }
            }
            lead += 1;
        }

        (m_copy, operations)
    }

    /// Computes the image of the linear application corresponding to a GF(2) matrix.
    /// 
    /// If the matrix is not in reduced row echelon form (RREF),
    /// it first converts the matrix to its RREF before computing the image.
    /// 
    /// # Returns
    /// A vector of row vectors, each representing a basis vector of the image.
    fn image(&self) -> Vec<Vec<u8>>{
        let mat = if !self.is_reduced_echelon() {
            let (m, _) = self.echelon_form();
            m
        } else {
            self.clone()
        };
        let mut image_base: Vec<Vec<u8>> = Vec::new();
        for i in 0..mat.nrows() {
            let row = mat.elements[i].clone();
            let piv = GF2Matrix::get_pivot(&row);
            if !piv.is_none() {
                image_base.push(row);
            }
        }
        image_base
    }

    /// Returns the pivot index of a row in a GF(2) matrix.
    ///
    /// The pivot is defined as the index of the first non-zero (1) element in the row.
    /// Returns `None` if the row contains only zeros.
    ///
    /// # Arguments
    ///
    /// * `row` - A reference to a vector representing a row in a binary matrix.```
    fn get_pivot(row: &Vec<u8>) -> Option<usize> {
        row.iter().position(|&x| x == 1)
    }

}

impl GF2Matrix  {

    fn apply_operations(operations: &Vec<(usize, usize)>, v: &Vec<u8>) -> Vec<u8> {
        let mut result = v.clone();
        for &(op1, op2) in operations.iter() {
            result[op1] = (result[op1] + result[op2]) % 2;
        }
        result
    }

    /// Returns the column of the matrix by column index.
    fn column(&self, idx: usize) -> Vec<u8> {
        self.elements
            .iter()
            .map(|row| row[idx])
            .collect()
    }

    /// Solves for X such that equation A*X = B where A  and B are GF2Matrix natrices.
    pub fn solve_matrix_system(&self, y: &GF2Matrix) -> GF2Matrix{
        let mut solution: Vec<Vec<u8>> = Vec::new();
        if self.rank() < self.ncols(){
            panic!("Matrix must have full rank");
        }
        let (ech, operations) = self.echelon_form();
        let to_remove = ech.nrows() - ech.rank();
        for i  in 0..y.ncols(){
            let mut solved_b = Self::apply_operations(&operations, &y.column(i));
            for _ in 0..to_remove{
                solved_b.remove(solved_b.len()-1);
            }
            solution.push(solved_b);
        }

        GF2Matrix::new(solution)

    }

    /// Solves for X such that equation A*x = b where A  is a GF2Matrix and b a Vec<u8>.
    pub fn solve(&self, b: &Vec<u8>) -> Vec<u8>{
        // et mut y_copy = y.clone();
        if self.rank() < self.ncols(){
            panic!("Matrix must have full rank");
        }
        let (ech, operations) = self.echelon_form();
        let to_remove = ech.nrows() - ech.rank();
        let mut solved_b= Self::apply_operations(&operations, b);
        for _ in 0..to_remove{
            solved_b.remove(solved_b.len()-1);
        }
        solved_b

    }


    fn swap_rows(&mut self, row1: usize, row2: usize) {
        self.elements.swap(row1, row2);
    }

    fn rank_echelon_form(&self) -> usize {
        let mut count = 0;
        let mut pivot_columns = std::collections::HashSet::new();

        for i in 0..self.nrows() {
            let p = GF2Matrix::get_pivot(&self.elements[i]);
            if let Some(col) = p {
                if pivot_columns.insert(col) {
                    count += 1;
                }
            }
        }
        count
    }


    fn kernel_echelon_form(&self) -> Vec<Vec<u8>> {
        let rows = self.nrows();
        let cols = self.ncols();

        let mut pivots = std::collections::HashMap::new();
        let mut kernel_base: Vec<Vec<u8>> = Vec::new();
        let mut free_columns: Vec<usize> = Vec::new();
        let mut row_index = 0;

        for j in 0..cols {
            if row_index < rows && self.elements[row_index][j] == 1 {
                pivots.insert(j, row_index);
                row_index += 1;
            } else {
                free_columns.push(j);
            }
        }

        for &free_col in &free_columns {
            let mut kernel_vector = vec![0; cols];
            kernel_vector[free_col] = 1;

            for (&p_index, &p_row) in &pivots {
                let mut sum = 0;

                for col in (0..cols).rev() {
                    if col != p_index {
                        sum = sum ^ (self.elements[p_row][col] * kernel_vector[col]);
                    }
                }

                kernel_vector[p_index] = sum;
            }

            kernel_base.push(kernel_vector);
        }

        kernel_base
    }
    
}