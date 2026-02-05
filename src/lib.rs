pub mod gf2_matrix;
pub mod matrix;
use gf2_matrix::GF2Matrix;
mod vector_space; 

#[cfg(test)]
mod tests {
    use crate::matrix::MatrixTrait;

    use super::*;

    #[test]
    fn instantiate_matrix() {
        let mat = GF2Matrix::new(vec![vec![1, 0], vec![1, 1]]);
        let el = mat.elements;
        assert_eq!(el, vec![vec![1, 0], vec![1, 1]]);
    }

    #[test]
    fn echelon_form_matrix() {
        let mat = GF2Matrix::new(vec![vec![1, 0], vec![1, 1]]);
        let (ech_form, _) = mat.echelon_form();
        assert_eq!(ech_form.elements, vec![vec![1, 0], vec![0, 1]]);
        let mat = GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1], vec![0, 1, 0, 1]]);
        let (ech_form, _) = mat.echelon_form();
        assert_eq!(
            ech_form.elements,
            vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1], vec![0, 0, 0, 0]]
        );
    }

    #[test]
    fn echelon_form_ops() {
        let mat = GF2Matrix::new(vec![vec![1, 0], vec![1, 1]]);
        let (_, ops) = mat.echelon_form();
        assert_eq!(ops, vec![(1, 0)])
    }

    #[test]
    fn ncols_nrows() {
        let mat = GF2Matrix::new(vec![vec![1, 0, 1], vec![1, 1, 0]]);
        assert_eq!(mat.nrows(), 2);
        assert_eq!(mat.ncols(), 3)
    }

    #[test]
    fn get_pivot() {
        assert_eq!(GF2Matrix::get_pivot(&vec![1, 0, 0, 1]).unwrap(), 0);
        assert_eq!(GF2Matrix::get_pivot(&vec![0, 1, 1, 0]).unwrap(), 1);
        assert_eq!(GF2Matrix::get_pivot(&vec![0, 0, 1, 1]).unwrap(), 2);
        assert_eq!(GF2Matrix::get_pivot(&vec![0, 0, 0, 1]).unwrap(), 3);
        assert!(GF2Matrix::get_pivot(&vec![0, 0, 0, 0]).is_none());
    }

    #[test]
    fn is_reduced_echelon() {
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![1, 1, 0, 1]]).is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![1, 1, 0, 1]]).is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 1, 0, 1], vec![0, 1, 0, 1]]).is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 1], vec![0, 1, 0, 1]]).is_reduced_echelon(),
            true
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![1, 1, 0, 1], vec![0, 0, 1, 1]])
                .is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 0], vec![0, 0, 1, 1]])
                .is_reduced_echelon(),
            true
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1], vec![0, 1, 0, 1]])
                .is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 0, 0, 0], vec![0, 1, 0, 1]])
                .is_reduced_echelon(),
            false
        );
    }

    #[test]
    fn is_reduced_echelon_zero_row() {
        assert_eq!(
            GF2Matrix::new(vec![vec![0, 0, 0, 0], vec![1, 0, 0, 0], vec![0, 1, 0, 1]])
                .is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 0, 0, 0], vec![0, 1, 0, 1]])
                .is_reduced_echelon(),
            false
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 1, 0], vec![0, 0, 0, 0]]).is_reduced_echelon(),
            true
        );
    }

    #[test]
    fn rank() {
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1]]).rank(),
            2
        );
        assert_eq!(
            GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![1, 0, 0, 0]]).rank(),
            1
        );
    }

    #[test]
    fn image() {
        let mat_1 = GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1]]);
        assert_eq!(mat_1.image(), mat_1.elements);
        let mat_2 = GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1], vec![0, 1, 0, 1]]);
        assert_eq!(mat_2.image(), mat_1.elements);
        let mat = GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 0, 0, 0]]);
        assert_eq!(mat.image(), vec![vec![1, 0, 0, 0]]);
    }

    #[test]
    fn kernel() {
        let mat = GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 1, 0, 1]]);
        assert_eq!(mat.kernel(), vec![vec![0, 0, 1, 0], vec![0, 1, 0, 1]]);
        let mat = GF2Matrix::new(vec![vec![1, 0, 0, 0], vec![0, 0, 0, 0], vec![0, 1, 0, 1]]);
        assert_eq!(mat.kernel(), vec![vec![0, 0, 1, 0], vec![0, 1, 0, 1]])
    }

    #[test]
    fn test_solve_matrix_system() {
        let left_elements = vec![vec![1, 0, 0], vec![0, 1, 1], vec![1, 0, 1]];

        let right_elements = vec![vec![0, 0, 1], vec![0, 1, 1], vec![1, 1, 1]];

        let m = GF2Matrix::new(left_elements);
        let rm = GF2Matrix::new(right_elements);
        let r = m.solve_matrix_system(&rm);
        assert_eq!(
            r.elements,
            vec![vec![0, 0, 1], vec![1, 0, 1], vec![1, 1, 0]]
        );
    }

    #[test]
    fn test_solve_matrix_system_linear_dependence() {
        let left_elements = vec![vec![1, 0, 0], vec![0, 1, 1], vec![1, 0, 1], vec![1, 1, 1]];

        let right_elements = vec![vec![0, 0, 1], vec![0, 1, 1], vec![1, 1, 1], vec![1, 1, 0]];

        let m = GF2Matrix::new(left_elements);
        let rm = GF2Matrix::new(right_elements);
        let r = m.solve_matrix_system(&rm);
        assert_eq!(
            r.elements,
            vec![vec![0, 0, 1], vec![1, 0, 1], vec![1, 1, 0]]
        );
    }

    #[test]
    #[should_panic(expected = "Matrix must have full rank")]
    fn test_solve_matrix_system_no_full_rank() {
        let left_elements = vec![vec![1, 0, 0], vec![0, 1, 1]];

        let right_elements = vec![vec![0, 0, 1], vec![0, 1, 1], vec![1, 1, 1]];

        let m = GF2Matrix::new(left_elements);
        let rm = GF2Matrix::new(right_elements);
        m.solve_matrix_system(&rm);
    }

    #[test]
    fn test_solve() {
        let left_elements = vec![vec![1, 0, 0], vec![0, 1, 1], vec![1, 0, 1]];

        let b = vec![0, 0, 1];

        let m = GF2Matrix::new(left_elements);
        let r = m.solve(&b);
        assert_eq!(r, vec![0, 1, 1]);
    }

    #[test]
    #[should_panic(expected = "Matrix must have full rank")]
    fn test_solve_no_full_rank() {
        let left_elements = vec![vec![1, 0, 0], vec![0, 1, 1]];

        let b = vec![0, 0, 1];

        let m = GF2Matrix::new(left_elements);
        m.solve(&b);
    }
}
