# lin_algebra
A Rust library for linear algebra.

## Features
All operations are currently implemented for matrices over GF(2):

- Compute the echelon form of a matrix along with the history of all the row operations applied.
- Compute the rank of the linear application represented by a matrix.
- Compute the kernel of the linear application represented by a matrix.
- Compute the image of the linear application represented by a matrix.
- Solve system of equations in GF(2).

## Installation

Add the dependency to your `Cargo.toml`:

```toml
[dependencies]
lin_algebra = "0.2.0"
```

## Usage
```Rust
use lin_algebra::gf2_matrix::GF2Matrix;
use lin_algebra::matrix::MatrixTrait;

fn main() {
    // Reduces echelon form

    let gf2_mat = gf2_matrix::GF2Matrix::new(
        vec![vec![1,0,0,0], vec![0,1,0,1], vec![0,1,0,1]]
    )
    let (educed_echelon_form, row_operations) = gf2_mat.echelon_form();

    println!("Reduces echelon form {:?}", educed_echelon_form);
    println!(
        "Row operation applied to reach the reduces echelon form {:?}", row_operations
    );

    // Kernel
    println!("Kenel {:?}", gf2_mat.kernel());

    // Image
    println!("Image {:?}", gf2_mat.image());

    // Rank
    println!("Rank {}", gf2_mat.rank());

    // Solve M x = b
    let b = vec![1,0,0,1]
    println!("x: {:?}", gf2_mat.solve(&b));

    // Solve A X = Y
    let y_matrix = gf2_matrix::GF2Matrix::new(
        vec![vec![1,0,0,1], vec![1,1,0,1], vec![0,1,0,1]]
    )
    println!("x: {:?}", gf2_mat.solve_matrix_system(&y_matrix));

}
```

## License
MIT