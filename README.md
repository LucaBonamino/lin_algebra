# lin_algebra
A Rust library for linear algebra.

## Features
Operations are currently implemented for matrices over GF(2), with both
explicit and bit-packed representations:

- `GF2Matrix`: a standard explicit matrix representation over GF(2).
- `PackedGF2Matrix`: a compact bit-packed representation for improved memory usage
  and efficient bitwise operations.

Supported operations include:

- Compute the echelon form of a matrix along with the history of all row operations applied.
- Compute the rank of the linear application represented by a matrix.
- Compute the kernel of the linear application represented by a matrix.
- Compute the image of the linear application represented by a matrix.
- Solve systems of equations in GF(2).
- Convert between packed and explicit GF(2) matrix representations.
- Multiply a bit-packed matrix by a bit-packed vector.

## Installation

Add the dependency to your `Cargo.toml`:

```toml
[dependencies]
lin_algebra = "0.3.1"
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
    println!("Kernel {:?}", gf2_mat.kernel());

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
```rust
use lin_algebra::packed_gf2_matrix::{BitOrder, PackedGF2Matrix};

fn main() {
    // Each integer encodes one row of the matrix in bit-packed form.
    let packed = PackedGF2Matrix::<u8>::new(
        vec![
            0b1000u8,
            0b0101u8,
            0b0101u8,
        ],
        4,
    );

    let (echelon, ops) = packed.echelon_form();

    println!("Packed echelon form: {:?}", echelon);
    println!("Row operations: {:?}", ops);

    // Convert to an explicit GF(2) matrix
    let explicit = packed.from_int_matrix_to_gf2_matrix(BitOrder::LSB);
    println!("Explicit matrix: {:?}", explicit);

    // Example of packed matrix-vector multiplication
    let v = 0b0011u8;
    let result = PackedGF2Matrix::matrix_by_vector(&packed, &v);
    println!("Packed matrix-vector product: {:?}", result);
}
```
## Matrix Representations

This crate currently provides two matrix representations over GF(2):

- `GF2Matrix`: stores entries explicitly as `0` and `1`.
  This is easier to inspect and manipulate directly.
- `PackedGF2Matrix<T>`: stores each row as a packed unsigned integer type
  such as `u8`, `u16`, `u32`, `u64`, or `u128`.
  This is more compact and allows efficient XOR-based row operations.

Use `GF2Matrix` when clarity is more important.
Use `PackedGF2Matrix` when performance or memory efficiency matters.

## License
MIT

## Python bindings

Python bindings are provided in the separate project
[`gf2_lin_algebra`](https://github.com/LucaBonamino/gf2_lin_algebra).

They are intentionally limited to GF(2) and focus on performance and simplicity. The Rust crate `lin_algebra` is designed to be more general.

## Roadmap

This project is under active development, with a focus on both research and practical use.

Current and planned directions include:
- improving the internal design of GF(2) matrices
- adding functionality for optimized bit-packed matrix representations
- generalizing the Rust core to support linear algebra over arbitrary fields

The Python bindings (`gf2_lin_algebra`) are intentionally limited to GF(2), while the Rust crate (`lin_algebra`) aims to be more general.

For a detailed list of planned improvements and areas where help is welcome, see [ROADMAP.md](ROADMAP.md).

## Contributing

Contributions are welcome!

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.