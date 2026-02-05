
struct GF2LinearCombination<'a, T> {
    element: &'a[T],
    idx: usize,
    previous_code: T,
    current_code: T
}


struct GF2VectorSpace<T> {
    elements: Vec<T>
}

trait GF2VectorSpaceTrait<T> {
    type GF2LinearCombinations<'a>: Iterator<Item = T> where Self:'a, T: 'a;
    fn span_iter<'a>(&self) -> Self::GF2LinearCombinations<'a>; 
}
