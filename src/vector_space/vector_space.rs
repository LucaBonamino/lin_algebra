struct LinearCombinations<'a, T> {
    basis: &'a [T],
    mask: T,
    limit: usize,
}


trait VectorSpaceTrait<T> {
    fn base(&self) -> Vec<T>;
    fn span(&self) -> Vec<T>;
}

struct VectorSpace<T> {
    elements: Vec<T>
}


trait VectorSpaceIterTrait<T> {
    type LinearCombinations<'a>: Iterator<Item = T> where Self:'a, T: 'a;
    fn span_iter<'a>(&self) -> Self::LinearCombinations<'a>; 
}
