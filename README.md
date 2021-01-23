# TNT

tnt is a set of Rust crates for algebraic number theory. It provides traits for working with various structures, as well as types for working with polynomials, finite fields, and number fields. tnt requires nightly for `feature(trait_alias)` until `trait_alias` gets stabilized.

## tnt-alg

The base crate `tnt-alg` provides many basic traits and trait aliases for abstract algebra structures. For example, to write a function that's generic over the type of Ring that's used in an algorithm, you could write:

```rust
use tnt_alg::Ring;

fn my_func<T: Ring + std::fmt::Display>(a: &T, b: &T) {
    println!("{} + {} = {}", a, b, a + b);
    println!("{} * {} = {}", a, b, a * b);
}
```

You'll notice that unlike in similar crates, you can add and multiply references to `T` when `T: Ring`. This is because this crate makes use of `trait_alias`, which as it stands is much more capable than typical traits. In particular, most of the traits in `tnt_alg` are just aliases for all the combinations of being able to add `T` and its references. This means that to implement `Ring` or the other algebraic traits, the majority of the work is in implementing all the `std::ops` traits, in addition to implementing some marker traits that indicate that certain operations are associative and commutative.

In addition to this, `tnt-alg` re-exposes the `Int` and `Rational` types from the `ramp` crate, and implements the appropriate marker traits for these types.

## tnt-polynomial

The crate `tnt-polynomial` has the `Polynomial<T>` type, which represents a polynomial with coefficients of type `T`. It provides most basic manipulations, as well as algorithms for dividing, computing remainders, and computing gcds of polynomials.

## tnt-number-field

The crate `tnt-number-field` has the `AlgebraicNumber<T>` and `AlgebraicNumberR<T, R>` types. The most basic use case is just using the `AlgebraicNumber<T>` type for some integral `T` type (i.e. `i64` or `Int` from `tnt-alg`.

The `AlgebraicNumberR<T, R>` type is generic over the type of its reference to its modulus. This is to allow for potentially more efficient operations (for example, `AlgebraicNumberR<T, Polynomial<T>>` is more time-efficient, but less space efficient than `AlgebraicNumberR<T, Arc<Polynomial<T>>>`, which is the default used for `AlgebraicNumber<T>`).
