#![feature(trait_alias)]

use std::ops::*;
pub use num_traits::identities::{One, Zero};
pub use ramp::Int;
pub use ramp::rational::Rational;

use static_assertions::assert_impl_all;

pub trait AssociativeAddition {}
pub trait CommutativeAddition {}
pub trait AdditiveIdentity = Zero;

pub trait ClosedAddition = where
    Self: Sized,
    Self: Add<Self, Output=Self>,
    for<'a> &'a Self: Add<Self, Output=Self>,
    for<'a> Self: Add<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Add<&'b Self, Output=Self>,
    Self: AddAssign<Self>,
    for<'a> Self: AddAssign<&'a Self>;

pub trait ClosedSubtraction = where
    Self: Sized,
    Self: Sub<Self, Output=Self>,
    for<'a> &'a Self: Sub<Self, Output=Self>,
    for<'a> Self: Sub<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Sub<&'b Self, Output=Self>,
    Self: SubAssign<Self>,
    for<'a> Self: SubAssign<&'a Self>;

pub trait AdditiveInverse = where
    Self: Sized,
    Self: ClosedSubtraction,
    Self: Neg<Output=Self>,
    for<'a> &'a Self: Neg<Output=Self>;

pub trait Magma = ClosedAddition;
pub trait Semigroup = AssociativeAddition + Magma;
pub trait Monoid = Semigroup + AdditiveIdentity;
pub trait UnitalMagma = Magma + AdditiveIdentity;
pub trait CommutativeMonoid = Monoid + CommutativeAddition;
pub trait Group = Monoid + AdditiveInverse;
pub trait Loop = Magma + AdditiveIdentity + AdditiveInverse;
pub trait AbelianGroup = Group + CommutativeAddition;

// Second operation

pub trait AssociativeMultiplication {}
pub trait CommutativeMultiplication {}
pub trait NonZeroProdProperty {}
pub trait MultiplicativeIdentity = One;

pub trait ClosedMultiplication = where
    Self: Sized,
    Self: Mul<Self, Output=Self>,
    for<'a> &'a Self: Mul<Self, Output=Self>,
    for<'a> Self: Mul<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Mul<&'b Self, Output=Self>,
    Self: MulAssign<Self>,
    for<'a> Self: MulAssign<&'a Self>;

pub trait ClosedDivision = where
    Self: Sized,
    Self: Div<Self, Output=Self>,
    for<'a> &'a Self: Div<Self, Output=Self>,
    for<'a> Self: Div<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Div<&'b Self, Output=Self>,
    Self: DivAssign<Self>,
    for<'a> Self: DivAssign<&'a Self>;

pub trait ExactDivision {}

pub trait Semiring = CommutativeMonoid + ClosedMultiplication + AssociativeMultiplication;
pub trait NonassociativeRing = AbelianGroup + ClosedMultiplication;
pub trait Rng = NonassociativeRing + AssociativeMultiplication;
pub trait Ring = Rng + MultiplicativeIdentity;
pub trait CommutativeRing = Ring + CommutativeMultiplication;
pub trait IntegralDomain = CommutativeRing + NonZeroProdProperty;
pub trait Field = IntegralDomain + ClosedDivision + ExactDivision;

// Implement our traits for integral types

macro_rules! impl_num_type {
    ($t:ty) => {
        impl AssociativeAddition for $t {}
        impl CommutativeAddition for $t {}
        impl AssociativeMultiplication for $t {}
        impl CommutativeMultiplication for $t {}
        impl NonZeroProdProperty for $t {}
    };
}

impl_num_type!(u8);
impl_num_type!(u16);
impl_num_type!(u32);
impl_num_type!(u64);
impl_num_type!(u128);
impl_num_type!(i8);
impl_num_type!(i16);
impl_num_type!(i32);
impl_num_type!(i64);
impl_num_type!(i128);
impl_num_type!(f32);
impl_num_type!(f64);

// Implementations for ramp's arbitrary sized ints and rationals
impl_num_type!(Int);
impl_num_type!(Rational);

impl ExactDivision for f32 {}
impl ExactDivision for f64 {}
impl ExactDivision for Rational {}

// Assertions

assert_impl_all!(f32: Field);
assert_impl_all!(f64: Field);
assert_impl_all!(Int: Ring);
assert_impl_all!(Rational: Field);


#[cfg(test)]
mod tests {
    use super::*;

    fn my_append<T: Semigroup>(a: &T, b: &T) -> T {
        a + b
    }

    #[test]
    fn main_test() {
        assert_eq!(my_append(&10, &12), 22);
    }
}
