#![feature(trait_alias)]

use std::ops::*;
use ordered_float::NotNan;
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

pub trait EuclideanFunction {
    type Order: std::cmp::Ord;

    fn norm(&self) -> Self::Order;
    fn gcd(&self, other: &Self) -> Self 
    where Self: DivRem<Output=Self> + Clone + Zero
    {
        let mut a = self.clone();
        let mut b = other.clone();
        while !b.is_zero() {
            let t = b.clone();
            b = DivRem::divrem(a, b).rem;
            a = t;
        }
        a
    }
}

pub fn gcd<T: EuclideanFunction>(a: &T, b: &T) -> T
    where T: DivRem<Output=T> + Clone + Zero {
    a.gcd(b)
}

pub trait EuclideanDomain = EuclideanFunction + IntegralDomain;

// Implement our traits for integral types

macro_rules! impl_num_type {
    ($t:ty) => {
        impl AssociativeAddition for $t {}
        impl CommutativeAddition for $t {}
        impl AssociativeMultiplication for $t {}
        impl CommutativeMultiplication for $t {}
        impl NonZeroProdProperty for $t {}

        impl EuclideanFunction for $t {
            type Order = $t;

            fn norm(&self) -> Self::Order {
                self.clone()
            }
        }
    };
}

macro_rules! impl_num_type_f {
    ($t:ty) => {
        impl AssociativeAddition for $t {}
        impl CommutativeAddition for $t {}
        impl AssociativeMultiplication for $t {}
        impl CommutativeMultiplication for $t {}
        impl NonZeroProdProperty for $t {}

        impl EuclideanFunction for $t {
            type Order = NotNan<$t>;

            fn norm(&self) -> Self::Order {
                NotNan::new(self.clone()).unwrap()
            }
        }
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
impl_num_type_f!(f32);
impl_num_type_f!(f64);

// Implementations for ramp's arbitrary sized ints and rationals
impl_num_type!(Int);
impl_num_type!(Rational);

impl ExactDivision for f32 {}
impl ExactDivision for f64 {}
impl ExactDivision for Rational {}

// Assertions

assert_impl_all!(f32: Field, EuclideanFunction);
assert_impl_all!(f64: Field, EuclideanFunction);
assert_impl_all!(Int: Ring, EuclideanFunction);
assert_impl_all!(Rational: Field, EuclideanFunction);

// DivRem

#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub struct DivRemResult<T> {
    pub div: T,
    pub rem: T,
}

#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub struct PseudoDivRemResult<T, U> {
    pub mul: U,
    pub div: T,
    pub rem: T,
}

pub trait PseudoDivRem {
    type Output;
    type MultType;

    fn pseudo_divrem(a: Self, b: Self) -> PseudoDivRemResult<Self::Output, Self::MultType>;
}

pub trait DivRem {
    type Output;

    fn divrem(a: Self, b: Self) -> DivRemResult<Self::Output>;
}

impl<V, T: Rem<T, Output=V> + Div<T, Output=V> + Clone> DivRem for T {
    type Output=V;

    fn divrem(a: T, b: T) -> DivRemResult<V> {
        DivRemResult {
            div: a.clone() / b.clone(),
            rem: a % b,
        }
    }
}

impl<V, T> PseudoDivRem for T
where T: DivRem<Output=V>,
      V: One {
    type Output = V;
    type MultType = V;

    fn pseudo_divrem(a: Self, b: Self) -> PseudoDivRemResult<V, V> {
        let dr = DivRem::divrem(a, b);
        PseudoDivRemResult {
            mul: One::one(),
            div: dr.div,
            rem: dr.rem,
        }
    }
}

assert_impl_all!(i64: DivRem);
assert_impl_all!(f64: DivRem);
assert_impl_all!(&u32: DivRem);
assert_impl_all!(Int: DivRem);
assert_impl_all!(&Int: DivRem);

// Multiplicative Exponentiation

pub trait Power {
    type Power;
    type Output;

    fn pow(self, val: &Self::Power) -> Self::Output;
}

impl<V: One + MulAssign<T> + MulAssign<V> + Clone, T: Mul<T, Output=V> + Clone> Power for T
where for<'a> V: MulAssign<&'a V>
{
    type Power = u64;
    type Output = V;

    fn pow(self, val: &u64) -> V {
        // A simple exponentiation by squaring algorithm. This is a little complicated because we
        // are general enough that V can be u32 and T can be &u32, for example
        let mut res: V = One::one();
        let mut e = *val;

        // Unrolled first iteration
        if e & 1 == 1 {
            res *= self.clone();
        }
        e >>= 1;
        
        let mut acc: V = self.clone() * self;

        while e != 0 {
            if e & 1 == 1 {
                res *= &acc;
            }
            acc *= acc.clone();
            e >>= 1;
        }
        res
    }
}

assert_impl_all!(u32: Power);
assert_impl_all!(&u32: Power);
assert_impl_all!(Int: Power);
assert_impl_all!(&Int: Power);
assert_impl_all!(f32: Power);

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

    #[test]
    fn test_exponent() {
        assert_eq!(2.pow(&4), 16);
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(&150, &27), 3);
        assert_eq!(gcd(&150, &100), 50);
        assert_eq!(gcd(&27, &25), 1);
    }
}
