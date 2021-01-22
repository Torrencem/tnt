//! This crate defines many algebraic traits and trait aliases that are useful for generic
//! mathematical applications. Most of the marker trait aliases are automatically implemented from
//! combinations of std traits. In particular, ClosedAddition is implemented for T where all
//! combinations of T and &T can be added together (in addition to AddAssign), and there are
//! similar traits for subtraction, multiplication and division.

#![feature(trait_alias)]

use std::ops::*;
use ordered_float::NotNan;
pub use num_traits::identities::{One, Zero};
pub use ramp::Int;
pub use ramp::rational::Rational;

use static_assertions::assert_impl_all;

/// A set whose addition operation is associative. This is a *free trait*.
pub trait AssociativeAddition {}
/// A set whose addition operation is commutative. This is a *free trait*.
pub trait CommutativeAddition {}
/// A set with an additive identity. This is an alias for the Zero trait, which should be
/// implemented instead.
pub trait AdditiveIdentity = Zero;

/// A set whose elements can be added together.
pub trait ClosedAddition = where
    Self: Sized,
    Self: Add<Self, Output=Self>,
    for<'a> &'a Self: Add<Self, Output=Self>,
    for<'a> Self: Add<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Add<&'b Self, Output=Self>,
    Self: AddAssign<Self>,
    for<'a> Self: AddAssign<&'a Self>;

/// A set whose elements can be subtracted.
pub trait ClosedSubtraction = where
    Self: Sized,
    Self: Sub<Self, Output=Self>,
    for<'a> &'a Self: Sub<Self, Output=Self>,
    for<'a> Self: Sub<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Sub<&'b Self, Output=Self>,
    Self: SubAssign<Self>,
    for<'a> Self: SubAssign<&'a Self>;

/// A set whose elements all have additive inverses.
pub trait AdditiveInverse = where
    Self: Sized,
    Self: ClosedSubtraction,
    Self: Neg<Output=Self>,
    for<'a> &'a Self: Neg<Output=Self>;

/// A set whose elements can be added, with no other restrictions.
pub trait Magma = ClosedAddition;
/// A Magma whose addition is associative.
pub trait Semigroup = AssociativeAddition + Magma;
/// A Semigroup with an additive identity element.
pub trait Monoid = Semigroup + AdditiveIdentity;
/// A Magma with an additive identity element.
pub trait UnitalMagma = Magma + AdditiveIdentity;
/// A Monoid whose addition is commutative.
pub trait CommutativeMonoid = Monoid + CommutativeAddition;
/// A Monoid where each element has an additive inverse.
pub trait Group = Monoid + AdditiveInverse;
/// A Group without associative addition.
pub trait Loop = Magma + AdditiveIdentity + AdditiveInverse;
/// A Group whose addition is commutative.
pub trait AbelianGroup = Group + CommutativeAddition;

// Second operation

/// A set whose multiplication operation is associative. This is a *free trait*
pub trait AssociativeMultiplication {}
/// A set whose multiplication operation is commutative. This is a *free trait*
pub trait CommutativeMultiplication {}
/// A set with multiplication with the property that for all elements a, b, ab is 0 if and only if
/// one of a or b was 0.
pub trait NonZeroProdProperty {}
/// A set with a multiplicative identity. This is an alias for the One trait, which should be
/// implemented instead.
pub trait MultiplicativeIdentity = One;

/// A set whose elements can be multiplied together.
pub trait ClosedMultiplication = where
    Self: Sized,
    Self: Mul<Self, Output=Self>,
    for<'a> &'a Self: Mul<Self, Output=Self>,
    for<'a> Self: Mul<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Mul<&'b Self, Output=Self>,
    Self: MulAssign<Self>,
    for<'a> Self: MulAssign<&'a Self>;

/// A set whose elements can be divided. This only shows that there is some "division-like
/// operation" on the set. To guarantee more, see the ExactDivision trait and the DivRem trait.
pub trait ClosedDivision = where
    Self: Sized,
    Self: Div<Self, Output=Self>,
    for<'a> &'a Self: Div<Self, Output=Self>,
    for<'a> Self: Div<&'a Self, Output=Self>,
    for<'a, 'b> &'a Self: Div<&'b Self, Output=Self>,
    Self: DivAssign<Self>,
    for<'a> Self: DivAssign<&'a Self>;

/// A set whose elements can be divided exactly, so that for all elements a and b, (a / b) * b ==
/// a. Note that this specifically is not how integer division works, for example. Hence, integer
/// types implement ClosedDivision, but not ExactDivision.
pub trait ExactDivision {}

/// A set that is a commutative monoid additively, with an associative multiplication operation.
pub trait Semiring = CommutativeMonoid + ClosedMultiplication + AssociativeMultiplication;
/// An additive group with a multiplication operation defined on it that isn't necessarily
/// multiplicatively associative.
pub trait NonassociativeRing = AbelianGroup + ClosedMultiplication;
/// A NonassociativeRing whose multiplication is associative.
pub trait Rng = NonassociativeRing + AssociativeMultiplication;
/// A Rng with a multiplicative identity.
pub trait Ring = Rng + MultiplicativeIdentity;
/// A Ring whose multiplication operation is commutative.
pub trait CommutativeRing = Ring + CommutativeMultiplication;
/// A CommutativeRing with the non zero product property; that is, for all elements a and b, a*b is
/// 0 if and only if one of a or b is zero.
pub trait IntegralDomain = CommutativeRing + NonZeroProdProperty;
/// An integral domain whose elements all have multiplicative inverses (so that exact division can
/// be defined).
pub trait Field = IntegralDomain + ClosedDivision + ExactDivision;

/// A set on which a Euclidean Function is defined. A Euclidean Function is a norm that can be used
/// to write a Euclidean algorithm for computing GCD's. This trait is useful for getting a
/// free implementation of GCD for integer-like types.
pub trait EuclideanFunction {
    type Order: std::cmp::Ord;

    fn norm(&self) -> Self::Order;
}

/// A set which has a GCD function. For integer types, this is automatically implemented if a
/// Euclidean function and a DivRem function are defined.
pub trait Gcd: EuclideanFunction {
    fn gcd(&self, other: &Self) -> Self;
}

impl<T> Gcd for T
where T: EuclideanFunction + DivRem<Output=T> + Clone + Zero {
    fn gcd(&self, other: &Self) -> Self {
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

/// A convenience function for computing the GCD of two things that implement the Gcd trait.
pub fn gcd<T: Gcd>(a: &T, b: &T) -> T {
    a.gcd(b)
}

/// An Integral domain with a Euclidean function defined on it.
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

/// Returned from the divrem function
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub struct DivRemResult<T> {
    pub div: T,
    pub rem: T,
}

/// Returned from the pseudo_divrem function
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub struct PseudoDivRemResult<T, U> {
    pub mul: U,
    pub div: T,
    pub rem: T,
}

pub trait PseudoDivRem {
    type Output;
    type MultType;
    
    /// Given two elements a and b, pseudo_divrem(a, b) computes values m, d, and r, such that
    /// ma = bd + r. This is a "pseudo-division" of a by b. It's possible to constrain m to be
    /// some other type of value. In particular, pseudo_divrem(a, b).mul is m, pseudo_divrem(a,
    /// b).div is d, and pseudo_divrem(a, b).rem is r. For Polynomials, this is useful for asserting
    /// that m is a constant polynomial.
    fn pseudo_divrem(a: Self, b: Self) -> PseudoDivRemResult<Self::Output, Self::MultType>;
}

/// A convenience function for performing exact division when possible. If you have two elements a
/// and b such that a = b*k, pdiv(b, a) should always return k.
pub fn pdiv<U, V, T: PseudoDivRem<Output=U, MultType=V>>(a: T, b: T) -> U {
    PseudoDivRem::pseudo_divrem(a, b).div
}

pub trait DivRem {
    type Output;

    /// Given two elements a and b, divrem(a, b) computes values d and r such that a = b*d + r. In
    /// particular, divrem(a, b).div is d, and divrem(a, b).rem is r.
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

/// A set whose elements can be taken to powers. As long as MulAssign, Mul, One and Clone are
/// defined for a type, Power will be automatically implemented with an exponentiation by squaring
/// algorithm.
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

        // Unrolled first iteration so that V and T work out.
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
