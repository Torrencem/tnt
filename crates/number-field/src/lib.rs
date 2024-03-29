//! This module defines the AlgebraicNumber type. Note: Be slightly careful using the One and Zero
//! trait implementations of AlgebraicNumber, since they're pretty ill-defined. See the code for
//! more details.

// Right now, this module doesn't use the r_i trick described on page 159-160 of the course book
// for multiplication, and so it could serve to be a little more time efficient.
// NOTE: Be careful using One::one() and Zero::zero(), since they are pretty ill-defined (no
// modulus). If you know what modulus you want you should be using that to avoid logic errors

use tnt_polynomial::*;
use tnt_alg::*;
use std::ops::Deref;
use std::sync::Arc;

use static_assertions::assert_impl_all;

/// A module containing the raw types that are generic over the container of their reference to
/// their fields.
pub mod raw {
    use super::*;
    use std::ops::*;

    /// An element of a number field isomorphic to Frac(T)\[x\]/<*modulus>
    /// This struct is generic over the way to references the modulus.
    #[derive(Clone, Debug)]
    pub struct AlgebraicNumberR<T, R> 
        // TODO: Make this bound a trait alias? Might be cool.
    where T: Ring + Gcd + PartialEq + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone,
    {
        pub(crate) value: Polynomial<T>,
        // Invariant: denom should never be negative!
        pub(crate) denom: T,
        // This value is 0 if this AlgebraicNumberR was generated from the One or Zero traits. If
        // that's the case, traits will use the "more defined" moduli in operations.
        pub(crate) modulus: R,
    }

    use std::fmt;
    impl<T, R> fmt::Display for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq + Clone + fmt::Display + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone,
    {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            let should_write_parens = self.value.coeffs().iter().filter(|x| !x.is_zero()).count() > 1;
            if should_write_parens {
                write!(f, "(")?;
            }
            self.value.format_with_var("a", false, f)?;
            if should_write_parens {
                write!(f, ")")?;
            }
            if !self.denom.is_one() {
                write!(f, " / {}", self.denom)?;
            }
            Ok(())
        }
    }

    impl<T, R> PartialEq for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone,
    {
        fn eq(&self, other: &Self) -> bool {
            self.value == other.value && self.denom == other.denom
        }
    }

    impl<T, R> Eq for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone,
    { }

    impl<T, R, U> AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq + PseudoDivRem<Output=T, MultType=U> + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone,
    {
        /// Create a new `AlgebraicNumber`. Use this method if you already have an R that you can
        /// pass in. This makes it possible for multiple `AlgebraicNumber`s to share a reference to
        /// the same modulus.
        ///
        /// Note that this function does not check if the modulus is irreducible, and so logic
        /// errors may occur if that is not the case.
        pub fn new_with_modulus(value: Polynomial<T>, denom: T, modulus: R) -> Self {
            AlgebraicNumberR {
                value, denom, modulus
            }.reduction()
        }

        /// Create a new `AlgebraicNumber`, automatically wrapping the modulus in an R using the
        /// assumption that `R: From<Polynomial<T>>`.
        ///
        /// Note that this function does not check if the modulus is irreducible, and so logic
        /// errors may occur if that is not the case.
        pub fn new(value: Polynomial<T>, denom: T, modulus: Polynomial<T>) -> Self {
            AlgebraicNumberR {
                value, denom,
                modulus: From::from(modulus),
            }.reduction()
        }

        pub(crate) fn reduce(&mut self) {
            let common_factor: T = gcd(&self.value.cont(), &self.denom);
            if common_factor.is_one() {
                return;
            }
            self.value.map_coeffs(|val| pdiv(val, common_factor.clone()));
            self.denom = pdiv(self.denom.clone(), common_factor);

            if self.denom < Zero::zero() {
                self.denom = -&self.denom;
                self.value = -&self.value;
            }
        }

        pub(crate) fn reduction(mut self) -> Self {
            self.reduce();
            self
        }

        pub fn mul_inverse(&self) -> Self {
            debug_assert!(!self.modulus.is_zero());
            // Compute B inverse mod T according to the method described on page 160

            let egcd = extended_gcd(self.value.clone(), Polynomial::clone(&*self.modulus));
            debug_assert!(egcd.gcd.coeffs().len() == 1);
            // TODO: Can some of these clones be avoided?
            AlgebraicNumberR {
                value: egcd.a * self.denom.clone(),
                denom: egcd.gcd.coeffs()[0].clone(),
                modulus: self.modulus.clone(),
            }.reduction()
        }
    }
    
    impl<'a, 'b, T, R, U> Add<&'b AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            let modulus = if self.modulus.is_zero() { &other.modulus } else { &self.modulus };
            // New denominator will be the lcm of the old denominators
            let new_denom = pdiv(&self.denom * &other.denom, gcd(&self.denom, &other.denom));
            let my_new_value = &self.value * &pdiv(new_denom.clone(), self.denom.clone());
            let other_new_value = &other.value * &pdiv(new_denom.clone(), other.denom.clone());
            let new_value = my_new_value + other_new_value;
            AlgebraicNumberR {
                value: new_value,
                denom: new_denom,
                modulus: modulus.clone(),
            }.reduction()
        }
    }
    
    // Boilerplate reference reimpls

    impl<'a, T, R, U> Add<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self + other
        }
    }

    impl<'a, T, R, U> Add<AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self + &other
        }
    }
    
    impl<T, R, U> Add<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self + &other
        }
    }

    impl<T, R, U> AddAssign<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn add_assign(&mut self, other: Self) {
            *self = &*self + other;
        }
    }

    impl<'a, T, R, U> AddAssign<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn add_assign(&mut self, other: &Self) {
            *self = &*self + other;
        }
    }
    
    // Boilerplate repeat, but for Sub instead of Add

    impl<'a, 'b, T, R, U> Sub<&'b AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            let modulus = if self.modulus.is_zero() { &other.modulus } else { &self.modulus };
            // New denominator will be the lcm of the old denominators
            let new_denom = pdiv(&self.denom * &other.denom, gcd(&self.denom, &other.denom));
            let my_new_value = &self.value * &pdiv(new_denom.clone(), self.denom.clone());
            let other_new_value = &other.value * &pdiv(new_denom.clone(), other.denom.clone());
            let new_value = my_new_value - other_new_value;
            AlgebraicNumberR {
                value: new_value,
                denom: new_denom,
                modulus: modulus.clone(),
            }.reduction()
        }
    }
    
    // Boilerplate reference reimpls

    impl<'a, T, R, U> Sub<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self - other
        }
    }

    impl<'a, T, R, U> Sub<AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self - &other
        }
    }
    
    impl<T, R, U> Sub<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self - &other
        }
    }

    impl<T, R, U> SubAssign<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn sub_assign(&mut self, other: Self) {
            *self = &*self - other;
        }
    }

    impl<'a, T, R, U> SubAssign<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn sub_assign(&mut self, other: &Self) {
            *self = &*self - other;
        }
    }
    
    // Multiplication impls

    impl<'a, 'b, T, R, U> Mul<&'b AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            let modulus = if self.modulus.is_zero() { &other.modulus } else { &self.modulus };
            let prod = &self.value * &other.value;
            AlgebraicNumberR {
                value: PseudoDivRem::pseudo_divrem(&prod, &*modulus).rem,
                denom: &self.denom * &other.denom,
                modulus: modulus.clone(),
            }.reduction()
        }
    }
    
    // Boilerplate reference reimpls

    impl<'a, T, R, U> Mul<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self * other
        }
    }

    impl<'a, T, R, U> Mul<AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * &other
        }
    }
    
    impl<T, R, U> Mul<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self * &other
        }
    }

    impl<T, R, U> MulAssign<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn mul_assign(&mut self, other: Self) {
            *self = &*self * other;
        }
    }

    impl<'a, T, R, U> MulAssign<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn mul_assign(&mut self, other: &Self) {
            *self = &*self * other;
        }
    }

    // Multiplicative Inverse implementation

    impl<'a, 'b, T, R, U> Div<&'a AlgebraicNumberR<T, R>> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
    
    impl<'b, T, R, U> Div<AlgebraicNumberR<T, R>> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
    
    impl<'a, T, R, U> Div<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
    
    impl<T, R, U> Div<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }

    impl<T, R, U> DivAssign<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn div_assign(&mut self, other: Self) {
            *self = &*self + other;
        }
    }

    impl<'a, T, R, U> DivAssign<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn div_assign(&mut self, other: &Self) {
            *self = &*self + other;
        }
    }

    macro_rules! impl_with_bounds {
        ($trait:ident) => {
            impl<T, R, U> $trait for AlgebraicNumberR<T, R>
            where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
                  R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone { }
        };
    }

    impl_with_bounds!(AssociativeAddition);
    impl_with_bounds!(CommutativeAddition);
    impl_with_bounds!(AssociativeMultiplication);
    impl_with_bounds!(CommutativeMultiplication);
    impl_with_bounds!(NonZeroProdProperty);
    impl_with_bounds!(ExactDivision);

    // Neg traits
    impl<'a, T, R, U> Neg for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn neg(self) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value = -this.value;
            this
        }
    }
    
    impl<T, R, U> Neg for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn neg(mut self) -> AlgebraicNumberR<T, R> {
            self.value = -self.value;
            self
        }
    }

    // Zero and One traits. These cheat by using the zero polynomial for modulus.

    impl<T, R, U> Zero for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn zero() -> Self {
            AlgebraicNumberR {
                value: Zero::zero(),
                denom: One::one(),
                modulus: R::from(Zero::zero()),
            }
        }

        fn is_zero(&self) -> bool {
            self.value.is_zero()
        }
    }
    
    impl<T, R, U> One for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn one() -> Self {
            AlgebraicNumberR {
                value: One::one(),
                denom: One::one(),
                modulus: R::from(Zero::zero()),
            }
        }

        fn is_one(&self) -> bool {
            self.value.is_one()
        }
    }

    // Scalar value operations

    impl<T, R, U> Add<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(mut self, other: T) -> AlgebraicNumberR<T, R> {
            self.value += other;
            self
        }
    }

    impl<'a, T, R, U> Add<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(mut self, other: &T) -> AlgebraicNumberR<T, R> {
            self.value += other;
            self
        }
    }

    impl<'b, T, R, U> Add<T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value += other;
            this
        }
    }

    impl<'a, 'b, T, R, U> Add<&'a T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: &T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value += other;
            this
        }
    }
    
    impl<T, R, U> Sub<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(mut self, other: T) -> AlgebraicNumberR<T, R> {
            self.value -= other;
            self
        }
    }

    impl<'a, T, R, U> Sub<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(mut self, other: &T) -> AlgebraicNumberR<T, R> {
            self.value -= other;
            self
        }
    }

    impl<'b, T, R, U> Sub<T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value -= other;
            this
        }
    }

    impl<'a, 'b, T, R, U> Sub<&'a T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: &T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value -= other;
            this
        }
    }
    
    impl<T, R, U> Mul<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(mut self, other: T) -> AlgebraicNumberR<T, R> {
            self.value *= other;
            self
        }
    }

    impl<'a, T, R, U> Mul<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(mut self, other: &T) -> AlgebraicNumberR<T, R> {
            self.value *= other;
            self
        }
    }

    impl<'b, T, R, U> Mul<T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value *= other;
            this
        }
    }

    impl<'a, 'b, T, R, U> Mul<&'a T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: &T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.value *= other;
            this
        }
    }
    
    impl<T, R, U> Div<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(mut self, other: T) -> AlgebraicNumberR<T, R> {
            self.denom *= other;
            self.reduction()
        }
    }

    impl<'a, T, R, U> Div<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(mut self, other: &T) -> AlgebraicNumberR<T, R> {
            self.denom *= other;
            self.reduction()
        }
    }

    impl<'b, T, R, U> Div<T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.denom *= other;
            this.reduction()
        }
    }

    impl<'a, 'b, T, R, U> Div<&'a T> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: &T) -> AlgebraicNumberR<T, R> {
            let mut this = self.clone();
            this.denom *= other;
            this.reduction()
        }
    }

    // *Assign traits for scalars
    
    impl<T, R, U> AddAssign<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn add_assign(&mut self, other: T) {
            self.value += other;
        }
    }

    impl<'a, T, R, U> AddAssign<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn add_assign(&mut self, other: &T) {
            self.value += other;
        }
    }
    
    impl<T, R, U> SubAssign<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn sub_assign(&mut self, other: T) {
            self.value -= other;
        }
    }

    impl<'a, T, R, U> SubAssign<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn sub_assign(&mut self, other: &T) {
            self.value -= other;
        }
    }
    
    impl<T, R, U> MulAssign<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn mul_assign(&mut self, other: T) {
            self.value *= other;
        }
    }

    impl<'a, T, R, U> MulAssign<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn mul_assign(&mut self, other: &T) {
            self.value *= other;
        }
    }
    
    impl<T, R, U> DivAssign<T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn div_assign(&mut self, other: T) {
            self.denom *= other;
            self.reduce()
        }
    }

    impl<'a, T, R, U> DivAssign<&'a T> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone + PartialOrd,
          R: Deref<Target=Polynomial<T>> + From<Polynomial<T>> + Clone {
        fn div_assign(&mut self, other: &T) {
            self.denom *= other;
            self.reduce()
        }
    }
}

/// An element of a number field isomorphic to Q\[x\]/<*modulus>. See [raw::AlgebraicNumberR] for
/// more functions and implementations.
pub type AlgebraicNumber<T> = raw::AlgebraicNumberR<T, Arc<Polynomial<T>>>;

impl<T, U> AlgebraicNumber<T>
where T: Ring + Gcd + PartialEq + PseudoDivRem<Output=T, MultType=U> + Clone + PartialOrd {
    // N/A yet! Nothing Arc specific.
}

assert_impl_all!(AlgebraicNumber<i32>: Field);
assert_impl_all!(AlgebraicNumber<Int>: Field);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        // Some select examples over Q(i)
        let modulus = Polynomial::from_coefficients(vec![1i64, 0, 1]);
        
        // Check that ((1 - i) / 2)^2 == -i / 2
        // 1 - i
        let val1 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![1, -1]),
            2,
            modulus.clone(),
        );
        // -i / 2
        let val2 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![0, -1]),
            2,
            modulus.clone(),
        );
        
        assert_eq!(&val1 * &val1, val2);

        // Check that (3 + 5i) / (-2 + 3i) == (9 - 19i) / 13
        // 3 + 5i
        let val1 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![3, 5]),
            1,
            modulus.clone(),
        );
        // -2 + 3i
        let val2 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![-2, 3]),
            1,
            modulus.clone(),
        );
        // (9 - 19i) / 13
        let val3 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![9, -19]),
            13,
            modulus.clone(),
        );

        assert_eq!(&val1 / &val2, val3);

        let modulus = Polynomial::from_coefficients(vec![Int::from(2), Int::from(-10), Int::from(0), Int::from(0), Int::from(0), Int::from(1)]);

        // Check that adding a root of x^5 - 10x + 2 works properly
        let val = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![Int::from(0), Int::from(1)]),
            Int::from(1),
            modulus.clone(),
        );
        let val10 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![Int::from(10)]),
            Int::from(1),
            modulus.clone(),
        );
        let val2 = AlgebraicNumber::new(
            Polynomial::from_coefficients(vec![Int::from(2)]),
            Int::from(1),
            modulus.clone(),
        );


        assert!(((&val).pow(&5) - &val * &val10 + &val2).is_zero());
        assert!(((&val).pow(&5) - &val * Int::from(10) + Int::from(2)).is_zero());

    }

    #[test]
    fn readme_test() {
        use std::sync::Arc;

        // Create the polynomial x^2 + 1
        let modulus = poly!(1;x^2 + 1;x^0);
        // Wrap it in an Arc
        let modulus = Arc::new(modulus);

        // Create the element (5 + 3i) / 2 in Q(i)
        let val_a = AlgebraicNumber::new_with_modulus(
            Polynomial::from_coefficients(vec![5, 3]), 2,
            Arc::clone(&modulus));

        // Create another element (5 - 3i) / 2
        let val_b = AlgebraicNumber::new_with_modulus(
            Polynomial::from_coefficients(vec![5, -3]), 2,
            Arc::clone(&modulus));

        assert_eq!((val_a * val_b).to_string(), "17 / 2");
    }
}
