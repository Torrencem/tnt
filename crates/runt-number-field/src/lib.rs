
// Right now, this module doesn't use the r_i trick described on page 159-160 of the course book
// for multiplication, and so it could serve to be a little more time efficient.

use runt_polynomial::*;
use runt_alg::*;
use std::ops::Deref;
use std::sync::Arc;

/// A module containing the raw types that are generic over the container of their reference to
/// their fields.
pub mod raw {
    use super::*;
    use std::ops::*;

    /// An element of a number field isomorphic to Q[x]/<*modulus>
    /// This struct is generic over the way to references the modulus. You could, for example, use
    /// AlgebraicNumberR<T, Arc<Polynomial<T>>> to keep track of a shared reference to a
    /// polynomial, or use AlgebraicNumberR<T, &'a Polynomial<T>> if you're able to bound the
    /// lifetime of the modulus for efficiency.
    #[derive(Clone, Debug)]
    pub struct AlgebraicNumberR<T, R> 
    where T: Ring + Gcd + PartialEq,
          R: Deref<Target=Polynomial<T>> + Clone,
    {
        pub(crate) value: Polynomial<T>,
        pub(crate) denom: T,
        pub(crate) modulus: R,
    }

    impl<T, R> PartialEq for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq,
          R: Deref<Target=Polynomial<T>> + Clone,
    {
        fn eq(&self, other: &Self) -> bool {
            self.value == other.value && self.denom == other.denom
        }
    }

    impl<T, R> Eq for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq,
          R: Deref<Target=Polynomial<T>> + Clone,
    { }

    // TODO: multiplicative inverse, Div impls

    impl<T, R, U> AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PartialEq + PseudoDivRem<Output=T, MultType=U> + Clone,
          R: Deref<Target=Polynomial<T>> + Clone,
    {
        pub fn new_raw(value: Polynomial<T>, denom: T, modulus: R) -> Self {
            AlgebraicNumberR {
                value, denom, modulus
            }.reduction()
        }

        pub fn reduce(&mut self) {
            let common_factor: T = gcd(&self.value.cont(), &self.denom);
            if common_factor.is_one() {
                return;
            }
            self.value.map_coeffs(|val| pdiv(val, common_factor.clone()));
            self.denom = pdiv(self.denom.clone(), common_factor);
        }

        pub fn reduction(mut self) -> Self {
            self.reduce();
            self
        }

        pub fn mul_inverse(&self) -> Self {
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
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            // New denominator will be the lcm of the old denominators
            let new_denom = pdiv(&self.denom * &other.denom, gcd(&self.denom, &other.denom));
            let my_new_value = &self.value * &pdiv(new_denom.clone(), self.denom.clone());
            let other_new_value = &other.value * &pdiv(new_denom.clone(), other.denom.clone());
            let new_value = my_new_value + other_new_value;
            AlgebraicNumberR {
                value: new_value,
                denom: new_denom,
                modulus: self.modulus.clone(),
            }.reduction()
        }
    }
    
    // Boilerplate reference reimpls

    impl<'a, T, R, U> Add<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self + other
        }
    }

    impl<'a, T, R, U> Add<AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self + &other
        }
    }
    
    impl<T, R, U> Add<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn add(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self + &other
        }
    }
    
    // Boilerplate repeat, but for Sub instead of Add

    impl<'a, 'b, T, R, U> Sub<&'b AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            // New denominator will be the lcm of the old denominators
            let new_denom = pdiv(&self.denom * &other.denom, gcd(&self.denom, &other.denom));
            let my_new_value = &self.value * &pdiv(new_denom.clone(), self.denom.clone());
            let other_new_value = &other.value * &pdiv(new_denom.clone(), other.denom.clone());
            let new_value = my_new_value - other_new_value;
            AlgebraicNumberR {
                value: new_value,
                denom: new_denom,
                modulus: self.modulus.clone(),
            }.reduction()
        }
    }
    
    // Boilerplate reference reimpls

    impl<'a, T, R, U> Sub<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self - other
        }
    }

    impl<'a, T, R, U> Sub<AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self - &other
        }
    }
    
    impl<T, R, U> Sub<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn sub(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self - &other
        }
    }
    
    // Multiplication impls

    impl<'a, 'b, T, R, U> Mul<&'b AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            let prod = &self.value * &other.value;
            AlgebraicNumberR {
                value: PseudoDivRem::pseudo_divrem(&prod, &*self.modulus).rem,
                denom: &self.denom * &other.denom,
                modulus: self.modulus.clone(),
            }.reduction()
        }
    }
    
    // Boilerplate reference reimpls

    impl<'a, T, R, U> Mul<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self * other
        }
    }

    impl<'a, T, R, U> Mul<AlgebraicNumberR<T, R>> for &'a AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * &other
        }
    }
    
    impl<T, R, U> Mul<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn mul(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            &self * &other
        }
    }

    // Multiplicative Inverse implementation

    impl<'a, 'b, T, R, U> Div<&'a AlgebraicNumberR<T, R>> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
    
    impl<'b, T, R, U> Div<AlgebraicNumberR<T, R>> for &'b AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
    
    impl<'a, T, R, U> Div<&'a AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: &AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
    
    impl<T, R, U> Div<AlgebraicNumberR<T, R>> for AlgebraicNumberR<T, R>
    where T: Ring + Gcd + PseudoDivRem<Output=T, MultType=U> + PartialEq + Clone,
          R: Deref<Target=Polynomial<T>> + Clone {
        type Output = AlgebraicNumberR<T, R>;

        fn div(self, other: AlgebraicNumberR<T, R>) -> AlgebraicNumberR<T, R> {
            self * other.mul_inverse()
        }
    }
}

pub type AlgebraicNumber<T> = raw::AlgebraicNumberR<T, Arc<Polynomial<T>>>;

impl<T, U> AlgebraicNumber<T>
where T: Ring + Gcd + PartialEq + PseudoDivRem<Output=T, MultType=U> + Clone {
    // TODO: Make this constuctor in a way that multiple AlgebraicNumbers can share a reference to
    // the same modulus
    // NOTE that this function doesn't assert that `modulus` is irreducible, and so logic errors
    // may occur if that's not the case.
    pub fn new(value: Polynomial<T>, denom: T, modulus: Polynomial<T>) -> Self {
        raw::AlgebraicNumberR {
            value, denom,
            modulus: Arc::new(modulus),
        }.reduction()
    }
}

// TODO: Complete this module, add test cases

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
            // This clone will probably be unnecessary eventually
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
    }
}
