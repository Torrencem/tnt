
use num_traits::{Zero, One};
use std::fmt::Write;
use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, MulAssign, SubAssign, DivAssign};
use std::cmp::max;

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct Polynomial<T> {
    // Coefficients of the polynomial, in order from lowest power of x to highest. Invariant:
    // there should never be trailing zeroes in this vector (so that the degree is always the
    // length of the vector
    cs: Vec<T>,
}

impl<T> Polynomial<T> {
    pub fn new() -> Self {
        Polynomial { cs: vec![] }
    }

    pub fn coefficients(&self) -> &[T] {
        &self.cs
    }
}

impl<T: Zero> Polynomial<T> {
    pub fn from_coefficients(mut coefficients: Vec<T>) -> Self {
        // Remove trailing zeroes from the coefficients
        let mut i = coefficients.len() - 1;
        loop {
            if coefficients[i].is_zero() {
                coefficients.pop();
            } else {
                break;
            }
            if i == 0 {
                break;
            } else {
                i -= 1;
            }
        }
        Polynomial {
            cs: coefficients
        }
    }
}

use std::fmt::Display;

impl<T: Display> Polynomial<T> {
    pub fn format<X: Display>(&self, var: X) -> String {
        let mut res = String::new();

        for (index, coeff) in self.cs.iter().rev().enumerate() {
            let exponent = self.cs.len() - 1 - index;
            
            if exponent == 0 {
                write!(res, " + {}", coeff).unwrap();
            } else if index == 0 {
                write!(res, "{}{}^{}", coeff, var, exponent).unwrap();
            } else {
                write!(res, " + {}{}^{}", coeff, var, exponent).unwrap();
            }
        }

        return res;
    }

    pub fn pretty_format<X: Display>(&self, var: X) -> String
    where T: Zero {
        let mut res = String::new();

        for (index, coeff) in self.cs.iter().rev().enumerate() {
            let exponent = self.cs.len() - 1 - index;

            if coeff.is_zero() {
                continue;
            }
            
            if exponent == 0 {
                write!(res, " + {}", coeff).unwrap();
            } else if index == 0 {
                write!(res, "{}{}^{}", coeff, var, exponent).unwrap();
            } else {
                write!(res, " + {}{}^{}", coeff, var, exponent).unwrap();
            }
        }

        return res;
    }
}

impl<T> Polynomial<T> {
    pub fn eval(&self, val: &T) -> T
    where
        T: for<'a> AddAssign<&'a T>,
        T: for<'a> MulAssign<&'a T>,
        T: Zero,
        T: Clone,
    {
        let mut res = Zero::zero();

        for coeff in self.cs.iter().rev() {
            res *= val;
            res += coeff;
        }

        res
    }
}

impl<T: Add<Output=T>> Add<Polynomial<T>> for Polynomial<T> {
    type Output = Polynomial<T>;
    
    fn add(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter().rev();
        let mut ocs = other.cs.into_iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval + oval);
                } else {
                    res.push(sval);
                    while let Some(sval) = scs.next() {
                        res.push(sval);
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(oval);
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

impl<'a, T: Clone> Add<Polynomial<T>> for &'a Polynomial<T>
where for<'b> &'b T: Add<T, Output=T> {
    type Output = Polynomial<T>;

    fn add(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter().rev();
        let mut ocs = other.cs.into_iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval + oval);
                } else {
                    res.push(sval.clone());
                    while let Some(sval) = scs.next() {
                        res.push(sval.clone());
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(oval);
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

impl<'a, T: Clone> Add<&'a Polynomial<T>> for Polynomial<T>
where for<'b> T: Add<&'b T, Output=T> {
    type Output = Polynomial<T>;

    fn add(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter().rev();
        let mut ocs = other.cs.iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval + oval);
                } else {
                    res.push(sval);
                    while let Some(sval) = scs.next() {
                        res.push(sval);
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(oval.clone());
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

impl<'a, 'b, T: Clone> Add<&'b Polynomial<T>> for &'a Polynomial<T>
where for<'c, 'd> &'c T: Add<&'d T, Output=T> {
    type Output = Polynomial<T>;

    fn add(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter().rev();
        let mut ocs = other.cs.iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval + oval);
                } else {
                    res.push(sval.clone());
                    while let Some(sval) = scs.next() {
                        res.push(sval.clone());
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(oval.clone());
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

// Boilerplate: repeat the same exact thing for subtraction

impl<T: Sub<Output=T> + Neg<Output=T>> Sub<Polynomial<T>> for Polynomial<T> {
    type Output = Polynomial<T>;
    
    fn sub(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter().rev();
        let mut ocs = other.cs.into_iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval - oval);
                } else {
                    res.push(sval);
                    while let Some(sval) = scs.next() {
                        res.push(sval);
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(-oval);
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

impl<'a, T: Clone + Neg<Output=T>> Sub<Polynomial<T>> for &'a Polynomial<T>
where for<'b> &'b T: Sub<T, Output=T> {
    type Output = Polynomial<T>;

    fn sub(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter().rev();
        let mut ocs = other.cs.into_iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval - oval);
                } else {
                    res.push(sval.clone());
                    while let Some(sval) = scs.next() {
                        res.push(sval.clone());
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(-oval);
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

impl<'a, T: Clone + Neg<Output=T>> Sub<&'a Polynomial<T>> for Polynomial<T>
where for<'b> T: Sub<&'b T, Output=T> {
    type Output = Polynomial<T>;

    fn sub(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter().rev();
        let mut ocs = other.cs.iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval - oval);
                } else {
                    res.push(sval);
                    while let Some(sval) = scs.next() {
                        res.push(sval);
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(-oval.clone());
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

impl<'a, 'b, T: Clone + Neg<Output=T>> Sub<&'b Polynomial<T>> for &'a Polynomial<T>
where for<'c, 'd> &'c T: Sub<&'d T, Output=T> {
    type Output = Polynomial<T>;

    fn sub(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter().rev();
        let mut ocs = other.cs.iter().rev();
        let mut res: Vec<T> = Vec::with_capacity(max(scs.len(), ocs.len()));
        loop {
            if let Some(sval) = scs.next() {
                if let Some(oval) = ocs.next() {
                    res.push(sval - oval);
                } else {
                    res.push(sval.clone());
                    while let Some(sval) = scs.next() {
                        res.push(sval.clone());
                    }
                    break;
                }
            } else {
                while let Some(oval) = ocs.next() {
                    res.push(-oval.clone());
                }
                break;
            }
        }
        Polynomial {
            cs: res
        }
    }
}

// Polynomial Multiplication algorithm

impl<'c, 'd, T: Zero> Mul<&'c Polynomial<T>> for &'d Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: &Polynomial<T>) -> Self::Output {
        if self.cs.len() == 0 || other.cs.len() == 0{
            return Polynomial { cs: vec![] };
        }
        let mut res = (0..self.cs.len() + other.cs.len() - 1)
            .map(|_| Zero::zero())
            .collect::<Vec<_>>();

        let mut ia = 0;
        let mut ib = 0;
        let l = max(self.cs.len(), other.cs.len());
        while ia <= l {
            loop {
                if let (Some(x), Some(y)) = (self.cs.get(ia), other.cs.get(ib)) {
                    res[ia + ib] += x * y;
                }
                if ia == 0 {
                    break;
                } else {
                    ia -= 1;
                    ib += 1;
                }
            }
            ia = ib + 1;
            ib = 0;
        }

        Polynomial::from_coefficients(res)
    }
}

impl<T: Zero> Mul<Polynomial<T>> for Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: Polynomial<T>) -> Self::Output {
        (&self) * (&other)
    }
}

impl<'c, T: Zero> Mul<Polynomial<T>> for &'c Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: Polynomial<T>) -> Self::Output {
        self * (&other)
    }
}

impl<'d, T: Zero> Mul<&'d Polynomial<T>> for Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: &Polynomial<T>) -> Self::Output {
        (&self) * other
    }
}

// The *Assign traits
// TODO: You should be able to add and mul T to Polynomial<T>

impl<T: Clone> AddAssign for Polynomial<T>
where for<'c, 'd> &'c T: Add<&'d T, Output=T> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}

impl<'a, T: Clone> AddAssign<&'a Polynomial<T>> for Polynomial<T>
where for<'c, 'd> &'c T: Add<&'d T, Output=T> {
    fn add_assign(&mut self, other: &Self) {
        *self = &*self + other;
    }
}

impl<T: Clone + Neg<Output=T>> SubAssign for Polynomial<T>
where for<'c, 'd> &'c T: Sub<&'d T, Output=T> {
    fn sub_assign(&mut self, other: Self) {
        *self = &*self - &other;
    }
}

impl<'a, T: Clone + Neg<Output=T>> SubAssign<&'a Polynomial<T>> for Polynomial<T>
where for<'c, 'd> &'c T: Sub<&'d T, Output=T> {
    fn sub_assign(&mut self, other: &Self) {
        *self = &*self - &other;
    }
}

impl<T: Zero + AddAssign> MulAssign for Polynomial<T>
where for<'c, 'd> &'c T: Mul<&'d T, Output=T> {
    fn mul_assign(&mut self, other: Self) {
        *self = &*self * &other;
    }
}

impl<'a, T: Zero + AddAssign> MulAssign<&'a Polynomial<T>> for Polynomial<T>
where for<'c, 'd> &'c T: Mul<&'d T, Output=T> {
    fn mul_assign(&mut self, other: &Self) {
        *self = &*self * other;
    }
}

// Traits for Polynomial<T> and T

impl<T> Add<T> for Polynomial<T>
where T: AddAssign {
    type Output = Polynomial<T>;

    fn add(mut self, other: T) -> Self {
        self.cs[0] += other;
        self
    }
}

impl<'a, T> Add<&'a T> for Polynomial<T>
where for<'b> T: AddAssign<&'b T> {
    type Output = Polynomial<T>;

    fn add(mut self, other: &T) -> Self {
        self.cs[0] += other;
        self
    }
}

impl<'a, T: Clone> Add<T> for &'a Polynomial<T>
where T: AddAssign<T> {
    type Output = Polynomial<T>;

    fn add(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] += other;
        p
    }
}

impl<'a, 'b, T: Clone> Add<&'a T> for &'b Polynomial<T>
where for<'c> T: AddAssign<&'c T> {
    type Output = Polynomial<T>;

    fn add(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] += other;
        p
    }
}

impl<T> Sub<T> for Polynomial<T>
where T: SubAssign {
    type Output = Polynomial<T>;

    fn sub(mut self, other: T) -> Self {
        self.cs[0] -= other;
        self
    }
}

impl<'a, T> Sub<&'a T> for Polynomial<T>
where for<'b> T: SubAssign<&'b T> {
    type Output = Polynomial<T>;

    fn sub(mut self, other: &T) -> Self {
        self.cs[0] -= other;
        self
    }
}

impl<'a, T: Clone> Sub<T> for &'a Polynomial<T>
where T: SubAssign<T> {
    type Output = Polynomial<T>;

    fn sub(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] -= other;
        p
    }
}

impl<'a, 'b, T: Clone> Sub<&'a T> for &'b Polynomial<T>
where for<'c> T: SubAssign<&'c T> {
    type Output = Polynomial<T>;

    fn sub(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] -= other;
        p
    }
}

impl<T> Mul<T> for Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(mut self, other: T) -> Self {
        for val in self.cs.iter_mut() {
            *val *= &other;
        }
        self
    }
}

impl<'a, T> Mul<&'a T> for Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(mut self, other: &T) -> Self {
        for val in self.cs.iter_mut() {
            *val *= other;
        }
        self
    }
}

impl<'a, T: Clone> Mul<T> for &'a Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val *= &other;
        }
        p
    }
}

impl<'a, 'b, T: Clone> Mul<&'a T> for &'b Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val *= other;
        }
        p
    }
}

impl<T> Div<T> for Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(mut self, other: T) -> Self {
        for val in self.cs.iter_mut() {
            *val /= &other;
        }
        self
    }
}

impl<'a, T> Div<&'a T> for Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(mut self, other: &T) -> Self {
        for val in self.cs.iter_mut() {
            *val /= other;
        }
        self
    }
}

impl<'a, T: Clone> Div<T> for &'a Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val /= &other;
        }
        p
    }
}

impl<'a, 'b, T: Clone> Div<&'a T> for &'b Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val /= other;
        }
        p
    }
}

impl<T: Add<Output=T>> Zero for Polynomial<T> {
    fn zero() -> Self {
        Polynomial { cs: vec![] }
    }
    fn is_zero(&self) -> bool {
        self.cs.is_empty()
    }
}

impl<T: One + Zero + AddAssign + PartialEq> One for Polynomial<T> 
where for<'c, 'd> &'c T: Mul<&'d T, Output=T> {
    fn one() -> Self {
        Polynomial { cs: vec![One::one()] }
    }
    fn is_one(&self) -> bool {
        self.cs.len() == 1 && self.cs[0].is_one()
    }
}

// Neg
impl<T: Neg<Output=T>> Neg for Polynomial<T> {
    type Output = Self;

    fn neg(mut self) -> Self {
        Polynomial { cs: self.cs.drain(..)
            .map(|val| -val)
            .collect() }
    }
}

impl<'a, T: Clone + Neg<Output=T>> Neg for &'a Polynomial<T> {
    type Output = Polynomial<T>;

    fn neg(self) -> Polynomial<T> {
        -(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiplication() {
        let p = Polynomial::from_coefficients(vec![1, 1, 1]);
        let q = Polynomial::from_coefficients(vec![1, 1]);
        assert_eq!((p * q).cs, vec![1, 2, 2, 1]);
        
        let p = Polynomial::from_coefficients(vec![1, 1, -20, 15, -5, 10]);
        assert_eq!(p.eval(&2), 283);
        
        let p = Polynomial::from_coefficients(vec![1, 1, -20]);
        assert_eq!(p.eval(&2), -77);
    }
}
