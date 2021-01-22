
#![feature(iterator_fold_self)]

use num_traits::{Zero, One};
use std::fmt::Write;
use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, MulAssign, SubAssign, DivAssign};
use std::cmp::max;
use std::borrow::Cow;
use runt_alg::*;

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

    pub fn coeffs(&self) -> &[T] {
        &self.cs
    }

    pub fn map_coeffs<F>(&mut self, f: F)
    where F: FnMut(T) -> T,
          T: Zero {
        let new_cs = self.cs.drain(..)
            .map(f)
            .collect();
        *self = Polynomial::from_coefficients(new_cs);
    }

    pub fn degree(&self) -> usize {
        self.cs.len().saturating_sub(1)
    }

    /// Get the leading coefficient of this polynomial. Note that this returns either a reference
    /// to the leading coefficient, or an owned zero value if this polynomial is zero, so you will
    /// need to call .into_owned() to get a T or .as_ref() to get a &T.
    pub fn lc<'a>(&'a self) -> Cow<'a, T> 
    where T: Clone + Zero
    {
        if self.cs.len() == 0 {
            Cow::Owned(Zero::zero())
        } else {
            Cow::Borrowed(&self.cs[self.cs.len() - 1])
        }
    }

    pub fn cont(&self) -> T
    where T: Gcd + Zero + Clone {
        self.cs.iter()
            .cloned()
            .fold_first(|a, b| {
                gcd(&a, &b)
            })
            .unwrap_or_else(|| Zero::zero())
    }

    pub fn pp<U>(&self) -> Polynomial<T>
    where T: Gcd + Zero + Clone + PseudoDivRem<Output=T, MultType=U>
    {
        let c = self.cont();
        let mut p = self.clone();
        for coeff in p.cs.iter_mut() {
            let val = PseudoDivRem::pseudo_divrem(coeff.clone(), c.clone()).div;
            *coeff = val;
        }
        p
    }
}

impl<T: Zero> Polynomial<T> {
    pub fn from_coefficients(mut coefficients: Vec<T>) -> Self {
        if coefficients.is_empty() {
            return Polynomial { cs: vec![] };
        }
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

    pub fn fix_coefficients(&mut self) {
        if self.cs.len() == 0 {
            return;
        }
        let mut i = self.cs.len() - 1;
        loop {
            if self.cs[i].is_zero() {
                self.cs.pop();
            } else {
                break;
            }
            if i == 0 {
                break;
            } else {
                i -= 1;
            }
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

impl<T: Add<Output=T> + Zero> Add<Polynomial<T>> for Polynomial<T> {
    type Output = Polynomial<T>;
    
    fn add(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter();
        let mut ocs = other.cs.into_iter();
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
        Polynomial::from_coefficients(res)
    }
}

impl<'a, T: Clone + Zero> Add<Polynomial<T>> for &'a Polynomial<T>
where for<'b> &'b T: Add<T, Output=T> {
    type Output = Polynomial<T>;

    fn add(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter();
        let mut ocs = other.cs.into_iter();
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
        Polynomial::from_coefficients(res)
    }
}

impl<'a, T: Clone + Zero> Add<&'a Polynomial<T>> for Polynomial<T>
where for<'b> T: Add<&'b T, Output=T> {
    type Output = Polynomial<T>;

    fn add(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter();
        let mut ocs = other.cs.iter();
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
        Polynomial::from_coefficients(res)
    }
}

impl<'a, 'b, T: Clone + Zero> Add<&'b Polynomial<T>> for &'a Polynomial<T>
where for<'c, 'd> &'c T: Add<&'d T, Output=T> {
    type Output = Polynomial<T>;

    fn add(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter();
        let mut ocs = other.cs.iter();
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
        Polynomial::from_coefficients(res)
    }
}

// Boilerplate: repeat the same exact thing for subtraction

impl<T: Sub<Output=T> + Neg<Output=T> + Zero> Sub<Polynomial<T>> for Polynomial<T> {
    type Output = Polynomial<T>;
    
    fn sub(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter();
        let mut ocs = other.cs.into_iter();
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
        Polynomial::from_coefficients(res)
    }
}

impl<'a, T: Clone + Neg<Output=T> + Zero> Sub<Polynomial<T>> for &'a Polynomial<T>
where for<'b> &'b T: Sub<T, Output=T> {
    type Output = Polynomial<T>;

    fn sub(self, other: Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter();
        let mut ocs = other.cs.into_iter();
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
        Polynomial::from_coefficients(res)
    }
}

impl<'a, T: Clone + Neg<Output=T> + Zero> Sub<&'a Polynomial<T>> for Polynomial<T>
where for<'b> T: Sub<&'b T, Output=T> {
    type Output = Polynomial<T>;

    fn sub(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.into_iter();
        let mut ocs = other.cs.iter();
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
        Polynomial::from_coefficients(res)
    }
}

impl<'a, 'b, T: Clone + Neg<Output=T> + Zero> Sub<&'b Polynomial<T>> for &'a Polynomial<T>
where for<'c, 'd> &'c T: Sub<&'d T, Output=T> {
    type Output = Polynomial<T>;

    fn sub(self, other: &Polynomial<T>) -> Self::Output {
        let mut scs = self.cs.iter();
        let mut ocs = other.cs.iter();
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
        Polynomial::from_coefficients(res)
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

        for i in 0..self.cs.len() {
            for j in 0..other.cs.len() {
                res[i + j] += &self.coeffs()[i] * &other.coeffs()[j];
            }
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

impl<T: Clone + Zero> AddAssign for Polynomial<T>
where for<'c, 'd> &'c T: Add<&'d T, Output=T> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
        self.fix_coefficients();
    }
}

impl<'a, T: Clone + Zero> AddAssign<&'a Polynomial<T>> for Polynomial<T>
where for<'c, 'd> &'c T: Add<&'d T, Output=T> {
    fn add_assign(&mut self, other: &Self) {
        *self = &*self + other;
        self.fix_coefficients();
    }
}

impl<T: Clone + Neg<Output=T> + Zero> SubAssign for Polynomial<T>
where for<'c, 'd> &'c T: Sub<&'d T, Output=T> {
    fn sub_assign(&mut self, other: Self) {
        *self = &*self - &other;
        self.fix_coefficients();
    }
}

impl<'a, T: Clone + Neg<Output=T> + Zero> SubAssign<&'a Polynomial<T>> for Polynomial<T>
where for<'c, 'd> &'c T: Sub<&'d T, Output=T> {
    fn sub_assign(&mut self, other: &Self) {
        *self = &*self - &other;
        self.fix_coefficients();
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
where T: AddAssign + Zero {
    type Output = Polynomial<T>;

    fn add(mut self, other: T) -> Self {
        self.cs[0] += other;
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Zero> Add<&'a T> for Polynomial<T>
where for<'b> T: AddAssign<&'b T> {
    type Output = Polynomial<T>;

    fn add(mut self, other: &T) -> Self {
        self.cs[0] += other;
        self.fix_coefficients();
        self
    }
}

impl<T> AddAssign<T> for Polynomial<T>
where T: AddAssign + Zero {
    fn add_assign(&mut self, other: T) {
        self.cs[0] += other;
        self.fix_coefficients();
    }
}

impl<'a, T: Zero> AddAssign<&'a T> for Polynomial<T>
where for<'b> T: AddAssign<&'b T> {
    fn add_assign(&mut self, other: &T) {
        self.cs[0] += other;
        self.fix_coefficients();
    }
}

impl<'a, T: Clone> Add<T> for &'a Polynomial<T>
where T: AddAssign<T> + Zero {
    type Output = Polynomial<T>;

    fn add(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] += other;
        p.fix_coefficients();
        p
    }
}

impl<'a, 'b, T: Clone + Zero> Add<&'a T> for &'b Polynomial<T>
where for<'c> T: AddAssign<&'c T> {
    type Output = Polynomial<T>;

    fn add(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] += other;
        p.fix_coefficients();
        p
    }
}

impl<T> Sub<T> for Polynomial<T>
where T: SubAssign + Zero {
    type Output = Polynomial<T>;

    fn sub(mut self, other: T) -> Self {
        self.cs[0] -= other;
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Zero> Sub<&'a T> for Polynomial<T>
where for<'b> T: SubAssign<&'b T> {
    type Output = Polynomial<T>;

    fn sub(mut self, other: &T) -> Self {
        self.cs[0] -= other;
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Clone> Sub<T> for &'a Polynomial<T>
where T: SubAssign<T> + Zero {
    type Output = Polynomial<T>;

    fn sub(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] -= other;
        p.fix_coefficients();
        p
    }
}

impl<'a, 'b, T: Clone + Zero> Sub<&'a T> for &'b Polynomial<T>
where for<'c> T: SubAssign<&'c T> {
    type Output = Polynomial<T>;

    fn sub(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        p.cs[0] -= other;
        p.fix_coefficients();
        p
    }
}

impl<T> SubAssign<T> for Polynomial<T>
where T: SubAssign + Zero {
    fn sub_assign(&mut self, other: T) {
        self.cs[0] -= other;
        self.fix_coefficients();
    }
}

impl<'a, T: Zero> SubAssign<&'a T> for Polynomial<T>
where for<'b> T: SubAssign<&'b T> + Zero {
    fn sub_assign(&mut self, other: &T) {
        self.cs[0] -= other;
        self.fix_coefficients();
    }
}

impl<T: Zero> Mul<T> for Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(mut self, other: T) -> Self {
        for val in self.cs.iter_mut() {
            *val *= &other;
        }
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Zero> Mul<&'a T> for Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(mut self, other: &T) -> Self {
        for val in self.cs.iter_mut() {
            *val *= other;
        }
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Clone + Zero> Mul<T> for &'a Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val *= &other;
        }
        p.fix_coefficients();
        p
    }
}

impl<'a, 'b, T: Clone + Zero> Mul<&'a T> for &'b Polynomial<T>
where for<'c> T: MulAssign<&'c T> {
    type Output = Polynomial<T>;

    fn mul(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val *= other;
        }
        p.fix_coefficients();
        p
    }
}

impl<T: Zero> MulAssign<T> for Polynomial<T>
where for<'b> T: MulAssign<&'b T> {
    fn mul_assign(&mut self, other: T) {
        for val in self.cs.iter_mut() {
            *val *= &other;
        }
        self.fix_coefficients();
    }
}

impl<'a, T: Zero> MulAssign<&'a T> for Polynomial<T>
where for<'b> T: MulAssign<&'b T> {
    fn mul_assign(&mut self, other: &T) {
        for val in self.cs.iter_mut() {
            *val *= other;
        }
        self.fix_coefficients();
    }
}

impl<T: Zero> Div<T> for Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(mut self, other: T) -> Self {
        for val in self.cs.iter_mut() {
            *val /= &other;
        }
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Zero> Div<&'a T> for Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(mut self, other: &T) -> Self {
        for val in self.cs.iter_mut() {
            *val /= other;
        }
        self.fix_coefficients();
        self
    }
}

impl<'a, T: Clone + Zero> Div<T> for &'a Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(self, other: T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val /= &other;
        }
        p.fix_coefficients();
        p
    }
}

impl<'a, 'b, T: Clone + Zero> Div<&'a T> for &'b Polynomial<T>
where for<'c> T: DivAssign<&'c T> {
    type Output = Polynomial<T>;

    fn div(self, other: &T) -> Polynomial<T> {
        let mut p = self.clone();
        for val in p.cs.iter_mut() {
            *val /= other;
        }
        p.fix_coefficients();
        p
    }
}

impl<T: Zero> DivAssign<T> for Polynomial<T>
where for<'b> T: DivAssign<&'b T> {
    fn div_assign(&mut self, other: T) {
        for val in self.cs.iter_mut() {
            *val /= &other;
        }
        self.fix_coefficients();
    }
}

impl<'a, T: Zero> DivAssign<&'a T> for Polynomial<T>
where for<'b> T: DivAssign<&'b T> {
    fn div_assign(&mut self, other: &T) {
        for val in self.cs.iter_mut() {
            *val /= other;
        }
        self.fix_coefficients();
    }
}

impl<T: Add<Output=T> + Zero> Zero for Polynomial<T> {
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

// Implementations for runt-algs traits
// TODO: Take a closer look at these weird requirements
impl<T: AssociativeAddition> AssociativeAddition for Polynomial<T> {}
impl<T: CommutativeAddition + AssociativeAddition> CommutativeAddition for Polynomial<T> {}
impl<T: AssociativeMultiplication + CommutativeAddition + AssociativeAddition> AssociativeMultiplication for Polynomial<T> {}
impl<T: CommutativeMultiplication + CommutativeAddition + AssociativeAddition> CommutativeMultiplication for Polynomial<T> {}
impl<T: NonZeroProdProperty> NonZeroProdProperty for Polynomial<T> {}

impl<'a, T: Ring + Power<Power=u64, Output=T> + Clone> PseudoDivRem for &'a Polynomial<T> {
    type Output = Polynomial<T>;
    type MultType = T;

    fn pseudo_divrem(u: &Polynomial<T>, v: &Polynomial<T>) -> PseudoDivRemResult<Polynomial<T>, T> {
        // Algorithm R from Knuth Vol 2
        // This has a tendancy to blow up mul and overflow, so there's an alternative algorithm
        // that does the divison without returning mul
        let m = u.degree();
        let n = v.degree();
        if m < n {
            return PseudoDivRemResult {
                mul: One::one(),
                div: Zero::zero(),
                rem: u.clone(),
            };
        }
        let mut u = u.coeffs().to_owned();
        let mut qs = vec![Zero::zero(); m - n + 1];
        for k in (0..=(m - n)).rev() {
            qs[k] = &u[n + k] * &v.lc().into_owned().pow(&(k as u64));
            for j in (0..=(n + k - 1)).rev() {
                if j < k {
                    u[j] = v.lc().as_ref() * &u[j];
                } else {
                    u[j] = &(v.lc().as_ref() * &u[j]) - &(&u[n + k] * &v.coeffs()[j - k]);
                }
            }
        }
        u.truncate(n);
        PseudoDivRemResult {
            mul: v.lc().into_owned().pow(&((m - n + 1) as u64)),
            div: Polynomial::from_coefficients(qs),
            rem: Polynomial::from_coefficients(u),
        }
    }
}

/// Finds the quotient of two polynomials, without returning the remainder. This is the same as
/// pseudo_divrem, except it doesn't compute mul, which prevents overflows.
fn poly_div<T: Ring + Power<Power=u64, Output=T> + Clone>(u: &Polynomial<T>, v: &Polynomial<T>) -> Polynomial<T> {
        // Algorithm R from Knuth Vol 2
        let m = u.degree();
        let n = v.degree();
        if m < n {
            return Zero::zero();
        }
        let mut u = u.coeffs().to_owned();
        let mut qs = vec![Zero::zero(); m - n + 1];
        for k in (0..=(m - n)).rev() {
            qs[k] = &u[n + k] * &v.lc().into_owned().pow(&(k as u64));
            for j in (0..=(n + k - 1)).rev() {
                if j < k {
                    u[j] = v.lc().as_ref() * &u[j];
                } else {
                    u[j] = &(v.lc().as_ref() * &u[j]) - &(&u[n + k] * &v.coeffs()[j - k]);
                }
            }
        }
        Polynomial::from_coefficients(qs)
}

impl<T> EuclideanFunction for Polynomial<T> {
    type Order = usize;

    fn norm(&self) -> usize {
        self.degree()
    }
}

impl<T, U> Gcd for Polynomial<T>
where T: PseudoDivRem<Output=T, MultType=U> + Gcd + Zero + Clone + Ring + Power,
      // for<'c> T: MulAssign<&'c T>,
      for<'c> &'c Polynomial<T>: PseudoDivRem<Output=Polynomial<T>, MultType=T>,
      // T: Display // Temporary
{
    fn gcd(&self, other: &Self) -> Polynomial<T> {
        let u = self.clone();
        let v = other.clone();
        // Algorithm E from Knuth vol 2
        // TODO: Replace this with algorithm C for efficiency
        let d = gcd(&u.cont(), &v.cont());
        let mut u = u.pp();
        let mut v = v.pp();
        
        loop {
            let pdr = PseudoDivRem::pseudo_divrem(&u, &v);
            let r = pdr.rem;
            if r.is_zero() {
                break v * d;
            } else if r.degree() == 0 {
                break Polynomial::from_coefficients(vec![d]);
            }
            u = v;
            v = r.pp();
        }
    }
}

// It should be the case that ua + vb = k*gcd(u, v) for some constant k
#[derive(Debug, Clone)]
pub struct ExtendedGcdResult<T> {
    pub a: Polynomial<T>,
    pub b: Polynomial<T>,
    pub gcd: Polynomial<T>,
}

pub fn extended_gcd<U, T: PseudoDivRem<Output=T, MultType=U>>(mut a: Polynomial<T>, mut b: Polynomial<T>) -> ExtendedGcdResult<T>
where T: Gcd + Zero + Clone + Ring + Power<Power=u64, Output=T> + One,
      for<'c> &'c Polynomial<T>: PseudoDivRem<Output=Polynomial<T>, MultType=T>,
      Polynomial<T>: Gcd,
      T: AddAssign,
{
    let mut swapped = false;
    if a.degree() < b.degree() {
        swapped = true;
        std::mem::swap(&mut a, &mut b);
    }

    let mut r_prev = a;
    let mut r = b;
    let mut s_prev: Polynomial<T> = Polynomial::from_coefficients(vec![One::one()]);
    let mut s: Polynomial<T> = Polynomial::from_coefficients(vec![]);
    let mut t_prev: Polynomial<T> = Polynomial::from_coefficients(vec![]);
    let mut t: Polynomial<T> = Polynomial::from_coefficients(vec![One::one()]);

    while r.degree() != 0 {
        let q: Polynomial<T> = poly_div(&r_prev, &r);
        let d: T = r.lc().into_owned();
        let e = r_prev.degree() - r.degree() + 1;
        let new_r = r_prev * d.clone().pow(&(e as u64)) - &q * &r;
        let new_s = s_prev * d.clone().pow(&(e as u64)) - &q * &s;
        let new_t = t_prev * d.clone().pow(&(e as u64)) - &q * &t;
        r_prev = r.clone();
        r = new_r;
        s_prev = s.clone();
        s = new_s;
        t_prev = t.clone();
        t = new_t;
    }

    if swapped {
        ExtendedGcdResult {
            a: t,
            b: s,
            gcd: r,
        }
    } else {
        ExtendedGcdResult {
            a: s,
            b: t,
            gcd: r,
        }
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

    #[test]
    fn test_misc() {
        let p = Polynomial::from_coefficients(vec![2, 1, -3]);
        let q = Polynomial::from_coefficients(vec![2, 0, 0, 1, -3]);
        println!("{}", (&p + &q).pretty_format("x"))
    }

    #[test]
    fn test_polynomial_division() {
        for u in &[
            vec![-5, 2, 8, -3, -3, 0, 1, 0, 1],
            vec![2, 3, 4, 10, -1, 2, 0, -2, -3],
            vec![1, 3, 2, 3, 2, 3, 4, 2, 2, 0, -1],
        ] {
            for v in &[
                vec![21, -9, -4, 0, 5, 0, 3],
                vec![0, 0, 0, 1],
                vec![2, 3, 4, 5, 4, 3],
            ] {
                let u = Polynomial::from_coefficients(u.clone());
                let v = Polynomial::from_coefficients(v.clone());
                // println!("u = {}, v = {}", u.pretty_format("x"), v.pretty_format("x"));
                let dr = PseudoDivRem::pseudo_divrem(&u, &v);
                // All results should always satisfy this relation
                let lhs = &u * dr.mul;
                let rhs = &dr.div * &v + &dr.rem;
                assert_eq!(lhs, rhs);
            }
        }
    }

    #[test]
    fn test_gcd() {
        let u = Polynomial::from_coefficients(vec![-5i64, 2, 8, -3, -3, 0, 1, 0, 1]);
        let v = Polynomial::from_coefficients(vec![21, -9, -4, 0, 5, 0, 3]);
        let g = gcd(&u, &v);
        assert!(g.coeffs() == &[1] || g.coeffs() == &[-1]);
    }

    #[test]
    fn test_extended_gcd() {
        for u in &[
            vec![-3, -3, 0, 1, 1, 1, 3],
            vec![1, 2, 3, 0, 0, 2, 3, -1],
            vec![0, 0, 0, 0, 0, 0, 0, 1],
        ] {
            for v in &[
                vec![0, 5, 0, 3],
                vec![0, 0, 0, 0, 1],
                vec![-2, -3, 0, 4],
                vec![1, 1, 1, 1],
            ] {
                let u = Polynomial::from_coefficients(u
                                                      .iter()
                                                      .cloned()
                                                      .map(|val| Int::from(val))
                                                      .collect());
                let v = Polynomial::from_coefficients(v
                                                      .iter()
                                                      .cloned()
                                                      .map(|val| Int::from(val))
                                                      .collect());
                let egcd = extended_gcd(u.clone(), v.clone());
                assert_eq!(&u * &egcd.a + &v * &egcd.b, egcd.gcd);
            }
        }
    }
}
