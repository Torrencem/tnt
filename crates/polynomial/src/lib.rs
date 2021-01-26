
#![feature(iterator_fold_self)]

use std::fmt::Write;
use std::ops::{Add, Sub, Mul, Div, Neg, AddAssign, MulAssign, SubAssign, DivAssign};
use std::cmp::max;
use std::borrow::Cow;
use tnt_alg::*;

// needs to be exported for macros
#[doc(hidden)]
pub use num_traits::{Zero, One};

/// A polynomial with coefficients in some type T.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct Polynomial<T> {
    // Coefficients of the polynomial, in order from lowest power of x to highest. Invariant:
    // there should never be trailing zeroes in this vector (so that the degree is always the
    // length of the vector
    cs: Vec<T>,
}

// Impl display strictly. It would be nice to have fallbacks that don't require so much (maybe with
// specializations on different trait bounds) but this isn't easy yet.
// It's also annoying that this doesn't support nested polynomials! But again, specialization might
// be the only solution.
impl<T> Polynomial<T>
where T: Zero + One + fmt::Display + Clone + PartialEq,
      for<'a> &'a T: Neg<Output=T>,
      for<'a> &'a T: PartialOrd {
    /// Format this polynomial with a given variable name `var`. This is meant to be used in a
    /// similar way to fmt::Display::fmt(..).
    pub fn format_with_var<C: fmt::Display, W: fmt::Write>(&self, var: C, with_parens: bool, f: &mut W) -> fmt::Result {
        if self.cs.len() == 0 {
            let z: T = Zero::zero();
            return write!(f, "{}", z);
        }
        let mut leading = true;
        for (index, coeff) in self.cs.iter()
                        .enumerate()
                        .rev()
                        .filter(|(_, val)| !val.is_zero())
        {
            let negative = coeff < &Zero::zero();
            let coeff = if negative && !leading { Cow::Owned(-coeff) } else { Cow::Borrowed(coeff) };
            let lead = if !leading { 
                if negative { " - " } else { " + " }
            } else { "" };
            let lparen = if with_parens { "(" } else { "" };
            let rparen = if with_parens { ")" } else { "" };
            if coeff.is_one() {
                if index == 0 {
                    write!(f, "{}{}{}{}", lead, lparen, coeff, rparen)?;
                } else if index == 1 {
                    write!(f, "{}{}", lead, var)?;
                } else {
                    write!(f, "{}{}^{}", lead, var, index)?;
                }
            } else {
                if index == 0 {
                    write!(f, "{}{}{}{}", lead, lparen, coeff, rparen)?;
                } else if index == 1 {
                    write!(f, "{}{}{}{}{}", lead, lparen, coeff, rparen, var)?;
                } else {
                    write!(f, "{}{}{}{}{}^{}", lead, lparen, coeff, rparen, var, index)?;
                }
            }
            leading = false;
        }
        Ok(())
    }
}

use std::fmt;
impl<T> fmt::Display for Polynomial<T>
where T: Zero + One + fmt::Display + Clone + PartialEq,
      for<'a> &'a T: Neg<Output=T>,
      for<'a> &'a T: PartialOrd {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.format_with_var("x", false, f)
    }
}

impl<T> Polynomial<T> {
    /// Returns a reference to the coefficients of the polynomial. The i-th index of this slice is
    /// the coefficient on x^i.
    pub fn coeffs(&self) -> &[T] {
        &self.cs
    }

    /// Apply a function to each of the coefficients of the polynomial.
    pub fn map_coeffs<F>(&mut self, f: F)
    where F: FnMut(T) -> T,
          T: Zero {
        let new_cs = self.cs.drain(..)
            .map(f)
            .collect();
        *self = Polynomial::from_coefficients(new_cs);
    }

    /// Apply a function to each of the coefficients of the polynomial. This also passes the index
    /// of the coefficient to the function (the index is the same as the exponent on x).
    pub fn map_enumerate_coeffs<F>(&mut self, mut f: F)
    where F: FnMut(usize, T) -> T,
          T: Zero {
        let new_cs = self.cs.drain(..)
            .enumerate()
            .map(|(i, c)| f(i, c))
            .collect();
        *self = Polynomial::from_coefficients(new_cs);
    }

    /// Get the degree of the polynomial.
    pub fn degree(&self) -> usize {
        self.cs.len().saturating_sub(1)
    }

    /// Get the leading coefficient of the polynomial. Note that this returns either a reference
    /// to the leading coefficient, or an owned zero value if this polynomial is zero, so you will
    /// need to call .lc().into_owned() to get a T or .lc().as_ref() to get a &T.
    pub fn lc<'a>(&'a self) -> Cow<'a, T> 
    where T: Clone + Zero
    {
        if self.cs.len() == 0 {
            Cow::Owned(Zero::zero())
        } else {
            Cow::Borrowed(&self.cs[self.cs.len() - 1])
        }
    }

    /// Get the content of the polynomial. The content is the gcd of all of the coefficients.
    pub fn cont(&self) -> T
    where T: Gcd + Zero + Clone {
        self.cs.iter()
            .cloned()
            .fold_first(|a, b| {
                gcd(&a, &b)
            })
            .unwrap_or_else(|| Zero::zero())
    }

    /// Reduce the coefficients of the polynomial by the polynomial's content. Given a polynomial
    /// p(x), this computes p(x) / cont(p(x)).
    pub fn pp<U>(&self) -> Polynomial<T>
    where T: Gcd + Zero + Clone + PseudoDivRem<Output=T, MultType=U>
    {
        let c = self.cont();
        let mut p = self.clone();
        for coeff in p.cs.iter_mut() {
            let val = pdiv(coeff.clone(), c.clone());
            *coeff = val;
        }
        p
    }
}

impl<T: Zero> Polynomial<T> {
    /// Create a new polynomial from a Vec of coefficients. In the vec, index i should be the
    /// coefficient on x^i.
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

    pub(crate) fn fix_coefficients(&mut self) {
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
        // TODO: This function could be a lot smarter. It could even use specialization to print
        // things better with One and Zero etc.
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
}

impl<T> Polynomial<T> {
    /// Evaluate the polynomial at a specific value.
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

// Implementations for alg traits
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
/// pseudo_divrem, except it doesn't compute mul, which prevents overflows sometimes.
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
      for<'c> &'c Polynomial<T>: PseudoDivRem<Output=Polynomial<T>, MultType=T>,
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


// Utility functions for polynomial macro
#[doc(hidden)]
pub const fn const_max__(a: usize, b: usize) -> usize {
    if a < b {
        b
    } else {
        a
    }
}

#[doc(hidden)]
#[macro_export]
macro_rules! max__ {
    ($e:expr, $($es:expr),+) => {
        $crate::const_max__($e, $crate::max__!($($es),+))
    };
    ($e:expr) => {
        $e
    };
}

/// Create a polynomial from coefficients. Note that the syntax of this macro includes semicolons
/// ";" after each coefficient expression, since otherwise syntax would be very ambiguous. In
/// addition, each exponent must be a valid usize literal.
///
/// # Examples
///
/// ```
/// // Full way to use this macro
/// use tnt_polynomial::{Polynomial, poly};
/// let val = 10;
/// let p = poly!(val*2;x^2 + 10;x^1 + -1;x^0);
/// assert_eq!(p, Polynomial::from_coefficients(vec![-1, 10, 20]));
///
/// // All possible abbreviated shortcuts for convenience
/// let p = poly!(x^5);
/// assert_eq!(p, Polynomial::from_coefficients(vec![0, 0, 0, 0, 0, 1]));
///
/// let p = poly!(x);
/// assert_eq!(p, Polynomial::from_coefficients(vec![0, 1]));
///
/// let p = poly!(3;x);
/// assert_eq!(p, Polynomial::from_coefficients(vec![0, 3]));
///
/// let p = poly!(100);
/// assert_eq!(p, Polynomial::from_coefficients(vec![100]));
/// ```
#[macro_export]
macro_rules! poly {
    ($a:expr ;x^ $b:literal $(+ $c:expr ;x^ $d:literal)*) => {{
        const SIZE: usize = $crate::max__!($b $(, $d)*) + 1;
        let mut v = vec![$crate::Zero::zero(); SIZE];
        v[$b] = $a;
        $(
            v[$d] = $c;
        )*
        Polynomial::from_coefficients(v)
    }};
    (x^$b:literal) => {{
        let mut v = vec![$crate::Zero::zero(); $b + 1];
        v[$b] = 1;
        Polynomial::from_coefficients(v)
    }};
    ($a:expr ;x) => {{
        Polynomial::from_coefficients(vec![$crate::Zero::zero(), $a])
    }};
    (x) => {{
        Polynomial::from_coefficients(vec![$crate::Zero::zero(), $crate::One::one()])
    }};
    ($a:expr) => {{
        Polynomial::from_coefficients(vec![$a])
    }};
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

    #[test]
    fn test_fmt() {
        let p = Polynomial::from_coefficients(vec![3, 2, 1]);
        assert_eq!(p.to_string(), "x^2 + 2x + 3");

        let p = Polynomial::from_coefficients(vec![1, 1, 2, 0, 3]);
        assert_eq!(p.to_string(), "3x^4 + 2x^2 + x + 1");
    }
}
