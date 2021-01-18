
use runt_alg::Zero;
use std::fmt::Write;

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

use std::ops::Add;
use std::cmp::max;

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

use std::ops::Sub;
use std::ops::Neg;

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

use std::ops::{Mul, AddAssign};
// Polynomial Multiplication algorithm

impl<'c, 'd, T: Zero> Mul<&'c Polynomial<T>> for &'d Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      for<'b> T: AddAssign {
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
      for<'b> T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: Polynomial<T>) -> Self::Output {
        (&self) * (&other)
    }
}

impl<'c, T: Zero> Mul<Polynomial<T>> for &'c Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      for<'b> T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: Polynomial<T>) -> Self::Output {
        self * (&other)
    }
}

impl<'d, T: Zero> Mul<&'d Polynomial<T>> for Polynomial<T> 
where for<'a, 'b> &'a T: Mul<&'b T, Output=T>,
      for<'b> T: AddAssign {
    type Output = Polynomial<T>;

    fn mul(self, other: &Polynomial<T>) -> Self::Output {
        (&self) * other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiplication() {
        let p = Polynomial::from_coefficients(vec![1, 1, 1]);
        let q = Polynomial::from_coefficients(vec![1, 1]);
        println!("({}) * ({}) = {}", p.pretty_format("x"), q.pretty_format("x"), (p * q).pretty_format("x"));
    }
}
