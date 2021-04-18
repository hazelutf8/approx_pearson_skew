#![cfg_attr(not(test), no_std)]
#![forbid(unsafe_code)]
#![deny(missing_docs)]
//! Pearson skew second coeffiecient, using mean and median assuming immutable byte slice
//!
//! References:
//! - [Wolfram Skewness](https://mathworld.wolfram.com/PearsonsSkewnessCoefficients.html)
//! - [Wikipedia Standard Deviation (Population)](https://en.wikipedia.org/wiki/Standard_deviation#Uncorrected_sample_standard_deviation)
//!
//! Usable on no_std due to use of approximate square root from [micromath](https://github.com/tarcieri/micromath)

use core::cmp::Ordering;
use micromath::F32Ext;

/// Slice word size
type Word = u8;
/// Real fraction type
type Rational = f32;

/// Mean/Average value of byte in slice
///
/// ```
/// # use crate::approx_pearson_skew::*;
/// let arr = [1, 1, 5, 6];
/// let avg = mean(&arr).unwrap();
/// assert_eq!(avg, 3.25);
/// ```
pub fn mean(slice: &[Word]) -> Option<Rational> {
    if slice.is_empty() {
        return None;
    }
    // For each elm, sum values, then divide by number of elements
    let avg = slice
        .iter()
        .fold(0.0 as Rational, |acc, &elm| acc + (elm as Rational))
        / (slice.len() as Rational);
    Some(avg)
}

/// Determine the next minimum index and count of occurances
///
/// If given a previous minimum index, find the first index and count occurances.
/// Returns None if slice length zero, or if previous value was max.
///
/// ```
/// # use crate::approx_pearson_skew::*;
/// let arr = [0, 2, 5, 7, 2, 1];
/// let found = next_min(&arr, Some(&1));
/// let (ind, occurance, value) = (1, 2, 2);
/// assert_eq!(found, Some((ind, occurance)));
/// assert_eq!(arr[ind], value);
/// ```
pub fn next_min(slice: &[Word], prev: Option<&Word>) -> Option<(usize, usize)> {
    if slice.is_empty() {
        return None;
    }
    let mut v_in = Word::max_value(); // Inclusive (default) minimum value
    let mut v_ind = 0usize; // Index of first found, always valid if Some(_) returned
    let mut c = 0usize; // Count of found value instances
    if let Some(l_ex) = prev {
        let mut v_found = false;
        for (i, e) in slice.iter().enumerate() {
            match (e.cmp(l_ex), e.cmp(&v_in)) {
                (Ordering::Greater, Ordering::Less) => {
                    v_in = *e;
                    v_ind = i;
                    c = 1usize;
                    v_found = true;
                }
                (_, Ordering::Equal) => {
                    c += 1usize;
                    if !v_found {
                        v_in = *e;
                        v_ind = i;
                        v_found = true;
                    }
                }
                _ => {}
            }
        }
        // If the "found minimum" == "previous minimum", previous was max
        // Moves the `e.cmp(l_ex) == Ordering::Greater` for `c+=1` outside the loop
        if l_ex == &v_in {
            c = 0usize;
        }
    } else {
        for (i, e) in slice.iter().enumerate() {
            match e.cmp(&v_in) {
                Ordering::Less => {
                    v_in = *e;
                    v_ind = i;
                    c = 1usize;
                }
                Ordering::Equal => {
                    c += 1usize;
                }
                _ => {}
            }
        }
    }
    // Case 0: `c == 0`                      -    None
    // Case 1: v_ind changed                 -    Valid `v_ind` and `c` written
    // Case 2: v_ind==0 && prev.is_none()    -    All values are <_>::max_value(), default of 0usize is valid
    // Case 3: v_ind==0 && prev.is_some()    -    Bool `v_found` handles case where remaining values are <_>::max_value()
    if c > 0 {
        Some((v_ind, c))
    } else {
        None
    }
}

/// Find kth value index from unsorted immutable slice
///
/// Using a O(1) size, but O(nk) time method, worst case O(n^2).
///
/// ```
/// # use crate::approx_pearson_skew::*;
/// let arr = [0, 2, 5, 7, 2, 1];
/// let ind = kth_ind(&arr, 3).unwrap();
/// assert_eq!(arr[ind], 2);
/// ```
pub fn kth_ind(slice: &[Word], k: usize) -> Option<usize> {
    if k >= slice.len() {
        return None;
    }
    let mut items = 0usize;
    let mut v_ind = Option::<usize>::None;
    while items <= k {
        let prev = {
            if let Some(v_i) = v_ind {
                Some(&slice[v_i])
            } else {
                None
            }
        };
        if let Some((ind, count)) = next_min(slice, prev) {
            items += count;
            v_ind = Some(ind);
        }
    }
    v_ind
}

/// Find median from unsorted immutable slice
///
/// Using a O(1) size, but O(nk) time method, worst case O(n^2) time.
///
/// ```
/// # use crate::approx_pearson_skew::*;
/// let arr = [1, 2, 6, 7, 6, 1];
/// let med = median(&arr).unwrap();
/// assert_eq!(med, 4.0);
/// ```
pub fn median(slice: &[Word]) -> Option<Rational> {
    if slice.is_empty() {
        return None;
    }

    // kth value of median
    let k = slice.len() / 2;

    // Find kth and possibly k-1 values in unsorted immutable slice
    // Alternately could call `kth_ind()` twice, but that would be slower
    let mut items = 0usize;
    let mut v_items = 0usize;
    let mut v_ind = Option::<usize>::None;
    let mut p_ind = Option::<usize>::None;
    while items <= k {
        p_ind = v_ind;
        let prev = {
            if let Some(v_i) = v_ind {
                Some(&slice[v_i])
            } else {
                None
            }
        };
        if let Some((ind, count)) = next_min(slice, prev) {
            v_items = count;
            items += count;
            v_ind = Some(ind);
        }
    }

    // First of possibly two middle values
    let mut med = slice[v_ind.unwrap()] as Rational;

    // Even slice length, need to average two middle values
    if (slice.len() % 2) == 0 {
        // Determine items covered by previous biggest value
        let p_total_items = items - v_items;
        let k_l = k - 1;

        // If previous biggest value is k-1
        if p_total_items > k_l {
            let k_l_v = slice[p_ind.unwrap()];
            med += k_l_v as Rational;
            med /= 2_f32;
        } else {
            // Both median middle (k and k-1) parts are the same value
            // No need to average identical values
        }
    }
    Some(med)
}

/// Immutable slice population standard deviation
///
/// The mean/average argument allows for value reuse if already known.
/// Uses the approximate square root to allow `no_std` use.
///
/// ```
/// # use crate::approx_pearson_skew::*;
/// let arr = [0, 0, 0, 5, 10];
/// let avg = mean(&arr).unwrap();
/// let std = std_dev_pop(&avg, &arr).unwrap();
/// assert_eq!(std, 4.0);
/// ```
pub fn std_dev_pop(avg: &Rational, slice: &[Word]) -> Option<Rational> {
    if slice.is_empty() {
        return None;
    }
    // Summation of (x_n - avg)^2 for all n elements in the slice
    let sq_sum = slice.iter().fold(0.0 as Rational, |acc, &elm| {
        let e = elm as Rational;
        let delta = e - avg;
        acc + (delta * delta)
    });
    let norm_sq_sum = sq_sum / (slice.len() as Rational);
    Some(F32Ext::sqrt(norm_sq_sum))
}

/// Pearson second skew coefficient, based on the difference between the average and median
///
/// Assumes immutable unsorted slice and uses approximate square root for `no_std` use.
/// Algorithm used for median optimized for size, not time complexity.
///
/// ```
/// # use crate::approx_pearson_skew::*;
/// let arr = [0, 0, 0, 5, 10];
/// let std = pearson_skew_median(&arr).unwrap();
/// assert_eq!(std, 2.25);
/// ```
pub fn pearson_skew_median(slice: &[Word]) -> Option<Rational> {
    let avg = mean(slice)?;
    let med = median(slice)?;
    let std = std_dev_pop(&avg, slice)?;
    Some((3.0 * (avg - med)) / std)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn next_min_check() {
        let test_vec = [
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 1, 1],
            [8, 1, 4, 5, 6, 3, 2, 7],
            [7, 4, 6, 7, 2, 3, 2, 2],
        ];

        // idx 0
        {
            assert_eq!(next_min(&test_vec[0], None), Some((0, 8)));
            assert_eq!(next_min(&test_vec[0], Some(&0)), None);
            assert_eq!(next_min(&test_vec[0], Some(&1)), None);
        }

        // idx 1
        {
            assert_eq!(next_min(&test_vec[1], None), Some((0, 4)));
            assert_eq!(next_min(&test_vec[1], Some(&0)), Some((1, 4)));
            assert_eq!(next_min(&test_vec[1], Some(&1)), None);
            assert_eq!(next_min(&test_vec[1], Some(&2)), None);
        }

        // idx 2
        {
            assert_eq!(next_min(&test_vec[2], None), Some((1, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&0)), Some((1, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&1)), Some((6, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&2)), Some((5, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&3)), Some((2, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&4)), Some((3, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&5)), Some((4, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&6)), Some((7, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&7)), Some((0, 1)));
            assert_eq!(next_min(&test_vec[2], Some(&8)), None);
            assert_eq!(next_min(&test_vec[2], Some(&9)), None);
        }

        //idx 3
        {
            assert_eq!(next_min(&test_vec[3], None), Some((4, 3)));
            assert_eq!(next_min(&test_vec[3], Some(&0)), Some((4, 3)));
            assert_eq!(next_min(&test_vec[3], Some(&1)), Some((4, 3)));
            assert_eq!(next_min(&test_vec[3], Some(&2)), Some((5, 1)));
            assert_eq!(next_min(&test_vec[3], Some(&3)), Some((1, 1)));
            assert_eq!(next_min(&test_vec[3], Some(&4)), Some((2, 1)));
            assert_eq!(next_min(&test_vec[3], Some(&5)), Some((2, 1)));
            assert_eq!(next_min(&test_vec[3], Some(&6)), Some((0, 2)));
            assert_eq!(next_min(&test_vec[3], Some(&7)), None);
            assert_eq!(next_min(&test_vec[3], Some(&8)), None);
        }
    }

    #[test]
    fn kth_ind_check() {
        let test_vec = [
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 1, 1],
            [8, 1, 4, 5, 6, 3, 2, 7],
            [7, 4, 6, 7, 2, 3, 2, 2],
        ];

        // If all values are the same, 0th index should be returned
        for k in 0..test_vec[0].len() {
            let ind = kth_ind(&test_vec[0], k).unwrap();
            assert_eq!(ind, 0);
        }

        // First four are 0s, second four are 1s
        for k in 0..test_vec[1].len() {
            let ind = kth_ind(&test_vec[1], k).unwrap();
            dbg!(k, ind);
            if k < 4 {
                assert_eq!(ind, 0);
            } else {
                assert_eq!(ind, 1);
            }
        }

        // Index of kth unique values in an unsorted slice
        {
            assert_eq!(kth_ind(&test_vec[2], 0), Some(1));
            assert_eq!(kth_ind(&test_vec[2], 1), Some(6));
            assert_eq!(kth_ind(&test_vec[2], 2), Some(5));
            assert_eq!(kth_ind(&test_vec[2], 3), Some(2));
            assert_eq!(kth_ind(&test_vec[2], 4), Some(3));
            assert_eq!(kth_ind(&test_vec[2], 5), Some(4));
            assert_eq!(kth_ind(&test_vec[2], 6), Some(7));
            assert_eq!(kth_ind(&test_vec[2], 7), Some(0));
            assert_eq!(kth_ind(&test_vec[2], 8), None);
            assert_eq!(kth_ind(&test_vec[2], 9), None);
        }

        // Index of kth values with duplicates in an unsorted slice
        // Index is first instance of this value in the slice
        {
            assert_eq!(kth_ind(&test_vec[3], 0), Some(4));
            assert_eq!(kth_ind(&test_vec[3], 1), Some(4));
            assert_eq!(kth_ind(&test_vec[3], 2), Some(4));
            assert_eq!(kth_ind(&test_vec[3], 3), Some(5));
            assert_eq!(kth_ind(&test_vec[3], 4), Some(1));
            assert_eq!(kth_ind(&test_vec[3], 5), Some(2));
            assert_eq!(kth_ind(&test_vec[3], 6), Some(0));
            assert_eq!(kth_ind(&test_vec[3], 7), Some(0));
            assert_eq!(kth_ind(&test_vec[3], 8), None);
            assert_eq!(kth_ind(&test_vec[3], 9), None);
        }
    }

    /// Mutable slice median
    ///
    /// Not the most efficient, but used to check immutable slice median answer
    ///
    /// ```
    /// # use crate::approx_pearson_skew::*;
    /// let mut arr = [1, 0, 6, 0, 5];
    /// let med = mean(mut &arr);
    /// assert_eq!(med, 1.0);
    /// ```
    fn median_mut(slice: &mut [Word]) -> Option<Rational> {
        if slice.is_empty() {
            return None;
        }
        slice.sort();
        let i = slice.len() / 2;
        let mut med = slice[i] as Rational;
        if (slice.len() % 2) == 0 {
            med = (med + (slice[i - 1] as Rational)) / 2.0;
        }
        Some(med)
    }

    #[test]
    fn median_odd_check() {
        let mut test_vec = [
            [0, 0, 0, 0, 0, 0, 0, 0, 1],
            [1, 1, 0, 1, 1, 0, 0, 0, 1],
            [9, 1, 8, 2, 7, 3, 6, 3, 5],
            [1, 1, 1, 1, 7, 3, 6, 3, 5],
            [9, 1, 8, 2, 7, 0, 0, 0, 0],
        ];
        for i in 0..test_vec.len() {
            let v = &mut test_vec[i];
            assert_eq!(v.len() % 2, 1);
            dbg!(i, &v);
            let a = median(v);
            let b = median_mut(v);
            assert_eq!(a, b);
        }
    }

    #[test]
    fn median_even_check() {
        let mut test_vec = [
            [0, 0, 0, 0, 0, 0, 0, 1],
            [1, 1, 0, 1, 0, 0, 0, 1],
            [9, 1, 8, 2, 3, 6, 3, 5],
            [1, 1, 1, 1, 3, 6, 3, 5],
            [9, 1, 8, 2, 0, 0, 0, 0],
        ];
        for i in 0..test_vec.len() {
            let v = &mut test_vec[i];
            assert_eq!(v.len() % 2, 0);
            dbg!(i, &v);
            let a = median(v);
            let b = median_mut(v);
            assert_eq!(a, b);
        }
    }

    #[test]
    fn skew_check() {
        // Python: Avg 3.000000, Median 0.000000, StdDevPop 4.000000, Skew 2.250000
        let arr = [0, 0, 0, 5, 10];
        let skew = pearson_skew_median(&arr).unwrap();
        assert_eq!(skew, 2.25)
    }
}
