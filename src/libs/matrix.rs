//! A *symmetric* scoring matrix to be used for clustering.

use std::collections::HashMap;

#[derive(Debug)]
pub struct ScoringMatrix<T> {
    size: usize,
    same: T,
    missing: T,
    data: HashMap<(usize, usize), T>,
}

impl<T> ScoringMatrix<T>
where
    T: Default + Copy,
{
    /// Constructs a new symmetric scoring matrix, setting the size and the default values
    ///
    /// # Examples
    ///
    /// ```
    /// # use hnsm::ScoringMatrix;
    /// let m: ScoringMatrix<i32> = ScoringMatrix::new(5, 0, 1);
    /// ```
    pub fn new(size: usize, same: T, missing: T) -> Self {
        let data: HashMap<(usize, usize), T> = HashMap::new();
        ScoringMatrix {
            size,
            same,
            missing,
            data,
        }
    }

    /// Returns the number of rows and columns of the matrix.
    ///
    /// ```
    /// # use hnsm::ScoringMatrix;
    /// let m: ScoringMatrix<i32> = ScoringMatrix::new(5, 0, 1);
    /// assert_eq!(m.size(), 5);
    /// ```
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns the value of the given cell.
    ///
    /// ```
    /// # use hnsm::ScoringMatrix;
    /// let mut m = ScoringMatrix::new(5, 0, 1);
    /// m.set(1, 2, 42);
    /// assert_eq!(m.get(1, 2), 42);
    /// assert_eq!(m.get(2, 1), 42);
    /// assert_eq!(m.get(3, 3), 0);
    /// assert_eq!(m.get(1, 3), 1);
    /// ```
    pub fn get(&self, row: usize, col: usize) -> T {
        if row < col {
            *self.data.get(&(row, col)).unwrap_or(&self.missing)
        } else if row == col {
            *self.data.get(&(row, col)).unwrap_or(&self.same)
        }
        else {
            *self.data.get(&(col, row)).unwrap_or(&self.missing)
        }
    }

    /// Sets the value of the given cell.
    pub fn set(&mut self, row: usize, col: usize, value: T) {
        if row <= col {
            self.data.insert((row, col), value);
        } else {
            self.data.insert((col, row), value);
        }
    }
}
