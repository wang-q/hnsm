#![feature(array_chunks)]
#![feature(slice_as_chunks)]
// Add these imports to use the stdsimd library
#![feature(portable_simd)]

#[macro_use]
extern crate lazy_static;

pub mod libs;

pub use crate::libs::alignment::*;
pub use crate::libs::dbscan::Dbscan;
pub use crate::libs::fas::*;
pub use crate::libs::hash::*;
pub use crate::libs::io::*;
pub use crate::libs::linalg::*;
pub use crate::libs::loc::*;
pub use crate::libs::matrix::ScoringMatrix;
pub use crate::libs::nt::*;
