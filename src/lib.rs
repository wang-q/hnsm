// #![feature(array_chunks)]
// #![feature(slice_as_chunks)]
// Add these imports to use the stdsimd library
#![feature(portable_simd)]

pub mod libs;

pub use crate::libs::hash::*;
