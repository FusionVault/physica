#![doc = include_str!("../README.md")]

mod gfc;
mod technotes;

pub use gfc::{Gfc, GfcError, julian_epoch, read_gfc};
pub use technotes::{Degree1, SlrC20C30, TechNoteError, read_tn13, read_tn14};
