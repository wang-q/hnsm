/// Standard IUB/IUPAC Nucleic Acid Codes
/// ```text
/// Code =>  Nucleic Acid(s)
///  A   =>  Adenine
///  C   =>  Cytosine
///  G   =>  Guanine
///  T   =>  Thymine
///  U   =>  Uracil
///  M   =>  A or C (amino)
///  R   =>  A or G (purine)
///  W   =>  A or T (weak)
///  S   =>  C or G (strong)
///  Y   =>  C or T (pyrimidine)
///  K   =>  G or T (keto)
///  V   =>  A or C or G
///  H   =>  A or C or T
///  D   =>  A or G or T
///  B   =>  C or G or T
///  N   =>  A or G or C or T (any)
/// ```
#[allow(dead_code)]
#[repr(usize)]
#[derive(Clone, Copy)]
pub enum Nt {
    A = 0,
    C = 1,
    G = 2,
    T = 3, // U
    N = 4,
    Invalid = 255,
}

#[allow(dead_code)]
impl Nt {
    pub const U: Nt = Nt::T;
}

/// Maps an ASCII chars to index
/// ```text
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
/// N => 4
/// Invalid => 255
/// ```
pub static NT_VAL: &[usize; 256] = &{
    let mut array = [255; 256];

    array[b'A' as usize] = 0;
    array[b'a' as usize] = 0;

    array[b'C' as usize] = 1;
    array[b'c' as usize] = 1;

    array[b'G' as usize] = 2;
    array[b'g' as usize] = 2;

    array[b'T' as usize] = 3;
    array[b't' as usize] = 3;
    array[b'U' as usize] = 3;
    array[b'u' as usize] = 3;

    array[b'M' as usize] = 4;
    array[b'm' as usize] = 4;
    array[b'R' as usize] = 4;
    array[b'r' as usize] = 4;
    array[b'W' as usize] = 4;
    array[b'w' as usize] = 4;
    array[b'S' as usize] = 4;
    array[b's' as usize] = 4;
    array[b'Y' as usize] = 4;
    array[b'y' as usize] = 4;
    array[b'K' as usize] = 4;
    array[b'k' as usize] = 4;
    array[b'V' as usize] = 4;
    array[b'v' as usize] = 4;
    array[b'H' as usize] = 4;
    array[b'h' as usize] = 4;
    array[b'D' as usize] = 4;
    array[b'd' as usize] = 4;
    array[b'B' as usize] = 4;
    array[b'b' as usize] = 4;
    array[b'N' as usize] = 4;
    array[b'n' as usize] = 4;

    array
};

pub fn to_nt(nt: u8) -> Nt {
    match NT_VAL[nt as usize] {
        0 => Nt::A,
        1 => Nt::C,
        2 => Nt::G,
        3 => Nt::T,
        4 => Nt::N,
        _ => Nt::Invalid,
    }
}

/// ```
/// assert!(hnsm::is_n(b'n'));
/// assert!(!hnsm::is_n(b'Z'));
/// assert!(!hnsm::is_n(b'-'));
/// assert!(!hnsm::is_n(b'A'));
/// ```
pub fn is_n(nt: u8) -> bool {
    NT_VAL[nt as usize] == Nt::N as usize
}

/// ```
/// assert!(hnsm::is_lower(b'n'));
/// assert!(!hnsm::is_lower(b'A'));
/// ```
pub fn is_lower(nt: u8) -> bool {
    (nt as char).is_ascii_lowercase()
}

/// ```
/// let dna = b"NCTAGTCGTATCGTAGCTAGNNC";
/// assert_eq!(hnsm::count_n(dna), 3);
/// ```
pub fn count_n(seq: &[u8]) -> usize {
    let mut n_cnt = 0;
    for c in seq {
        if is_n(*c) {
            n_cnt += 1;
        }
    }

    n_cnt
}

/// convert IUPAC ambiguous codes to 'N'
/// ```
/// assert_eq!(hnsm::to_n(b'a'), b'a');
/// assert_eq!(hnsm::to_n(b'M'), b'N');
/// ```
pub fn to_n(nt: u8) -> u8 {
    if is_n(nt) {
        b'N'
    } else {
        nt
    }
}

/// A static lookup table for nucleotide complements.
///
/// The table maps each ASCII byte to its complement. For example:
/// - `A` maps to `T`
/// - `C` maps to `G`
/// - `G` maps to `C`
/// - `T` maps to `A`
/// - Lowercase letters are also supported (e.g., `a` maps to `t`).
/// - Non-nucleotide characters (e.g., ` `, `-`) map to themselves.
/// - IUPAC codes are supported (e.g., `M` maps to `K`, `R` maps to `Y`).
///
/// # Examples
///
/// ```
/// // Test the lookup table directly
/// assert_eq!(hnsm::NT_COMP[b'A' as usize], b'T');
/// assert_eq!(hnsm::NT_COMP[b'C' as usize], b'G');
/// assert_eq!(hnsm::NT_COMP[b'G' as usize], b'C');
/// assert_eq!(hnsm::NT_COMP[b'T' as usize], b'A');
/// assert_eq!(hnsm::NT_COMP[b'a' as usize], b't');
/// assert_eq!(hnsm::NT_COMP[b'c' as usize], b'g');
/// assert_eq!(hnsm::NT_COMP[b' ' as usize], b' ');
/// assert_eq!(hnsm::NT_COMP[b'-' as usize], b'-');
///
/// // Test IUPAC codes
/// assert_eq!(hnsm::NT_COMP[b'M' as usize], b'K'); // M (A or C) -> K (T or G)
/// assert_eq!(hnsm::NT_COMP[b'R' as usize], b'Y'); // R (A or G) -> Y (T or C)
/// assert_eq!(hnsm::NT_COMP[b'W' as usize], b'W'); // W (A or T) -> W (T or A)
/// assert_eq!(hnsm::NT_COMP[b'S' as usize], b'S'); // S (C or G) -> S (G or C)
/// assert_eq!(hnsm::NT_COMP[b'Y' as usize], b'R'); // Y (C or T) -> R (G or A)
/// assert_eq!(hnsm::NT_COMP[b'K' as usize], b'M'); // K (G or T) -> M (C or A)
/// assert_eq!(hnsm::NT_COMP[b'V' as usize], b'B'); // V (A or C or G) -> B (T or G or C)
/// assert_eq!(hnsm::NT_COMP[b'H' as usize], b'D'); // H (A or C or T) -> D (T or G or A)
/// assert_eq!(hnsm::NT_COMP[b'D' as usize], b'H'); // D (A or G or T) -> H (T or C or A)
/// assert_eq!(hnsm::NT_COMP[b'B' as usize], b'V'); // B (C or G or T) -> V (G or C or A)
/// assert_eq!(hnsm::NT_COMP[b'N' as usize], b'N'); // N (any base) -> N (any base)
/// assert_eq!(hnsm::NT_COMP[b'U' as usize], b'A'); // U (Uracil) -> A (Adenine)
/// ```
pub static NT_COMP: &[u8; 256] = &{
    let mut array = [255; 256];

    array[b' ' as usize] = b' ';
    array[b'-' as usize] = b'-';

    array[b'A' as usize] = b'T';
    array[b'C' as usize] = b'G';
    array[b'G' as usize] = b'C';
    array[b'T' as usize] = b'A';
    array[b'M' as usize] = b'K';
    array[b'R' as usize] = b'Y';
    array[b'W' as usize] = b'W';
    array[b'S' as usize] = b'S';
    array[b'Y' as usize] = b'R';
    array[b'K' as usize] = b'M';
    array[b'V' as usize] = b'B';
    array[b'H' as usize] = b'D';
    array[b'D' as usize] = b'H';
    array[b'B' as usize] = b'V';
    array[b'N' as usize] = b'N';
    array[b'U' as usize] = b'A';

    array[b'a' as usize] = b't';
    array[b'c' as usize] = b'g';
    array[b'g' as usize] = b'c';
    array[b't' as usize] = b'a';
    array[b'm' as usize] = b'k';
    array[b'r' as usize] = b'y';
    array[b'w' as usize] = b'w';
    array[b's' as usize] = b's';
    array[b'y' as usize] = b'r';
    array[b'k' as usize] = b'm';
    array[b'v' as usize] = b'b';
    array[b'h' as usize] = b'd';
    array[b'd' as usize] = b'h';
    array[b'b' as usize] = b'v';
    array[b'n' as usize] = b'n';
    array[b'u' as usize] = b'a';

    array
};

/// Computes the complement of a nucleotide sequence, returning an iterator.
///
/// # Arguments
///
/// * `seq` - A slice of bytes representing the nucleotide sequence.
///
/// # Examples
///
/// ```
/// // Test complement
/// let seq = b"ACGT";
/// let complemented: Vec<u8> = hnsm::complement(seq).collect();
/// assert_eq!(complemented, b"TGCA");
///
/// let seq_with_iupac = b"MRWSYKVHDBN";
/// let complemented_iupac: Vec<u8> = hnsm::complement(seq_with_iupac).collect();
/// assert_eq!(complemented_iupac, b"KYWSRMBDHVN");
///
/// // Test reverse complement by reversing the input first
/// let rev_complemented: Vec<u8> = hnsm::complement(seq).rev().collect();
/// assert_eq!(rev_complemented, b"ACGT"); // Reverse complement of "ACGT" is "ACGT"
/// ```
pub fn complement<'a>(seq: &'a [u8]) -> impl DoubleEndedIterator<Item = u8> + 'a {
    seq.iter().copied().map(move |b| NT_COMP[b as usize])
}

/// Computes the reverse complement of a nucleotide sequence, returning an iterator.
///
/// # Arguments
///
/// * `seq` - A slice of bytes representing the nucleotide sequence.
///
/// # Examples
///
/// ```
/// let seq = b"ACGT";
/// let rev_complemented: Vec<u8> = hnsm::rev_comp(seq).collect();
/// assert_eq!(rev_complemented, b"ACGT"); // Reverse complement of "ACGT" is "ACGT"
///
/// let seq_with_iupac = b"MRWSYKVHDBN";
/// let rev_complemented_iupac: Vec<u8> = hnsm::rev_comp(seq_with_iupac).collect();
/// assert_eq!(rev_complemented_iupac, b"NVHDBMRSWYK"); // Reverse complement of "MRWSYKVHDBN"
/// ```
pub fn rev_comp<'a>(seq: &'a [u8]) -> impl Iterator<Item = u8> + 'a {
    seq.iter().rev().map(|&b| NT_COMP[b as usize])
}

/// block -> row -> column
pub static AA_TAB: &[[[char; 4]; 4]; 4] = &[
    [
        ['K', 'N', 'K', 'N'], // AAA, AAC, AAG, AAU/AAT
        ['T', 'T', 'T', 'T'], // ACA, ACC, ACG, ACU/ACT
        ['R', 'S', 'R', 'S'], // AGA, AGC, AGG, AGU/AGT
        ['I', 'I', 'M', 'I'], // AUA/ATA, AUC/ATC, AUG/ATG, AUU/ATT
    ],
    [
        ['Q', 'H', 'Q', 'H'], // CAA, CAC, CAG, CAU/CAT
        ['P', 'P', 'P', 'P'], // CCA, CCC, CCG, CCU/CCT
        ['R', 'R', 'R', 'R'], // CGA, CGC, CGG, CGU/CGT
        ['L', 'L', 'L', 'L'], // CUA/CTA, CUC/CTC, CUG/CTG, CUU/CTT
    ],
    [
        ['E', 'D', 'E', 'D'], // GAA, GAC, GAG, GAU/GAT
        ['A', 'A', 'A', 'A'], // GCA, GCC, GCG, GCU/GCT
        ['G', 'G', 'G', 'G'], // GGA, GGC, GGG, GGU/GGT
        ['V', 'V', 'V', 'V'], // GUA/GTA, GUC/GTC, GUG/GTG, GUU/GTT
    ],
    [
        ['*', 'Y', '*', 'Y'], // UAA/TAA, UAC/TAC, UAG/TAG, UAU/TAT
        ['S', 'S', 'S', 'S'], // UCA/TCA, UCC/TCC, UCG/TCG, UCU/TCT
        ['*', 'C', 'W', 'C'], // UGA/TGA, UGC/TGC, UGG/TGG, UGU/TGT
        ['L', 'F', 'L', 'F'], // UUA/TTA, UUC/TTC, UUG/TTG, UUU/TTT
    ],
];

/// ```
/// let dna = b"GCTAGTCGTATCGTAGCTAGTC";
/// assert_eq!(&hnsm::translate(dna), "ASRIVAS");
///
/// let rna = b"GCUAGUCGUAUCGUAGCUAGUC";
/// assert_eq!(&hnsm::translate(rna), "ASRIVAS");
///
/// // To shift the reading frame, pass in a slice
/// assert_eq!(&hnsm::translate(&dna[1..]), "LVVS*LV");
/// assert_eq!(&hnsm::translate(&dna[2..]), "*SYRS*");
/// ```
// https://github.com/dweb0/protein-translate/blob/master/src/lib.rs
pub fn translate(seq: &[u8]) -> String {
    let mut peptide = String::with_capacity(seq.len() / 3);

    'outer: for triplet in seq.chunks_exact(3) {
        for c in triplet {
            if !c.is_ascii() {
                peptide.push('X');
                continue 'outer;
            }
            if NT_VAL[*c as usize] == Nt::N as usize {
                peptide.push('X');
                continue 'outer;
            }
            if NT_VAL[*c as usize] == Nt::Invalid as usize {
                peptide.push('X');
                continue 'outer;
            }
        }

        let c1 = NT_VAL[triplet[0] as usize];
        let c2 = NT_VAL[triplet[1] as usize];
        let c3 = NT_VAL[triplet[2] as usize];

        let amino_acid = AA_TAB[c1][c2][c3];

        peptide.push(amino_acid);
    }
    peptide
}

/// Detect ORFs in a translated protein sequence
///
/// # Examples
///
/// ```
/// let protein = "MGGMGG*AGG";
/// let orfs = hnsm::find_orfs(protein);
/// assert_eq!(orfs, vec![
///     ("MGGMGG*".to_string(), 0, 7),
///     ("AGG".to_string(), 7, 10)
/// ]);
/// ```
pub fn find_orfs(protein: &str) -> Vec<(String, usize, usize)> {
    let mut orfs = Vec::new();
    let mut start = 0;

    while start < protein.len() {
        // Find the start aa (not just M)
        if let Some(orf_start) = protein[start..].find(|c: char| c != 'X' && c != 'x') {
            let orf_start_pos = start + orf_start;
            // Find the stop codon '*'
            if let Some(orf_end) = protein[orf_start_pos..].find('*') {
                let orf_end_pos = orf_start_pos + orf_end;
                let orf_seq = &protein[orf_start_pos..=orf_end_pos];
                // Include the stop codon
                orfs.push((orf_seq.to_string(), orf_start_pos, orf_end_pos + 1));
                start = orf_end_pos + 1; // Continue searching for the next ORF
            } else {
                // If no stop codon is found, treat the remaining sequence as an ORF
                let orf_seq = &protein[orf_start_pos..];
                orfs.push((orf_seq.to_string(), orf_start_pos, protein.len()));
                break;
            }
        } else {
            break; // No start codon found, end the search
        }
    }

    // eprintln!("orfs = {:#?}", orfs);

    orfs
}
