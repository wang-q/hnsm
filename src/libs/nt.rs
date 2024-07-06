/// Standard IUB/IUPAC Nucleic Acid Codes
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
#[allow(dead_code)]
#[repr(usize)]
enum Nt {
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
///
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
/// N => 4
/// Invalid => 255
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
