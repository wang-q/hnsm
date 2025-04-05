use std::fs::File;
use std::io::{Read, Write};

//----------------------------
// AsmEntry
//----------------------------
#[derive(Default, Clone)]
pub struct AsmEntry {
    name: String,
    list: Vec<f32>,
}

impl AsmEntry {
    // Immutable accessors
    pub fn name(&self) -> &String {
        &self.name
    }
    pub fn list(&self) -> &Vec<f32> {
        &self.list
    }

    pub fn new() -> Self {
        Self {
            name: String::new(),
            list: vec![],
        }
    }

    /// Constructed from range and seq
    ///
    /// ```
    /// # use hnsm::AsmEntry;
    /// let name = "Es_coli_005008_GCF_013426115_1".to_string();
    /// let list : Vec<f32> = vec![1.0,5.0,2.0,7.0,6.0,6.0];
    /// let entry = AsmEntry::from(&name, &list);
    /// # assert_eq!(*entry.name(), "Es_coli_005008_GCF_013426115_1");
    /// # assert_eq!(*entry.list().get(1).unwrap(), 5f32);
    /// ```
    pub fn from(name: &str, vector: &[f32]) -> Self {
        Self {
            name: name.to_owned(),
            list: Vec::from(vector),
        }
    }

    /// ```
    /// # use hnsm::AsmEntry;
    /// let line = "Es_coli_005008_GCF_013426115_1\t1,5,2,7,6,6".to_string();
    /// let entry = AsmEntry::parse(&line);
    /// # assert_eq!(*entry.name(), "Es_coli_005008_GCF_013426115_1");
    /// # assert_eq!(*entry.list().get(1).unwrap(), 5f32);
    /// ```
    pub fn parse(line: &str) -> AsmEntry {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            let name = fields[0].to_string();
            let parts: Vec<&str> = fields[1].split(',').collect();
            let list: Vec<f32> = parts.iter().map(|e| e.parse::<f32>().unwrap()).collect();
            Self::from(&name, &list)
        } else {
            Self::new()
        }
    }
}

impl std::fmt::Display for AsmEntry {
    /// To string
    ///
    /// ```
    /// # use hnsm::AsmEntry;
    /// let name = "Es_coli_005008_GCF_013426115_1".to_string();
    /// let list : Vec<f32> = vec![1.0,5.0,2.0,7.0,6.0,6.0];
    /// let entry = AsmEntry::from(&name, &list);
    /// assert_eq!(entry.to_string(), "Es_coli_005008_GCF_013426115_1\t1,5,2,7,6,6\n");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\n",
            self.name(),
            self.list
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<_>>()
                .join(","),
        )?;
        Ok(())
    }
}

//----------------------------
// Seq types
//----------------------------
pub fn is_fq(input: &str) -> bool {
    let path = std::path::Path::new(input);

    // Create a buffer to store the first two bytes
    let mut buffer = [0; 2];
    {
        let mut file = match std::fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };
        file.read_exact(&mut buffer).unwrap();
    }

    // Check if the file is in Gzip format
    let is_fq;
    if buffer[0] == 0x1f && buffer[1] == 0x8b {
        let mut decoder = flate2::read::GzDecoder::new(File::open(path).unwrap());
        let mut buffer = [0; 2]; // Recreate the buffer
        decoder.read_exact(&mut buffer).unwrap();

        // Determine the format of the decompressed file
        match buffer[0] as char {
            '>' => is_fq = false,
            '@' => is_fq = true,
            _ => unreachable!("Unknown file format"),
        }
    } else {
        // The file is in plain text format, determine the format
        match buffer[0] as char {
            '>' => is_fq = false,
            '@' => is_fq = true,
            _ => unreachable!("Unknown file format"),
        }
    }

    is_fq
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use tempfile::tempdir;

    #[test]
    fn test_is_fq_plain_text() {
        let dir = tempdir().unwrap();

        // Create a plain text FASTQ file
        let fq_file_path = dir.path().join("test.fq");
        {
            let mut file = File::create(&fq_file_path).unwrap();
            writeln!(file, "@SEQ_ID").unwrap(); // FASTQ format
        }
        assert!(is_fq(fq_file_path.to_str().unwrap()));

        // Create a plain text FASTA file
        let fasta_file_path = dir.path().join("test.fasta");
        {
            let mut file = File::create(&fasta_file_path).unwrap();
            writeln!(file, ">SEQ_ID").unwrap(); // FASTA format
        }
        assert!(!is_fq(fasta_file_path.to_str().unwrap()));
    }

    #[test]
    fn test_is_fq_gzip() {
        let dir = tempdir().unwrap();

        // Create a Gzip FASTQ file
        let fq_file_path = dir.path().join("test.fq.gz");
        {
            let file = File::create(&fq_file_path).unwrap();
            let mut encoder = GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, "@SEQ_ID").unwrap(); // FASTQ format
            encoder.finish().unwrap();
        }
        assert!(is_fq(fq_file_path.to_str().unwrap()));

        // Create a Gzip FASTA file
        let fasta_file_path = dir.path().join("test.fasta.gz");
        {
            let file = File::create(&fasta_file_path).unwrap();
            let mut encoder = GzEncoder::new(file, flate2::Compression::default());
            writeln!(encoder, ">SEQ_ID").unwrap(); // FASTA format
            encoder.finish().unwrap();
        }
        assert!(!is_fq(fasta_file_path.to_str().unwrap()));
    }
}

pub fn pause() {
    let mut stdin = std::io::stdin();
    let mut stdout = std::io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}
