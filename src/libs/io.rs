use std::collections::HashMap;
use std::io::BufRead;

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
    pub fn from(name: &String, vector: &[f32]) -> Self {
        Self {
            name: name.clone(),
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

pub fn load_file(infile: &str) -> Vec<((String, String), f32)> {
    let mut pair_scores = Vec::new();
    let reader = intspan::reader(infile);
    for line in reader.lines().map_while(Result::ok) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let n1 = fields[0].to_string();
            let n2 = fields[1].to_string();
            let score: f32 = fields[2].parse::<f32>().unwrap();
            pair_scores.push(((n1, n2), score));
        }
    }
    pair_scores
}

pub fn populate(
    pair_scores: &Vec<((String, String), f32)>,
) -> (crate::ScoringMatrix<f32>, Vec<String>) {
    // Create a mapping from string identifiers to indices
    let mut index_map = HashMap::new();
    let mut index_name = vec![];
    let mut current_index = 0usize;

    for ((n1, n2), _) in pair_scores {
        if !index_map.contains_key(n1) {
            index_map.insert(n1.clone(), current_index);
            current_index += 1;
            index_name.push(n1.clone());
        }
        if !index_map.contains_key(n2) {
            index_map.insert(n2.clone(), current_index);
            current_index += 1;
            index_name.push(n2.clone());
        }
    }

    // Determine the size of the matrix
    let size = index_map.len();

    // Create a new scoring matrix
    let mut matrix: crate::ScoringMatrix<f32> = crate::ScoringMatrix::new(size, 0.0, 1.0);

    // Populate the scoring matrix with the pair scores
    for ((n1, n2), score) in pair_scores {
        let i = index_map[n1];
        let j = index_map[n2];
        matrix.set(i, j, *score);
    }

    (matrix, index_name)
}
