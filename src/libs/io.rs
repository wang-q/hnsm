use std::collections::HashMap;
use std::io::{BufRead, Write};

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

pub fn load_pair_scores(infile: &str) -> (Vec<((usize, usize), f32)>, Vec<String>) {
    let mut pair_scores = Vec::new();
    // Create a mapping from string identifiers to indices
    let mut index_map = HashMap::new();
    let mut index_name = vec![];
    let mut current_index = 0usize;

    let reader = intspan::reader(infile);
    for line in reader.lines().map_while(Result::ok) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let n1 = fields[0].to_string();
            let n2 = fields[1].to_string();
            let score: f32 = fields[2].parse::<f32>().unwrap();

            if !index_map.contains_key(&n1) {
                index_map.insert(n1.clone(), current_index);
                current_index += 1;
                index_name.push(n1.clone());
            }
            if !index_map.contains_key(&n2) {
                index_map.insert(n2.clone(), current_index);
                current_index += 1;
                index_name.push(n2.clone());
            }

            pair_scores.push(((index_map[&n1], index_map[&n2]), score));
        }
    }
    (pair_scores, index_name)
}

pub fn populate_matrix(
    pair_scores: &Vec<((usize, usize), f32)>,
    index_name: &Vec<String>,
    same: f32,
    missing: f32,
) -> crate::ScoringMatrix<f32> {
    let size = index_name.len();
    let mut matrix: crate::ScoringMatrix<f32> = crate::ScoringMatrix::new(size, same, missing);
    for ((i, j), score) in pair_scores {
        matrix.set(*i, *j, *score);
    }
    matrix
}
