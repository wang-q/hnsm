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
