use std::collections::BTreeMap;
use std::io::{self, BufRead, Read, Write};
use std::{ffi, fs, path};

// https://www.ginkgobioworks.com/2023/03/17/even-more-rapid-retrieval-from-very-large-files-with-rust/
pub fn create_loc(infile: &str, locfile: &str) -> anyhow::Result<()> {
    let mut reader = reader_text(infile);

    let mut writer: Box<dyn io::Write> =
        Box::new(io::BufWriter::new(fs::File::create(locfile).unwrap()));

    let mut record_size = 0; // including header, sequence, newlines
    let mut offset = 0;
    let mut line = String::new();
    while let Ok(num) = reader.by_ref().read_line(&mut line) {
        if num == 0 {
            break;
        }

        if line.starts_with('>') {
            if record_size > 0 {
                // the size of the previous record
                writer.write_fmt(format_args!("\t{}\n", record_size))?;
            }
            // reset size counter for new record
            record_size = 0;

            //current record name
            let name = &line[1..]
                .splitn(2, |c: char| c.is_ascii_whitespace())
                .next()
                .unwrap();
            writer.write_fmt(format_args!("{}\t{}", name, offset))?;
        }

        record_size += num;
        offset += num;
        line.clear();
    }
    if record_size > 0 {
        writer.write_fmt(format_args!("\t{}\n", record_size))?;
    }

    Ok(())
}

pub fn reader_text(infile: &str) -> Box<dyn BufRead> {
    let mut reader: Box<dyn io::BufRead> = {
        let path = path::Path::new(infile);
        let file = match fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        Box::new(io::BufReader::new(file))
    };

    reader
}

pub fn load_loc(locfile: &String) -> anyhow::Result<BTreeMap<String, (u64, usize)>> {
    let mut reader = reader_text(locfile);

    let mut loc_of: BTreeMap<String, (u64, usize)> = BTreeMap::new();
    let mut line = String::new();
    while let Ok(num) = reader.by_ref().read_line(&mut line) {
        if num == 0 {
            break;
        }
        let mut fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 3 {
            continue;
        }

        loc_of.insert(
            fields[0].to_string(),
            (
                fields[1].parse::<u64>().unwrap(),
                fields[2].parse::<usize>().unwrap(),
            ),
        );

        line.clear();
    }

    Ok(loc_of)
}
