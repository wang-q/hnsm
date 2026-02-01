use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_fasta as fasta;
use std::io::{BufRead, Read, Seek, SeekFrom};

pub enum Input {
    Buf(Box<dyn BufRead>),
    File(std::fs::File),
    Bgzf(bgzf::io::IndexedReader<std::fs::File>),
}

pub fn create_loc(infile: &str, locfile: &str, is_bgzf: bool) -> anyhow::Result<()> {
    let mut reader = if is_bgzf {
        // http://www.htslib.org/doc/bgzip.html
        // Bgzip will attempt to ensure BGZF blocks end on a newline when the input is a text file.
        // The exception to this is where a single line is larger than a BGZF block (64Kb).
        Input::Bgzf(bgzf::io::indexed_reader::Builder::default().build_from_path(infile)?)
    } else {
        Input::Buf(reader_buf(infile))
    };

    let mut writer: Box<dyn std::io::Write> =
        Box::new(std::io::BufWriter::new(std::fs::File::create(locfile)?));

    // https://www.ginkgobioworks.com/2023/03/17/even-more-rapid-retrieval-from-very-large-files-with-rust/
    let mut record_size = 0; // including header, sequence, newlines
    let mut offset = 0;
    let mut line = String::new();
    while let Ok(num) = match &mut reader {
        Input::Buf(rdr) => rdr.read_line(&mut line),
        Input::Bgzf(rdr) => rdr.read_line(&mut line),
        &mut Input::File(_) => unreachable!(),
    } {
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

pub fn reader_buf(infile: &str) -> Box<dyn BufRead> {
    let reader: Box<dyn BufRead> = {
        let path = std::path::Path::new(infile);
        let file = match std::fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        Box::new(std::io::BufReader::new(file))
    };

    reader
}

pub fn load_loc(loc_file: &str) -> anyhow::Result<IndexMap<String, (u64, usize)>> {
    let mut reader = reader_buf(loc_file);

    let mut loc_of: IndexMap<String, (u64, usize)> = IndexMap::new();
    let mut line = String::new();
    while let Ok(num) = reader.by_ref().read_line(&mut line) {
        if num == 0 {
            break;
        }
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() != 3 {
            continue;
        }

        loc_of.insert(
            fields[0].to_string(),
            (fields[1].parse::<u64>()?, fields[2].parse::<usize>()?),
        );

        line.clear();
    }

    Ok(loc_of)
}

pub fn record_rg(
    reader: &mut Input,
    loc_of: &IndexMap<String, (u64, usize)>,
    rg: &str,
) -> anyhow::Result<fasta::Record> {
    let (offset, size) = loc_of.get(rg).unwrap();

    let data_buf = read_offset(reader, *offset, *size)?;
    let mut fa_in = fasta::io::Reader::new(&data_buf[..]);

    fa_in.read_definition(&mut String::new())?;
    let mut buf = Vec::new();
    fa_in.read_sequence(&mut buf)?;

    let definition = fasta::record::Definition::new(rg, None);
    let sequence = fasta::record::Sequence::from(buf);
    let record = fasta::Record::new(definition, sequence);

    Ok(record)
}

pub fn records_offset(
    reader: &mut Input,
    offset: u64,
    size: usize,
) -> anyhow::Result<Vec<fasta::Record>> {
    let mut records = Vec::new();

    let data_buf = read_offset(reader, offset, size)?;
    let mut fa_in = fasta::io::Reader::new(&data_buf[..]);

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;
        records.push(record);
    }

    Ok(records)
}

pub fn read_offset(reader: &mut Input, offset: u64, size: usize) -> anyhow::Result<Vec<u8>> {
    let mut data_buf = vec![0; size];

    match reader {
        Input::File(rdr) => {
            rdr.seek(SeekFrom::Start(offset))?;
            rdr.read_exact(&mut data_buf)?;
        }
        Input::Bgzf(rdr) => {
            rdr.seek(SeekFrom::Start(offset))?;
            rdr.read_exact(&mut data_buf)?;
        }
        Input::Buf(_) => unreachable!(),
    }

    Ok(data_buf)
}
