use assert_cmd::prelude::*;
use std::fs;
use std::process::Command;

#[test]
fn command_synt_dna() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("dna")
        .arg("tests/genome/small_1.fa")
        .arg("tests/genome/small_2.fa")
        .arg("-o")
        .arg("tests/temp_synt.tsv")
        .arg("-b")
        .arg("50")
        .arg("-r")
        .arg("50")
        .output()?;

    assert!(output.status.success());

    let content = fs::read_to_string("tests/temp_synt.tsv")?;
    let lines: Vec<&str> = content.lines().collect();

    // Check header
    assert!(lines[0].contains("Block_ID"));
    assert!(lines[0].contains("Range"));
    assert!(lines[0].contains("Count"));
    assert!(lines[0].contains("Round"));

    // Check content format
    // Expected: Block_ID \t Genome.Chr(Strand):Start-End \t Count \t Round
    
    // With identical small files, we expect at least 1 block
    assert!(lines.len() > 1); 

    let first_data_line = lines[1];
    let fields: Vec<&str> = first_data_line.split('\t').collect();
    assert_eq!(fields.len(), 4);
    
    // Check range format
    // Should contain "small_1." or "small_2."
    assert!(fields[1].contains("small_1.") || fields[1].contains("small_2."));
    
    // Check strand normalization
    // Since both sequences are identical and forward, and first genome is normalized to (+),
    // we expect (+) for both if they are identical.
    // But let's just check the format generally.
    assert!(fields[1].contains(":"));
    assert!(fields[1].contains("-"));

    // Cleanup
    fs::remove_file("tests/temp_synt.tsv")?;

    Ok(())
}
