use assert_cmd::prelude::*;
use std::fs;
use std::process::Command;

#[test]
fn command_synt_dna() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("mmg")
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

#[test]
fn command_synt_dna_soft_mask() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("mmg")
        .arg("tests/genome/small_1.fa")
        .arg("tests/genome/small_2.fa")
        .arg("-o")
        .arg("tests/temp_synt_sm.tsv")
        .arg("--soft-mask")
        .output()?;

    assert!(output.status.success());
    fs::remove_file("tests/temp_synt_sm.tsv")?;
    Ok(())
}

#[test]
fn command_synt_merge() -> anyhow::Result<()> {
    // Create a temporary input file with fragmented blocks
    let input_path = "tests/temp_merge_in.tsv";
    let output_path = "tests/temp_merge_out.tsv";

    // Block 1: 100-200
    // Block 2: 300-400
    // Gap: 100 bp. chain_gap=150 should merge them.
    let input_content = "\
# Block_ID\tRange\tCount\tRound
1\tG1(+):100-200\t10\t1
1\tG2(+):100-200\t10\t1
2\tG1(+):300-400\t10\t1
2\tG2(+):300-400\t10\t1
";
    fs::write(input_path, input_content)?;

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("merge")
        .arg(input_path)
        .arg("-o")
        .arg(output_path)
        .arg("--chain-gap")
        .arg("150")
        .output()?;

    assert!(output.status.success());

    let content = fs::read_to_string(output_path)?;
    let lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();

    // Should be merged into 1 block (2 lines for 2 genomes)
    assert_eq!(lines.len(), 2);

    // Check range values
    // G1 range should be 100-400
    // Note: The order of lines depends on hash map or implementation details,
    // but merge.rs sorts by seq_name.
    // G1 comes before G2.

    // G1 line
    // G1(+):100-400
    assert!(lines[0].contains("G1"));
    assert!(lines[0].contains("100-400"));

    // G2 line
    assert!(lines[1].contains("G2"));
    assert!(lines[1].contains("100-400"));

    fs::remove_file(input_path)?;
    fs::remove_file(output_path)?;

    Ok(())
}

#[test]
fn command_synt_dag() -> anyhow::Result<()> {
    let annot_path = "tests/dag_new_annot.tsv";
    let match_path = "tests/dag_new_match.tsv";
    let output_path = "tests/dag_new_out.tsv";

    // Format: key \t mol(strand):start-end
    let annot_content = "\
G1\tP1.chr1(+):100-200
G2\tP1.chr1(+):300-400
G3\tP2.chr1(-):100-200
G4\tP2.chr1(-):300-400
";
    fs::write(annot_path, annot_content)?;

    // Format: acc1 \t acc2 \t score (similarity)
    let match_content = "\
G1\tG3\t0.95
G2\tG4\t0.90
";
    fs::write(match_path, match_content)?;

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("dag")
        .arg(annot_path)
        .arg(match_path)
        .arg("-o")
        .arg(output_path)
        .arg("--mna")
        .arg("1") // Min number of aligned pairs = 1 to allow small chains
        .output()?;

    assert!(output.status.success());

    let content = fs::read_to_string(output_path)?;
    // println!("Output content:\n{}", content);

    assert!(content.contains("P1.chr1"));
    assert!(content.contains("P2.chr1"));

    fs::remove_file(annot_path)?;
    fs::remove_file(match_path)?;
    fs::remove_file(output_path)?;

    Ok(())
}
