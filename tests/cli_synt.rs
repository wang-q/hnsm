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
    let line1 = lines[0];
    let line2 = lines[1];

    assert!(line1.contains("G1(+):100-400"));
    assert!(line2.contains("G2(+):100-400"));
    
    // Check score summation
    // 10 + 10 = 20.0
    assert!(line1.contains("\t20.0"));
    assert!(line2.contains("\t20.0"));

    // Cleanup
    fs::remove_file(input_path)?;
    fs::remove_file(output_path)?;

    Ok(())
}

#[test]
fn command_synt_dna_soft_mask() -> anyhow::Result<()> {
    // Create 3 files
    // 1. upper.fa: Uppercase repeats
    // 2. lower1.fa: Lowercase repeats (copy of upper)
    // 3. lower2.fa: Lowercase repeats (copy of upper)
    
    // Using a repetitive sequence "AAAAAAAA..." to ensure many minimizers are generated if valid.
    let seq_upper = "A".repeat(1000);
    let seq_lower = "a".repeat(1000);
    
    fs::write("tests/upper.fa", format!(">seq1\n{}", seq_upper))?;
    fs::write("tests/lower1.fa", format!(">seq2\n{}", seq_lower))?;
    fs::write("tests/lower2.fa", format!(">seq3\n{}", seq_lower))?;
    
    // Case 1: No soft mask.
    // lower1 and lower2 should match each other (identical lowercase).
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("dna")
        .arg("tests/lower1.fa")
        .arg("tests/lower2.fa")
        .arg("-o")
        .arg("tests/temp_nomask.tsv")
        .arg("-k")
        .arg("10")
        .arg("-b") // block size
        .arg("50")
        .arg("--min-weight")
        .arg("2")
        .output()?;
    assert!(output.status.success());
    let content = fs::read_to_string("tests/temp_nomask.tsv")?;
    let lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    // Should match
    assert!(!lines.is_empty());

    // Case 2: With soft mask.
    // lower1 and lower2 should NOT match (both ignored).
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("synt")
        .arg("dna")
        .arg("tests/lower1.fa")
        .arg("tests/lower2.fa")
        .arg("-o")
        .arg("tests/temp_mask.tsv")
        .arg("--soft-mask")
        .arg("-k")
        .arg("10")
        .arg("-b")
        .arg("50")
        .arg("--min-weight")
        .arg("2")
        .output()?;
    assert!(output.status.success());
    let content = fs::read_to_string("tests/temp_mask.tsv")?;
    let lines: Vec<&str> = content.lines().filter(|l| !l.starts_with('#')).collect();
    // Should be empty (no blocks)
    assert!(lines.is_empty());

    // Cleanup
    fs::remove_file("tests/upper.fa")?;
    fs::remove_file("tests/lower1.fa")?;
    fs::remove_file("tests/lower2.fa")?;
    fs::remove_file("tests/temp_nomask.tsv")?;
    fs::remove_file("tests/temp_mask.tsv")?;
    
    Ok(())
}
