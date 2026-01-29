use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_interleave() -> anyhow::Result<()> {
    // count empty seqs
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa.gz")
        .arg("tests/fasta/ufasta.fa")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "").count(), 10);

    // count empty seqs (single)
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "").count(), 5);

    // count empty seqs (single)
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("interleave")
        .arg("tests/fasta/ufasta.fa")
        .arg("--fq")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "").count(), 10);

    Ok(())
}

#[test]
fn command_interleave_fq() -> anyhow::Result<()> {
    // fq
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("interleave")
        .arg("--fq")
        .arg("tests/fastq/R1.fq.gz")
        .arg("tests/fastq/R2.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "!").count(), 0);
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "+").count(), 50);
    assert_eq!(
        stdout
            .lines()
            .into_iter()
            .filter(|e| e.ends_with("/1"))
            .count(),
        25
    );
    assert_eq!(
        stdout
            .lines()
            .into_iter()
            .filter(|e| e.ends_with("/2"))
            .count(),
        25
    );

    // fq (single)
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("interleave")
        .arg("--fq")
        .arg("tests/fastq/R1.fq.gz")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "!").count(), 25);
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "+").count(), 50);
    assert_eq!(
        stdout
            .lines()
            .into_iter()
            .filter(|e| e.ends_with("/1"))
            .count(),
        25
    );
    assert_eq!(
        stdout
            .lines()
            .into_iter()
            .filter(|e| e.ends_with("/2"))
            .count(),
        25
    );

    Ok(())
}

#[test]
fn command_fq2fa() -> anyhow::Result<()> {
    // Basic conversion test
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd.arg("fq2fa").arg("tests/fastq/R1.fq.gz").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Verify output format
    assert_eq!(
        stdout
            .lines()
            .into_iter()
            .filter(|e| e.starts_with(">"))
            .count(),
        25
    );
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "").count(), 0);
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "+").count(), 0);
    assert_eq!(stdout.lines().into_iter().filter(|e| *e == "!").count(), 0);

    // Test file output
    let mut cmd = Command::cargo_bin("hnsm")?;
    let temp_file = "tests/fastq/temp.fa";
    cmd.arg("fq2fa")
        .arg("tests/fastq/R1.fq.gz")
        .arg("-o")
        .arg(temp_file)
        .assert()
        .success();

    // Read and verify output file
    let output = std::fs::read_to_string(temp_file)?;
    assert_eq!(
        output
            .lines()
            .into_iter()
            .filter(|e| e.starts_with(">"))
            .count(),
        25
    );

    // Cleanup temporary file
    std::fs::remove_file(temp_file)?;

    Ok(())
}
