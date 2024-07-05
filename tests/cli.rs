use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::TempDir;

#[test]
fn command_invalid() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("foobar");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("recognized"));

    Ok(())
}

#[test]
fn file_doesnt_provided() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("size");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("not provided"));

    Ok(())
}

#[test]
fn file_doesnt_exist() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("size").arg("tests/file/doesnt/exist");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("could not open"));

    Ok(())
}

#[test]
fn command_size() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("size")
        .arg("tests/fasta/ufasta.fa")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 50);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    let mut sum = 0;
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            sum += fields[1].parse::<i32>().unwrap();
        }
    }
    assert_eq!(sum, 9317, "sum length");

    Ok(())
}

#[test]
fn command_size_gz() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("size")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/ufasta.fa.gz")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    Ok(())
}
