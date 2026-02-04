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
fn command_sixframe() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd.arg("sixframe").arg("tests/fasta/trans.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 16);
    assert!(stdout.contains(">seq1(+):1-15|frame=0"));
    assert!(stdout.contains("MGMG*"));
    assert!(stdout.contains(">seq1(-):3-26|frame=2"));
    assert!(stdout.contains("TIYLYPIP"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("sixframe")
        .arg("tests/fasta/trans.fa")
        .arg("--len")
        .arg("3")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 12);

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("sixframe")
        .arg("tests/fasta/trans.fa")
        .arg("--len")
        .arg("3")
        .arg("--end")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 4);

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("sixframe")
        .arg("tests/fasta/trans.fa")
        .arg("--len")
        .arg("3")
        .arg("--start")
        .arg("--end")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);

    Ok(())
}
