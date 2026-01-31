use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_mat_phylip() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("nwr")?;
    let output = cmd
        .arg("mat")
        .arg("phylip")
        .arg("tests/mat/IBPA.fa.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 11);
    assert!(stdout.contains("IBPA_ECOLI\t0\t0.0669"));

    Ok(())
}
