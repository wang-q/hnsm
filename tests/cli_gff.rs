use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_rg_default() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("gff")
        .arg("rg")
        .arg("tests/gff_rg/test.gff")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("gene1\ttest.chr1(+):1000-2000"));
    assert!(stdout.contains("prefix:gene2\ttest.chr1(-):3000-4000"));
    assert!(!stdout.contains("mRNA1"));

    Ok(())
}

#[test]
fn command_rg_tag() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("gff")
        .arg("rg")
        .arg("tests/gff_rg/test.gff")
        .arg("--tag")
        .arg("mRNA")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("mRNA1\ttest.chr1(+):1000-2000"));
    assert!(!stdout.contains("gene1"));

    Ok(())
}

#[test]
fn command_rg_asm() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("gff")
        .arg("rg")
        .arg("tests/gff_rg/test.gff")
        .arg("--asm")
        .arg("Human")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("gene1\tHuman.chr1(+):1000-2000"));

    Ok(())
}

#[test]
fn command_rg_simplify() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("gff")
        .arg("rg")
        .arg("tests/gff_rg/test.gff")
        .arg("--simplify")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("gene2\ttest.chr1(-):3000-4000"));
    assert!(!stdout.contains("prefix:gene2"));

    Ok(())
}

#[test]
fn command_rg_case_insensitive() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("gff")
        .arg("rg")
        .arg("tests/gff_rg/test.gff")
        .arg("--tag")
        .arg("MRNA")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("mRNA1\ttest.chr1(+):1000-2000"));

    Ok(())
}
