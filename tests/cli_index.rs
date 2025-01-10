use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::TempDir;

#[test]
fn command_gz() -> anyhow::Result<()> {
    match which::which("bgzip") {
        Err(_) => return Ok(()),
        Ok(_) => {}
    }

    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("gz")
        .arg("tests/index/final.contigs.fa")
        .arg("-o")
        .arg(format!("{}/ctg.fa", tempdir_str))
        .assert()
        .success()
        .stdout(predicate::str::is_empty());

    assert!(&tempdir.path().join("ctg.fa.gz").exists());
    assert!(&tempdir.path().join("ctg.fa.gz.gzi").exists());

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_range() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("range")
        .arg("tests/index/final.contigs.fa.gz")
        .arg("k81_130")
        .arg("k81_130:11-20")
        .arg("k81_170:304-323")
        .arg("k81_170(-):1-20")
        .arg("k81_158:70001-70020")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">k81_130\nAGTTTCAACT"));
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));
    assert!(stdout.contains(">k81_170:304-323\nAGTTAAAAACCTGATTTATT\n"));
    assert!(stdout.contains(">k81_170(-):1-20\nATTAACCTGTTGTAGGTGTT\n"));
    assert!(stdout.contains(">k81_158:70001-70020\nTGGCTATAACCTAATTTTGT\n"));

    Ok(())
}

#[test]
fn command_range_r() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("range")
        .arg("tests/index/final.contigs.fa.gz")
        .arg("-r")
        .arg("tests/index/sample.rg")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 12);
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));

    Ok(())
}
