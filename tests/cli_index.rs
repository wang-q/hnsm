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

    let tempdir = TempDir::new().unwrap();
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
