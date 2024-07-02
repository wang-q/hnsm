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
    cmd.arg("sixframe");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("not provided"));

    Ok(())
}

#[test]
fn file_doesnt_exist() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("sixframe").arg("tests/file/doesnt/exist");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("could not open"));

    Ok(())
}
