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

#[test]
fn command_range_update() -> anyhow::Result<()> {
    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    // Copy test file to temp directory
    std::fs::copy(
        "tests/index/final.contigs.fa",
        format!("{}/test.fa", tempdir_str),
    )?;

    // First run, create index
    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("range")
        .arg(format!("{}/test.fa", tempdir_str))
        .arg("k81_130:11-20")
        .output()?;

    // Get index file's modification time
    let loc_file = format!("{}/test.fa.loc", tempdir_str);
    let first_modified = std::fs::metadata(&loc_file)?.modified()?;

    // Wait for a second to ensure different timestamps
    std::thread::sleep(std::time::Duration::from_secs(1));

    // Force update index with --update
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("range")
        .arg(format!("{}/test.fa", tempdir_str))
        .arg("k81_130:11-20")
        .arg("--update")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Verify output content
    assert!(stdout.contains(">k81_130:11-20\nGGTGAATCAA\n"));

    // Verify index file was updated
    let second_modified = std::fs::metadata(&loc_file)?.modified()?;
    assert!(second_modified > first_modified);

    tempdir.close()?;
    Ok(())
}
