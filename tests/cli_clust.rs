use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::TempDir;

#[test]
fn command_similarity() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("similarity")
        .arg("tests/clust/domain.tsv")
        .arg("--mode")
        .arg("jaccard")
        .arg("--bin")
        .arg("--parallel")
        .arg("1")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout
        .contains("Acin_baum_1326584_GCF_025854095_1\tAcin_baum_1326584_GCF_025854095_1\t1.0000"));

    Ok(())
}

#[test]
fn command_distance() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("distance")
        .arg("tests/clust/IBPA.fa")
        .arg("-k")
        .arg("7")
        .arg("-w")
        .arg("1")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.0669\t0.4556\t0.6260"));

    Ok(())
}

#[test]
fn command_cluster() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("cluster")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("--mode")
        .arg("matrix")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("IBPA_ECOLI\t0\t0.0669"));

    Ok(())
}
