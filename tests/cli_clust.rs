use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_similarity() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("similarity")
        .arg("tests/clust/domain.tsv")
        .arg("--mode")
        .arg("jaccard")
        .arg("--bin")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

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
        .arg("--zero")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.0669\t0.4556\t0.6260"));

    Ok(())
}

#[test]
fn command_distance_genome() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("distance")
        .arg("tests/genome/sakai.fa.gz")
        .arg("tests/genome/mg1655.fa.gz")
        .arg("-k")
        .arg("21")
        .arg("-w")
        .arg("5")
        .arg("--hasher")
        .arg("mod")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("NC_002695\tNC_000913\t0."));
    assert!(stdout.contains("NC_002128\tNC_000913\t0."));

    Ok(())
}

#[test]
fn command_distance_merge() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("distance")
        .arg("tests/clust/IBPA.fa")
        .arg("-k")
        .arg("7")
        .arg("-w")
        .arg("1")
        .arg("--merge")
        .arg("--hasher")
        .arg("murmur")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("tests/clust/IBPA.fa\ttests/clust/IBPA.fa\t763"));

    Ok(())
}

#[test]
fn command_clust_dbscan() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("dbscan")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("--eps")
        .arg("0.05")
        .arg("--min_points")
        .arg("2")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ESCF3\tA0A192CFC5_ECO25"));

    Ok(())
}

#[test]
fn command_clust_dbscan_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("dbscan")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("--eps")
        .arg("0.05")
        .arg("--min_points")
        .arg("2")
        .arg("--format")
        .arg("pair")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Each line contains a representative-member pair
    assert!(stdout.lines().count() > 0);
    assert!(
        stdout.contains("IBPA_ECOLI\tIBPA_ECOLI\n") || stdout.contains("IBPA_ESCF3\tIBPA_ESCF3\n")
    );
    assert!(
        stdout.contains("IBPA_ECOLI\tIBPA_ESCF3\n") || stdout.contains("IBPA_ESCF3\tIBPA_ECOLI\n")
    );

    Ok(())
}

#[test]
fn command_clust_cc() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("cc")
        .arg("tests/clust/IBPA.fa.05.tsv")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains("A0A192CFC5_ECO25\tIBPA_ECOLI\tIBPA_ESCF3\nIBPA_ECOLI_GA_LV"));

    Ok(())
}
