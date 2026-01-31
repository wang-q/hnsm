use assert_cmd::prelude::*;
use std::process::Command;

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
    assert!(stdout.contains("A0A192CFC5_ECO25\tIBPA_ECOLI\tIBPA_ESCF3"));

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

#[test]
fn command_clust_cc_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("cc")
        .arg("tests/clust/IBPA.fa.05.tsv")
        .arg("--format")
        .arg("pair")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains("A0A192CFC5_ECO25\tIBPA_ECOLI"));
    assert!(stdout.contains("IBPA_ECOLI_GA_LV\tIBPA_ECOLI_GA_LV"));

    Ok(())
}

#[test]
fn command_clust_mcl() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("mcl")
        .arg("tests/clust/mcl_test.tsv")
        .arg("--inflation")
        .arg("2.0")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("A\tB\tC"));
    assert!(stdout.contains("D\tE"));

    Ok(())
}

#[test]
fn command_clust_mcl_complex() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("mcl")
        .arg("tests/clust/mcl_complex.tsv")
        .arg("--inflation")
        .arg("2.0")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    // Cluster 1: n1, n2, n3, n4
    assert!(stdout.contains("n1\tn2\tn3\tn4"));
    // Cluster 2: n5, n6, n7
    assert!(stdout.contains("n5\tn6\tn7"));

    Ok(())
}

#[test]
fn command_clust_mcl_args() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("mcl")
        .arg("tests/clust/mcl_test.tsv")
        .arg("--prune")
        .arg("1e-3")
        .arg("--max_iter")
        .arg("50")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("A\tB\tC"));
    assert!(stdout.contains("D\tE"));

    Ok(())
}

#[test]
fn command_clust_mcl_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("mcl")
        .arg("tests/clust/mcl_test.tsv")
        .arg("--format")
        .arg("pair")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Cluster 1 (size 3) + Cluster 2 (size 2) = 5 pairs
    assert_eq!(stdout.lines().count(), 5);

    // Check representative output
    assert!(stdout.contains("A\tA"));
    assert!(stdout.contains("A\tB"));
    assert!(stdout.contains("A\tC"));
    assert!(stdout.contains("D\tD"));
    assert!(stdout.contains("D\tE"));

    Ok(())
}

#[test]
fn command_clust_kmedoids() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("km")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("-k")
        .arg("2")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Output should contain at least 2 lines (clusters)
    assert!(stdout.lines().count() >= 2);

    Ok(())
}

#[test]
fn command_clust_kmedoids_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("clust")
        .arg("k-medoids")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("-k")
        .arg("2")
        .arg("--format")
        .arg("pair")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Should contain tab-separated pairs
    assert!(stdout.contains("\t"));
    assert!(stdout.lines().count() >= 2);

    Ok(())
}
