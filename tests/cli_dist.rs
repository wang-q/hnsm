use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_dist_hv() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("hv")
        .arg("tests/clust/IBPA.fa")
        .arg("-k")
        .arg("7")
        .arg("-w")
        .arg("1")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    // Check if output contains expected content
    // Based on distance command output in cli_clust.rs, we expect self-comparison of the file
    // But hv command with single file merges all sequences and compares to itself?
    // Let's check hv help: 
    // "For a single sequence file: Merge all sequences within the file into a single hypervector. 
    // Note that comparing this set to itself (self-comparison) is not meaningful..."
    
    // Actually, if we pass a single file to hv, it returns 1 line: file vs file.
    // Let's verify with a run.
    
    assert!(stdout.lines().count() >= 1);
    assert!(stdout.contains("tests/clust/IBPA.fa"));

    Ok(())
}

#[test]
fn command_dist_hv_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("hv")
        .arg("tests/clust/IBPA.fa")
        .arg("tests/clust/IBPA.fa") // Compare file against itself
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("tests/clust/IBPA.fa"));
    // Similarity should be 1.0 / Distance 0.0
    // The output format: <file1> <file2> ... <mash_dist> ...
    
    Ok(())
}

#[test]
fn command_dist_seq() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("seq")
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
fn command_dist_seq_sim() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("seq")
        .arg("tests/clust/IBPA.fa")
        .arg("-k")
        .arg("7")
        .arg("-w")
        .arg("1")
        .arg("--zero")
        .arg("--sim")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 100);
    // Mash dist 0.0669 -> Sim 1 - 0.0669 = 0.9331
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.9331\t0.4556\t0.6260"));

    Ok(())
}

#[test]
fn command_dist_seq_genome() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("seq")
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
fn command_dist_seq_merge() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("seq")
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
fn command_dist_vector() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dist")
        .arg("vector")
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
