use assert_cmd::prelude::*;
use std::process::Command;

#[test]
fn command_mat_pair() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("mat")
        .arg("pair")
        .arg("tests/clust/IBPA.mat")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 55);
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI\t0\n"));
    assert!(stdout.contains("IBPA_ECOLI\tIBPA_ECOLI_GA\t0.058"));

    Ok(())
}

#[test]
fn command_mat_phylip_full() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("mat")
        .arg("phylip")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("--mode")
        .arg("full")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 11);
    assert!(stdout.contains("IBPA_ECOLI\t0\t0.0669"));

    Ok(())
}

#[test]
fn command_mat_phylip_lower() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("mat")
        .arg("phylip")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("--mode")
        .arg("lower")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 11);
    assert!(stdout.contains("IBPA_ECOLI\n"));
    assert!(stdout.contains("IBPA_ECOLI_GA\t0.0669\n"));

    Ok(())
}

#[test]
fn command_mat_phylip_strict() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("mat")
        .arg("phylip")
        .arg("tests/clust/IBPA.fa.tsv")
        .arg("--mode")
        .arg("strict")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 11);

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[0].trim(), "10"); // 序列数量行

    // 检查第一个序列的格式
    let first_seq = lines[1];
    assert!(first_seq.starts_with("IBPA_ECOLI"));
    assert_eq!(first_seq.chars().take(10).count(), 10); // 名称长度限制
    assert!(first_seq.contains(" 0.000000")); // 格式化的距离值

    Ok(())
}
