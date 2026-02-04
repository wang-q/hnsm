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
fn command_mask() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("mask")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/mask.json")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("read0\ntcgtttaacccaaatcaagg"), "read0");
    assert!(stdout.contains("read2\natagcaagct"), "read2");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("mask")
        .arg("--hard")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/mask.json")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(stdout.contains("read0\nNNNNNNNNNNNNNNNNNNNN"), "read0");
    assert!(stdout.contains("read2\nNNNNNNNNNN"), "read2");

    Ok(())
}




// faops filter -l 0 -a 10 -z 50 tests/fasta/ufasta.fa stdout
// faops filter -l 0 -a 1 -u <(cat tests/fasta/ufasta.fa tests/fasta/ufasta.fa) stdout
#[test]
fn command_filter() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/ufasta.fa")
        .arg("-a")
        .arg("10")
        .arg("-z")
        .arg("50")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 12);
    assert!(!stdout.contains(">read0"), "read0");
    assert!(stdout.contains(">read20"), "read20");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/ufasta.fa.gz")
        .arg("--uniq")
        .arg("-a")
        .arg("1")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 90);

    Ok(())
}

#[test]
fn command_filter_fmt() -> anyhow::Result<()> {
    // faops filter -N tests/fasta/filter.fa stdout
    // faops treats '-' as N, which is incorrect
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--iupac")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">iupac\nAMRG"), "iupac");
    assert!(stdout.contains(">iupac\nANNG"), "iupac");
    assert!(stdout.contains(">dash\nA-NG"), "dash not changed");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--dash")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">dash\nA-RG"), "dash");
    assert!(stdout.contains(">dash\nARG"), "dash");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--upper")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">upper\nAtcG"), "upper");
    assert!(stdout.contains(">upper\nATCG"), "upper");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--simplify")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert!(!stdout.contains(">read.1 simplify\nAGGG"), "simplify");
    assert!(stdout.contains(">read\nAGGG"), "simplify");

    Ok(())
}

#[test]
fn command_dedup() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd.arg("dedup").arg("tests/fasta/dedup.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 8);
    assert!(!stdout.contains(">read0 some text"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--desc")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">read0 some text"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 6);
    assert!(!stdout.contains(">read1"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--case")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 4);
    assert!(!stdout.contains(">read2"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--both")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);
    assert!(!stdout.contains(">read3"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--both")
        .arg("--file")
        .arg("stdout")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains(">read0"));
    assert!(stdout.contains("read0\tread3"));

    Ok(())
}



#[test]
fn command_sixframe() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd.arg("sixframe").arg("tests/fasta/trans.fa").output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 16);
    assert!(stdout.contains(">seq1(+):1-15|frame=0"));
    assert!(stdout.contains("MGMG*"));
    assert!(stdout.contains(">seq1(-):3-26|frame=2"));
    assert!(stdout.contains("TIYLYPIP"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("sixframe")
        .arg("tests/fasta/trans.fa")
        .arg("--len")
        .arg("3")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 12);

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("sixframe")
        .arg("tests/fasta/trans.fa")
        .arg("--len")
        .arg("3")
        .arg("--end")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 4);

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("sixframe")
        .arg("tests/fasta/trans.fa")
        .arg("--len")
        .arg("3")
        .arg("--start")
        .arg("--end")
        .output()?;
    let stdout = String::from_utf8(output.stdout)?;

    assert_eq!(stdout.lines().count(), 2);

    Ok(())
}
