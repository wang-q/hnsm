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
    cmd.arg("size");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("not provided"));

    Ok(())
}

#[test]
fn file_doesnt_exist() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("size").arg("tests/file/doesnt/exist");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("could not open"));

    Ok(())
}

#[test]
fn command_size() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("size")
        .arg("tests/fasta/ufasta.fa")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 50);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    let mut sum = 0;
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            sum += fields[1].parse::<i32>().unwrap();
        }
    }
    assert_eq!(sum, 9317, "sum length");

    Ok(())
}

#[test]
fn command_size_gz() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("size")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/ufasta.fa.gz")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 100);
    assert!(stdout.contains("read0\t359"), "read0");
    assert!(stdout.contains("read1\t106"), "read1");

    Ok(())
}

#[test]
fn command_some() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("some")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/list.txt")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("read0\n"), "read0");
    assert!(stdout.contains("read12\n"), "read12");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("some")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/list.txt")
        .arg("-i")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 91);
    assert!(!stdout.contains("read0\n"), "read0");
    assert!(!stdout.contains("read12\n"), "read12");

    Ok(())
}

#[test]
fn command_order() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("order")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/list.txt")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 4);
    assert!(stdout.contains("read12\n"), "read12");
    assert!(stdout.contains("read0\n"), "read0");

    Ok(())
}

#[test]
fn command_one() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("one")
        .arg("tests/fasta/ufasta.fa")
        .arg("read12")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("read12\n"), "read12");

    Ok(())
}

#[test]
fn command_masked() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("masked")
        .arg("tests/fasta/ufasta.fa")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.contains("read46:3-4"), "read46");

    Ok(())
}

#[test]
fn command_rc() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd.arg("rc").arg("tests/fasta/ufasta.fa").output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.contains("GgacTgcggCTagAA"), "read46");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("rc")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/list.txt")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.contains(">RC_read12"), "read12");
    assert!(!stdout.contains(">RC_read46"), "read46");
    assert!(!stdout.contains("GgacTgcggCTagAA"), "read46");

    Ok(())
}

#[test]
fn command_count() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("count")
        .arg("tests/fasta/ufasta.fa")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.contains("read45\t0\t0"), "empty");
    assert!(stdout.contains("total\t9317\t2318"), "total");

    Ok(())
}

#[test]
fn command_replace() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("replace")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/replace.tsv")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 95);
    assert!(stdout.contains(">359"), "read0");
    assert!(!stdout.contains(">read0"), "read0");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("replace")
        .arg("tests/fasta/ufasta.fa")
        .arg("tests/fasta/replace.tsv")
        .arg("--some")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 6);
    assert!(stdout.contains(">359"), "read0");
    assert!(!stdout.contains(">read0"), "read0");

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
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

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
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

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
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(!stdout.contains(">iupac\nAMRG"), "iupac");
    assert!(stdout.contains(">iupac\nANNG"), "iupac");
    assert!(stdout.contains(">dash\nA-NG"), "dash not changed");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--dash")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(!stdout.contains(">dash\nA-RG"), "dash");
    assert!(stdout.contains(">dash\nARG"), "dash");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--upper")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(!stdout.contains(">upper\nAtcG"), "upper");
    assert!(stdout.contains(">upper\nATCG"), "upper");

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("filter")
        .arg("tests/fasta/filter.fa")
        .arg("--simplify")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(!stdout.contains(">read.1 simplify\nAGGG"), "simplify");
    assert!(stdout.contains(">read\nAGGG"), "simplify");

    Ok(())
}

#[test]
fn command_dedup() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 8);
    assert!(!stdout.contains(">read0 some text"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--desc")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 10);
    assert!(stdout.contains(">read0 some text"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 6);
    assert!(!stdout.contains(">read1"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--case")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 4);
    assert!(!stdout.contains(">read2"));

    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("dedup")
        .arg("tests/fasta/dedup.fa")
        .arg("--seq")
        .arg("--both")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

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
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 7);
    assert!(stdout.contains(">read0"));
    assert!(stdout.contains("read0\tread3"));

    Ok(())
}

#[test]
fn command_split_name() -> anyhow::Result<()> {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("split")
        .arg("name")
        .arg("tests/fasta/ufasta.fa")
        .arg("-o")
        .arg(tempdir_str)
        .assert()
        .success()
        .stdout(predicate::str::is_empty());

    assert!(&tempdir.path().join("read0.fa").is_file());
    assert!(!&tempdir.path().join("000.fa").exists());

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_split_about() -> anyhow::Result<()> {
    let tempdir = TempDir::new().unwrap();
    let tempdir_str = tempdir.path().to_str().unwrap();

    let mut cmd = Command::cargo_bin("hnsm")?;
    cmd.arg("split")
        .arg("about")
        .arg("tests/fasta/ufasta.fa")
        .arg("-c")
        .arg("2000")
        .arg("-o")
        .arg(tempdir_str)
        .assert()
        .success()
        .stdout(predicate::str::is_empty());

    assert!(!&tempdir.path().join("read0.fa").is_file());
    assert!(&tempdir.path().join("000.fa").exists());
    assert!(&tempdir.path().join("004.fa").exists());
    assert!(!&tempdir.path().join("005.fa").exists());

    tempdir.close()?;
    Ok(())
}

#[test]
fn command_n50() -> anyhow::Result<()> {
    // display header
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("N50\t314"), "line 1");

    // doesn't display header
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 1);
    assert!(!stdout.contains("N50\t314"), "line 1");
    assert!(stdout.contains("314"), "line 1");

    // set genome size (NG50)
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .arg("-g")
        .arg("10000")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 1);
    assert!(stdout.contains("297"), "line 1");

    // sum and average of size
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .arg("-S")
        .arg("-A")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.contains("314\n9317\n186.34"), "line 1,2,3");

    // N10, N90, E-size
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-H")
        .arg("-E")
        .arg("-N")
        .arg("10")
        .arg("-N")
        .arg("90")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 3);
    assert!(stdout.contains("516\n112\n314.70\n"), "line 1,2,3");

    // transposed
    let mut cmd = Command::cargo_bin("hnsm")?;
    let output = cmd
        .arg("n50")
        .arg("tests/fasta/ufasta.fa")
        .arg("-E")
        .arg("-N")
        .arg("10")
        .arg("-N")
        .arg("90")
        .arg("-t")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("N10\tN90\tE\n"), "line 1");
    assert!(stdout.contains("516\t112\t314.70\n"), "line 2");

    Ok(())
}
