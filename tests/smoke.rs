//! Integration smoke test: `--synthetic` runs end-to-end and writes a
//! well-formed JSON document.

use std::process::Command;

#[test]
fn synthetic_demo_runs_end_to_end() {
    let tmp = tempdir();
    let out = tmp.join("tail_synthetic.json");

    let exe = env!("CARGO_BIN_EXE_evt-cli");
    let status = Command::new(exe)
        .args([
            "--synthetic",
            "--synthetic-n", "2000",
            "--output", out.to_str().unwrap(),
        ])
        .status()
        .expect("spawn evt-cli");
    assert!(status.success(), "evt-cli exited with {status:?}");

    let txt = std::fs::read_to_string(&out).expect("read output");
    let v: serde_json::Value = serde_json::from_str(&txt).expect("parse JSON");
    assert_eq!(v["source"], "synthetic");
    assert_eq!(v["labels"].as_array().unwrap().len(), 4);
    assert_eq!(v["chi_matrix"].as_array().unwrap().len(), 4);
    let fits = v["fits"].as_array().unwrap();
    assert_eq!(fits.len(), 4);
    for f in fits {
        let xi = f["xi"].as_f64().unwrap();
        let beta = f["beta"].as_f64().unwrap();
        assert!(xi.is_finite() && beta > 0.0);
    }
}

fn tempdir() -> std::path::PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("tail-evt-test-{}", std::process::id()));
    std::fs::create_dir_all(&p).unwrap();
    p
}
