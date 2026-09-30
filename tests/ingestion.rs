#[path = "../src/ingestion.rs"]
mod ingestion;
use ingestion::*;
fn asset(s: &str) -> Asset {
    Asset::from_bytes(s.as_bytes().to_vec(), usize::MAX).unwrap()
}
#[test]
fn released_rows_are_preserved_not_selected_or_deduplicated() {
    let a = asset(
        "CID IDSURVEY zHD zCMB zHEL m_b_corr EXTRA\nA 51 .009 .011 .012 17 keep\nA 52 .020 .022 .023 18 also\n",
    );
    let hash = a.sha256().to_owned();
    let t = pantheon_plus(a, 2).unwrap();
    assert_eq!(t.event_ids, vec!["A", "A"]);
    assert_ne!(t.measurement_ids[0], t.measurement_ids[1]);
    assert_eq!(t.measurement_ids[0], format!("{hash}:row:0"));
    assert_eq!(t.zhd, vec![0.009, 0.020]);
    assert_eq!(t.zcmb, vec![0.011, 0.022]);
    assert_eq!(t.zhel, vec![0.012, 0.023]);
    assert_eq!(t.original_fields[1][6], "also");
    assert_eq!(t.values, vec![17., 18.]);
    assert_eq!(t.quality, vec![0, 0]);
    assert_eq!(t.original_columns[0], "CID");
    assert_eq!(t.asset.sha256(), hash);
    let c =
        pantheon_covariance(asset("2\n4 1.00000003 1 9\n"), t.measurement_ids.clone(), 4).unwrap();
    assert_eq!(c.dimension, 2);
    assert!(c.ordering_provenance.contains("not independently verified"));
    assert_eq!(c.values[1], 1.00000003);
    assert_eq!(c.values[2], 1.);
    assert_eq!(c.axis_ids, t.measurement_ids);
    assert!(!c.asset.bytes().is_empty());
}
#[test]
fn strict_structural_failures() {
    for s in [
        "CID CID\nA A",
        "CID IDSURVEY zHD zCMB zHEL m_b_corr\nA 1 2 3 4",
        "CID IDSURVEY zHD zCMB zHEL m_b_corr\nA 1 2 3 4 x",
    ] {
        assert!(pantheon_plus(asset(s), 10).is_err());
    }
    for s in ["2 1 2 3", "2 1 2 3 4 5", "18446744073709551615"] {
        assert!(pantheon_covariance(asset(s), vec!["a".into(), "b".into()], 4).is_err());
    }
}

#[test]
fn bounded_read_hash_changes_when_path_bytes_change() {
    let path = std::env::temp_dir().join(format!("irred-ingestion-{}", std::process::id()));
    std::fs::write(&path, b"abc").unwrap();
    assert!(read_asset(&path, 2).is_err());
    let a = read_asset(&path, 3).unwrap();
    assert_eq!(
        a.sha256(),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    std::fs::write(&path, b"abd").unwrap();
    let b = read_asset(&path, 3).unwrap();
    assert_ne!(a.sha256(), b.sha256());
    assert_eq!(a.bytes(), b"abc");
    std::fs::remove_file(path).unwrap();
}
#[test]
fn synthetic_profile_is_explicit_and_retains_permutation() {
    let t = gaussian_fixture(asset("ROW_ID EVENT_ID MAG\nB event 3\nA event 2\n"), 2).unwrap();
    assert_eq!(t.measurement_ids, vec!["B", "A"]);
    assert_eq!(t.event_ids, vec!["event", "event"]);
    assert!(t.zhd.is_empty());
    assert!(gaussian_fixture(asset("ROW_ID EVENT_ID LENGTH\nA event 2"), 2).is_err());
}
