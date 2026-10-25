//! The live editor over real HTTP.

use std::io::{Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};

const BOARD: &str = r#"{
  "format": "fluid-project/1",
  "name": "board",
  "params": {
    "w": { "value": "10 mm", "tol": "0.1 mm", "owner": "electrical" },
    "i": { "value": "15 A", "owner": "electrical" },
    "temp": { "value": "85 degC", "owner": "thermal" },
    "rth": { "value": "30 K/W", "owner": "thermal" }
  },
  "derived": { "rise": { "expr": "trace.power * rth", "owner": "thermal", "unit": "K" } },
  "analyses": {
    "trace": { "kind": "pcb_trace", "owner": "electrical",
      "inputs": { "width": "w", "thickness": "35[um]", "length": "40[mm]", "current": "i", "temperature": "temp" },
      "materials": { "conductor": "copper" } }
  },
  "requirements": {
    "DROP": { "text": "drop", "of": "trace.voltage_drop", "max": "50 mV", "owner": "electrical" },
    "RISE": { "text": "local rise", "of": "rise", "max": "30 K", "owner": "thermal" }
  }
}"#;

fn http(addr: SocketAddr, method: &str, path: &str, body: &str) -> String {
    let mut s = TcpStream::connect(addr).unwrap();
    write!(s, "{method} {path} HTTP/1.1\r\nHost: x\r\nContent-Length: {}\r\n\r\n{body}", body.len()).unwrap();
    let mut out = String::new();
    s.read_to_string(&mut out).unwrap();
    out.split_once("\r\n\r\n").map(|(_, b)| b.to_string()).unwrap_or_default()
}

#[test]
fn edits_are_checked_saved_and_broadcast() {
    let dir = std::env::temp_dir().join(format!("project-serve-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("project.json");
    std::fs::write(&file, BOARD).unwrap();
    let addr = project::serve::start(&file, TcpListener::bind("127.0.0.1:0").unwrap()).unwrap();

    let page = http(addr, "GET", "/", "");
    assert!(page.contains(r#"class="edit" data-id="w""#), "editable parameters");
    assert!(page.contains("let version = 0;"));
    let state = |a| http(a, "GET", "/api/state", "");
    assert!(state(addr).contains(r#""version":0"#));

    // Electrical narrows the trace: thermal's requirement fails, and the
    // reply says so.
    let r = http(addr, "POST", "/api/set", r#"{"id":"w","value":"4 mm","who":"Ana (electrical)"}"#);
    assert!(r.contains(r#""ok":true"#), "{r}");
    assert!(r.contains("RISE [thermal] met → FAILED"), "{r}");
    assert!(std::fs::read_to_string(&file).unwrap().contains(r#""value": "4 mm""#));
    assert!(state(addr).contains(r#""version":1"#));
    let page = http(addr, "GET", "/", "");
    assert!(page.contains("Ana (electrical)") && page.contains("set w (electrical) 10 mm → 4 mm"), "activity feed");

    // A wrong unit is refused, and nothing is written.
    let saved = std::fs::read_to_string(&file).unwrap();
    let r = http(addr, "POST", "/api/set", r#"{"id":"w","value":"4 degC","who":"x"}"#);
    assert!(r.contains(r#""ok":false"#), "{r}");
    let r = http(addr, "POST", "/api/set", r#"{"id":"nope","value":"1 mm","who":"x"}"#);
    assert!(r.contains("not a parameter"), "{r}");
    assert_eq!(std::fs::read_to_string(&file).unwrap(), saved);
    assert!(state(addr).contains(r#""version":1"#));

    // An edit made outside the editor (a git pull, a text editor) is seen.
    std::thread::sleep(std::time::Duration::from_millis(20));
    std::fs::write(&file, saved.replace(r#""value": "30 K/W""#, r#""value": "10 K/W""#)).unwrap();
    let mut seen = false;
    for _ in 0..50 {
        if state(addr).contains(r#""version":2"#) {
            seen = true;
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    assert!(seen, "external edit picked up");
    assert!(http(addr, "GET", "/", "").contains("changed the project file outside the editor"));

    // The contributions behind a requirement, on demand.
    let why = http(addr, "GET", "/api/why?req=DROP", "");
    assert!(why.contains("Uncertainty comes from") && why.contains("w <span>(electrical)"), "{why}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn base64_matches_the_standard() {
    // RFC 4648 test vectors.
    for (raw, enc) in [("", ""), ("f", "Zg=="), ("fo", "Zm8="), ("foo", "Zm9v"), ("foob", "Zm9vYg=="), ("fooba", "Zm9vYmE="), ("foobar", "Zm9vYmFy")] {
        assert_eq!(project::render::base64(raw.as_bytes()), enc);
    }
}
