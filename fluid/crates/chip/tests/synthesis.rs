//! Synthesis is checked against the RTL, exhaustively where the space is
//! small enough to be a proof.

use chip::equiv::check;
use chip::rtl::parse;
use chip::synth::synthesize;

fn proven(src: &str) -> chip::synth::Netlist {
    let m = parse(src).unwrap_or_else(|e| panic!("{e}"));
    let n = synthesize(&m).unwrap_or_else(|e| panic!("{e}"));
    let eq = check(&m, &n, 200_000).unwrap();
    assert!(eq.mismatch.is_none(), "{}", eq.mismatch.unwrap());
    println!("{}: {:?}, {} transistors, {} points{}", n.name, n.cells(), n.transistors(), eq.points, if eq.exhaustive { " (exhaustive)" } else { "" });
    n
}

#[test]
fn accumulator_alu_is_proven_equivalent() {
    let n = proven(
        "module acc8(input clk, input rst, input [7:0] x, input [1:0] op,
                     output [7:0] y, output reg [7:0] acc);
           wire [7:0] sum = acc + x;
           assign y = (op == 2'd0) ? sum : (op == 2'd1 ? acc & x : (op == 2'd2 ? acc ^ x : acc - x));
           always @(posedge clk) begin
             if (rst) acc <= 8'd0;
             else acc <= y;
           end
         endmodule",
    );
    assert_eq!(n.flops.len(), 8);
    assert!(n.cells().get("FA").copied().unwrap_or(0) >= 7, "the adder maps to full-adder cells");
}

#[test]
fn arithmetic_shifts_comparisons_and_case() {
    proven(
        "module mul4(input [3:0] a, input [3:0] b, output [7:0] p);
           assign p = a * b;
         endmodule",
    );
    proven(
        "module shift(input [7:0] a, input [2:0] s, input dir, output [7:0] y);
           assign y = dir ? a >> s : a << s;
         endmodule",
    );
    proven(
        "module cmp(input [5:0] a, input [5:0] b, output lt, output le, output gt, output ge, output eq, output [5:0] mx);
           assign lt = a < b; assign le = a <= b; assign gt = a > b; assign ge = a >= b; assign eq = a == b;
           assign mx = (a > b) ? a : b;
         endmodule",
    );
    proven(
        "module bits(input [7:0] a, output p, output z, output all, output [3:0] hi, output [11:0] cat, output [7:0] neg);
           assign p = ^a; assign z = !a; assign all = &a; assign hi = a[7:4];
           assign cat = {a[3:0], {2{a[1:0]}}, 2'b10}; assign neg = -a;
         endmodule",
    );
    // A traffic-light controller: case on the state, a timer, a request.
    proven(
        "module light(input clk, input rst, input req, output [1:0] lamp, output reg [1:0] st, output reg [3:0] t);
           localparam GREEN = 2'd0, AMBER = 2'd1, RED = 2'd2;
           assign lamp = st;
           always @(posedge clk) begin
             if (rst) begin st <= GREEN; t <= 4'd0; end
             else begin
               t <= t + 4'd1;
               case (st)
                 GREEN: if (req && t >= 4'd5) begin st <= AMBER; t <= 4'd0; end
                 AMBER: if (t == 4'd2) begin st <= RED; t <= 4'd0; end
                 RED:   if (t == 4'd7) begin st <= GREEN; t <= 4'd0; end
                 default: st <= GREEN;
               endcase
             end
           end
         endmodule",
    );
}

#[test]
fn a_wrong_netlist_is_caught() {
    // Swap one gate's kind in a proven design: the checker must notice.
    let m = parse(
        "module f(input [3:0] a, input [3:0] b, output [3:0] s, output c);
           assign {c, s} = a + b;
         endmodule",
    )
    .unwrap();
    let mut n = synthesize(&m).unwrap();
    assert!(check(&m, &n, 0).unwrap().mismatch.is_none());
    use semiconductor::gates::CellKind;
    let i = n.comb.insts.iter().position(|i| i.kind == CellKind::Xor2 || i.kind == CellKind::Nand2).expect("a two-input gate");
    n.comb.insts[i].kind = if n.comb.insts[i].kind == CellKind::Xor2 { CellKind::Nor2 } else { CellKind::Nor2 };
    let eq = check(&m, &n, 0).unwrap();
    assert!(eq.mismatch.is_some(), "a mutated netlist passed");
}

#[test]
fn unsupported_constructs_are_refused() {
    let bad = [
        ("module a(input clk, input r, output reg q); always @(posedge clk or posedge r) q <= 1'b0; endmodule", "asynchronous"),
        ("module a(input x, output reg q); always @* q = x; endmodule", "combinational always"),
        ("module a(input c1, input c2, input x, output reg p, output reg q); always @(posedge c1) p <= x; always @(posedge c2) q <= x; endmodule", "second clock"),
        ("module a(input x, output y); initial y = 0; endmodule", "initial"),
        ("module a(input clk, input x, output reg q); always @(posedge clk) q = x; endmodule", "non-blocking"),
        ("module a(input [3:0] x, output y); assign y = x == 4'bx01z; endmodule", "x/z"),
        ("module a(input x, output y); assign y = z; endmodule", "unknown signal"),
        ("module a(input x, output y, output w); assign y = x; endmodule", "never assigned"),
        ("module a(input x, output y); wire w; assign w = y; assign y = w; endmodule", "loop"),
    ];
    for (src, want) in bad {
        let r = parse(src).and_then(|m| synthesize(&m).map(|_| ()));
        let e = r.expect_err(src);
        assert!(e.contains(want), "{src}\n  → {e}");
    }
}
