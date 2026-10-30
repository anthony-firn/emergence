// mini8: an 8-bit accumulator CPU with four registers and a 16-word ROM.
//
// Instruction word {op[3:0], arg[7:0]}:
//   1 LDI k   acc <= k            2 LD r   acc <= r
//   3 ST r    r <= acc            4 ADD r  acc <= acc + r
//   5 SUB r   acc <= acc - r      6 AND r  7 OR r  8 XOR r
//   9 OUT     port <= acc         A JMP a  pc <= a
//   B JZ a    if (zero) pc <= a   0 NOP
// The ROM holds a Fibonacci program: port shows 1, 2, 3, 5, 8, 13, ...
module mini8(input clk, input rst, output [7:0] out, output [3:0] at);
  reg [3:0] pc;
  reg [7:0] acc, r0, r1, r2, r3, port;
  reg zero;
  wire [11:0] ins =
      (pc == 4'd0)  ? 12'h100 :   // LDI 0
      (pc == 4'd1)  ? 12'h300 :   // ST r0
      (pc == 4'd2)  ? 12'h101 :   // LDI 1
      (pc == 4'd3)  ? 12'h301 :   // ST r1
      (pc == 4'd4)  ? 12'h200 :   // LD r0
      (pc == 4'd5)  ? 12'h401 :   // ADD r1
      (pc == 4'd6)  ? 12'h900 :   // OUT
      (pc == 4'd7)  ? 12'h302 :   // ST r2
      (pc == 4'd8)  ? 12'h201 :   // LD r1
      (pc == 4'd9)  ? 12'h300 :   // ST r0
      (pc == 4'd10) ? 12'h202 :   // LD r2
      (pc == 4'd11) ? 12'h301 :   // ST r1
      (pc == 4'd12) ? 12'hA04 :   // JMP 4
                      12'h000;    // NOP
  wire [3:0] op = ins[11:8];
  wire [7:0] arg = ins[7:0];
  wire [1:0] sel = ins[1:0];
  wire [7:0] rv = (sel == 2'd0) ? r0 : (sel == 2'd1) ? r1 : (sel == 2'd2) ? r2 : r3;
  wire [7:0] alu = (op == 4'h1) ? arg :
                   (op == 4'h2) ? rv :
                   (op == 4'h4) ? acc + rv :
                   (op == 4'h5) ? acc - rv :
                   (op == 4'h6) ? (acc & rv) :
                   (op == 4'h7) ? (acc | rv) : (acc ^ rv);
  wire writes_acc = (op == 4'h1) || (op == 4'h2) || ((op >= 4'h4) && (op <= 4'h8));
  assign out = port;
  assign at = pc;
  always @(posedge clk) begin
    if (rst) begin
      pc <= 4'd0; acc <= 8'd0; zero <= 1'b0; port <= 8'd0;
      r0 <= 8'd0; r1 <= 8'd0; r2 <= 8'd0; r3 <= 8'd0;
    end else begin
      pc <= pc + 4'd1;
      if (writes_acc) begin
        acc <= alu;
        zero <= (alu == 8'd0);
      end
      case (op)
        4'h3: case (sel)
                2'd0: r0 <= acc;
                2'd1: r1 <= acc;
                2'd2: r2 <= acc;
                default: r3 <= acc;
              endcase
        4'h9: port <= acc;
        4'hA: pc <= arg[3:0];
        4'hB: if (zero) pc <= arg[3:0];
      endcase
    end
  end
endmodule
