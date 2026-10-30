/*
 * Copyright (c) 2026 HADJAL001
 * SPDX-License-Identifier: Apache-2.0
 *
 * tt_um_hadjal001_motor_controller: dual-channel 8-bit PWM motor controller,
 * from https://github.com/HADJAL001/tt-motor-controller at commit
 * 0f078f59e94b7dfc94b4718f272cc4e9d9191f84 (src/project.v), a Tiny Tapeout
 * design. Adapted to the emergence Verilog subset (analysis/chip.md):
 *   - the asynchronous active-low reset is synchronous (the subset refuses
 *     asynchronous resets): registers clear on the next clock edge with
 *     rst_n low, instead of immediately;
 *   - `default_nettype and the lint-only _unused wire are removed;
 *   - per-bit output assigns are one concatenation; `input wire` is `input`.
 * Behaviour on every clock edge is otherwise unchanged; the equivalence check
 * below compares this RTL with its synthesized netlist, not with the original.
 *
 *   ui_in[7:0]  duty value      uio_in[0] chan_select  uio_in[1] load_strobe
 *   uio_in[2]   dir_a           uio_in[3] dir_b        uio_in[4] enable
 *   uo_out[0]   pwm_a  [1] pwm_b  [2] dir_a  [3] dir_b  [4] enable
 * PWM frequency = clk / 256.
 */
module tt_um_hadjal001_motor_controller(input [7:0] ui_in, output [7:0] uo_out, input [7:0] uio_in,
    output [7:0] uio_out, output [7:0] uio_oe, input ena, input clk, input rst_n);
  wire chan_select = uio_in[0];
  wire load_strobe = uio_in[1];
  wire dir_a_in = uio_in[2];
  wire dir_b_in = uio_in[3];
  wire enable_in = uio_in[4];
  reg [7:0] duty_a;
  reg [7:0] duty_b;
  reg [7:0] pwm_counter;
  reg load_strobe_prev;
  wire load_strobe_rising = load_strobe & ~load_strobe_prev;
  always @(posedge clk) begin
    if (!rst_n) begin
      duty_a <= 8'd0;
      duty_b <= 8'd0;
      pwm_counter <= 8'd0;
      load_strobe_prev <= 1'b0;
    end else begin
      load_strobe_prev <= load_strobe;
      if (load_strobe_rising) begin
        if (chan_select) duty_b <= ui_in;
        else duty_a <= ui_in;
      end
      pwm_counter <= pwm_counter + 8'd1;
    end
  end
  wire pwm_a = enable_in && (pwm_counter < duty_a);
  wire pwm_b = enable_in && (pwm_counter < duty_b);
  assign uo_out = {3'b000, enable_in, dir_b_in, dir_a_in, pwm_b, pwm_a};
  assign uio_out = 8'd0;
  assign uio_oe = 8'd0;
endmodule
