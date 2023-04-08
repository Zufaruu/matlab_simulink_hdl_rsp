// -------------------------------------------------------------
// 
// File Name: hdl_prj\hdlsrc\pemancar_bpsk\Sine_Wave1.v
// Created: 2023-04-05 16:09:07
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: Sine_Wave1
// Source Path: pemancar_bpsk/Transmitter/Sine Wave1
// Hierarchy Level: 1
// 
// Sine Wave
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module Sine_Wave1
          (clk,
           reset,
           enb,
           Sine_Wave_out1);


  input   clk;
  input   reset;
  input   enb;
  output  signed [15:0] Sine_Wave_out1;  // sfix16_En14


  reg  address_cnt1;  // ufix1
  wire signed [15:0] tableData0;  // sfix16_En14
  wire signed [15:0] tableData1;  // sfix16_En14
  wire signed [15:0] lut1out1;  // sfix16_En14


  // Count limited, Unsigned Counter
  //  initial value   = 0
  //  step value      = 1
  //  count to value  = 1
  always @(posedge clk or posedge reset)
    begin : Sine_Wave_addrcnt_temp_process1_process
      if (reset == 1'b1) begin
        address_cnt1 <= 1'b0;
      end
      else begin
        if (enb) begin
          address_cnt1 <=  ~ address_cnt1;
        end
      end
    end



  assign tableData0 = 16'sb0100000000000000;



  assign tableData1 = 16'sb1100000000000000;



  assign lut1out1 = (address_cnt1 == 1'b0 ? tableData0 :
              tableData1);



  assign Sine_Wave_out1 = lut1out1;

endmodule  // Sine_Wave1

