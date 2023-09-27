// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\RADIX22FFT_SDNF1_3_block9.v
// Created: 2023-09-27 04:10:36
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: RADIX22FFT_SDNF1_3_block9
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/RADIX22FFT_SDNF1_3
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module RADIX22FFT_SDNF1_3_block9
          (clk,
           reset,
           enb,
           twdlXdin_19_re,
           twdlXdin_19_im,
           twdlXdin_23_re,
           twdlXdin_23_im,
           twdlXdin_1_vld,
           dout_21_re,
           dout_21_im,
           dout_22_re,
           dout_22_im);


  input   clk;
  input   reset;
  input   enb;
  input   signed [14:0] twdlXdin_19_re;  // sfix15_En10
  input   signed [14:0] twdlXdin_19_im;  // sfix15_En10
  input   signed [14:0] twdlXdin_23_re;  // sfix15_En10
  input   signed [14:0] twdlXdin_23_im;  // sfix15_En10
  input   twdlXdin_1_vld;
  output  signed [14:0] dout_21_re;  // sfix15_En10
  output  signed [14:0] dout_21_im;  // sfix15_En10
  output  signed [14:0] dout_22_re;  // sfix15_En10
  output  signed [14:0] dout_22_im;  // sfix15_En10


  reg signed [15:0] Radix22ButterflyG1_NF_btf1_re_reg;  // sfix16
  reg signed [15:0] Radix22ButterflyG1_NF_btf1_im_reg;  // sfix16
  reg signed [15:0] Radix22ButterflyG1_NF_btf2_re_reg;  // sfix16
  reg signed [15:0] Radix22ButterflyG1_NF_btf2_im_reg;  // sfix16
  reg  Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  reg signed [15:0] Radix22ButterflyG1_NF_btf1_re_reg_next;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_btf1_im_reg_next;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_btf2_re_reg_next;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_btf2_im_reg_next;  // sfix16_En10
  reg  Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next;
  reg signed [14:0] dout_21_re_1;  // sfix15_En10
  reg signed [14:0] dout_21_im_1;  // sfix15_En10
  reg signed [14:0] dout_22_re_1;  // sfix15_En10
  reg signed [14:0] dout_22_im_1;  // sfix15_En10
  reg  dout_21_vld;
  reg signed [15:0] Radix22ButterflyG1_NF_add_cast;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_add_cast_0;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_sub_cast;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_sub_cast_0;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_add_cast_1;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_add_cast_2;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_sub_cast_1;  // sfix16_En10
  reg signed [15:0] Radix22ButterflyG1_NF_sub_cast_2;  // sfix16_En10


  // Radix22ButterflyG1_NF
  always @(posedge clk or posedge reset)
    begin : Radix22ButterflyG1_NF_process
      if (reset == 1'b1) begin
        Radix22ButterflyG1_NF_btf1_re_reg <= 16'sb0000000000000000;
        Radix22ButterflyG1_NF_btf1_im_reg <= 16'sb0000000000000000;
        Radix22ButterflyG1_NF_btf2_re_reg <= 16'sb0000000000000000;
        Radix22ButterflyG1_NF_btf2_im_reg <= 16'sb0000000000000000;
        Radix22ButterflyG1_NF_dinXtwdl_vld_dly1 <= 1'b0;
      end
      else begin
        if (enb) begin
          Radix22ButterflyG1_NF_btf1_re_reg <= Radix22ButterflyG1_NF_btf1_re_reg_next;
          Radix22ButterflyG1_NF_btf1_im_reg <= Radix22ButterflyG1_NF_btf1_im_reg_next;
          Radix22ButterflyG1_NF_btf2_re_reg <= Radix22ButterflyG1_NF_btf2_re_reg_next;
          Radix22ButterflyG1_NF_btf2_im_reg <= Radix22ButterflyG1_NF_btf2_im_reg_next;
          Radix22ButterflyG1_NF_dinXtwdl_vld_dly1 <= Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next;
        end
      end
    end

  always @(Radix22ButterflyG1_NF_btf1_im_reg, Radix22ButterflyG1_NF_btf1_re_reg,
       Radix22ButterflyG1_NF_btf2_im_reg, Radix22ButterflyG1_NF_btf2_re_reg,
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, twdlXdin_19_im, twdlXdin_19_re,
       twdlXdin_1_vld, twdlXdin_23_im, twdlXdin_23_re) begin
    Radix22ButterflyG1_NF_add_cast = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_add_cast_0 = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_sub_cast = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_sub_cast_0 = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_add_cast_1 = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_add_cast_2 = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_sub_cast_1 = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_sub_cast_2 = 16'sb0000000000000000;
    Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_btf1_re_reg;
    Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_btf1_im_reg;
    Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_btf2_re_reg;
    Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_btf2_im_reg;
    Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next = twdlXdin_1_vld;
    if (twdlXdin_1_vld) begin
      Radix22ButterflyG1_NF_add_cast = {twdlXdin_19_re[14], twdlXdin_19_re};
      Radix22ButterflyG1_NF_add_cast_0 = {twdlXdin_23_re[14], twdlXdin_23_re};
      Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_add_cast + Radix22ButterflyG1_NF_add_cast_0;
      Radix22ButterflyG1_NF_sub_cast = {twdlXdin_19_re[14], twdlXdin_19_re};
      Radix22ButterflyG1_NF_sub_cast_0 = {twdlXdin_23_re[14], twdlXdin_23_re};
      Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_sub_cast - Radix22ButterflyG1_NF_sub_cast_0;
      Radix22ButterflyG1_NF_add_cast_1 = {twdlXdin_19_im[14], twdlXdin_19_im};
      Radix22ButterflyG1_NF_add_cast_2 = {twdlXdin_23_im[14], twdlXdin_23_im};
      Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_add_cast_1 + Radix22ButterflyG1_NF_add_cast_2;
      Radix22ButterflyG1_NF_sub_cast_1 = {twdlXdin_19_im[14], twdlXdin_19_im};
      Radix22ButterflyG1_NF_sub_cast_2 = {twdlXdin_23_im[14], twdlXdin_23_im};
      Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_sub_cast_1 - Radix22ButterflyG1_NF_sub_cast_2;
    end
    dout_21_re_1 = Radix22ButterflyG1_NF_btf1_re_reg[14:0];
    dout_21_im_1 = Radix22ButterflyG1_NF_btf1_im_reg[14:0];
    dout_22_re_1 = Radix22ButterflyG1_NF_btf2_re_reg[14:0];
    dout_22_im_1 = Radix22ButterflyG1_NF_btf2_im_reg[14:0];
    dout_21_vld = Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  end



  assign dout_21_re = dout_21_re_1;

  assign dout_21_im = dout_21_im_1;

  assign dout_22_re = dout_22_re_1;

  assign dout_22_im = dout_22_im_1;

endmodule  // RADIX22FFT_SDNF1_3_block9
