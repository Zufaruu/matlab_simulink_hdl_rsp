// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\RADIX22FFT_SDNF1_3_block14.v
// Created: 2023-09-29 06:09:51
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: RADIX22FFT_SDNF1_3_block14
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/RADIX22FFT_SDNF1_3
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module RADIX22FFT_SDNF1_3_block14
          (clk,
           reset,
           enb,
           twdlXdin_28_re,
           twdlXdin_28_im,
           twdlXdin_32_re,
           twdlXdin_32_im,
           twdlXdin_1_vld,
           dout_31_re,
           dout_31_im,
           dout_32_re,
           dout_32_im);


  input   clk;
  input   reset;
  input   enb;
  input   signed [14:0] twdlXdin_28_re;  // sfix15_En10
  input   signed [14:0] twdlXdin_28_im;  // sfix15_En10
  input   signed [14:0] twdlXdin_32_re;  // sfix15_En10
  input   signed [14:0] twdlXdin_32_im;  // sfix15_En10
  input   twdlXdin_1_vld;
  output  signed [14:0] dout_31_re;  // sfix15_En10
  output  signed [14:0] dout_31_im;  // sfix15_En10
  output  signed [14:0] dout_32_re;  // sfix15_En10
  output  signed [14:0] dout_32_im;  // sfix15_En10


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
  reg signed [14:0] dout_31_re_1;  // sfix15_En10
  reg signed [14:0] dout_31_im_1;  // sfix15_En10
  reg signed [14:0] dout_32_re_1;  // sfix15_En10
  reg signed [14:0] dout_32_im_1;  // sfix15_En10
  reg  dout_31_vld;
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
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, twdlXdin_1_vld, twdlXdin_28_im,
       twdlXdin_28_re, twdlXdin_32_im, twdlXdin_32_re) begin
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
      Radix22ButterflyG1_NF_add_cast = {twdlXdin_28_re[14], twdlXdin_28_re};
      Radix22ButterflyG1_NF_add_cast_0 = {twdlXdin_32_re[14], twdlXdin_32_re};
      Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_add_cast + Radix22ButterflyG1_NF_add_cast_0;
      Radix22ButterflyG1_NF_sub_cast = {twdlXdin_28_re[14], twdlXdin_28_re};
      Radix22ButterflyG1_NF_sub_cast_0 = {twdlXdin_32_re[14], twdlXdin_32_re};
      Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_sub_cast - Radix22ButterflyG1_NF_sub_cast_0;
      Radix22ButterflyG1_NF_add_cast_1 = {twdlXdin_28_im[14], twdlXdin_28_im};
      Radix22ButterflyG1_NF_add_cast_2 = {twdlXdin_32_im[14], twdlXdin_32_im};
      Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_add_cast_1 + Radix22ButterflyG1_NF_add_cast_2;
      Radix22ButterflyG1_NF_sub_cast_1 = {twdlXdin_28_im[14], twdlXdin_28_im};
      Radix22ButterflyG1_NF_sub_cast_2 = {twdlXdin_32_im[14], twdlXdin_32_im};
      Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_sub_cast_1 - Radix22ButterflyG1_NF_sub_cast_2;
    end
    dout_31_re_1 = Radix22ButterflyG1_NF_btf1_re_reg[14:0];
    dout_31_im_1 = Radix22ButterflyG1_NF_btf1_im_reg[14:0];
    dout_32_re_1 = Radix22ButterflyG1_NF_btf2_re_reg[14:0];
    dout_32_im_1 = Radix22ButterflyG1_NF_btf2_im_reg[14:0];
    dout_31_vld = Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  end



  assign dout_31_re = dout_31_re_1;

  assign dout_31_im = dout_31_im_1;

  assign dout_32_re = dout_32_re_1;

  assign dout_32_im = dout_32_im_1;

endmodule  // RADIX22FFT_SDNF1_3_block14

