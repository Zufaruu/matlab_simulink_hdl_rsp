// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\RADIX22FFT_SDNF1_1.v
// Created: 2023-09-27 04:10:35
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: RADIX22FFT_SDNF1_1
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/RADIX22FFT_SDNF1_1
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module RADIX22FFT_SDNF1_1
          (clk,
           reset,
           enb,
           twdlXdin_1_re,
           twdlXdin_1_im,
           twdlXdin_17_re,
           twdlXdin_17_im,
           twdlXdin_1_vld,
           dout_1_re,
           dout_1_im,
           dout_2_re,
           dout_2_im,
           dout_1_vld);


  input   clk;
  input   reset;
  input   enb;
  input   signed [12:0] twdlXdin_1_re;  // sfix13_En10
  input   signed [12:0] twdlXdin_1_im;  // sfix13_En10
  input   signed [12:0] twdlXdin_17_re;  // sfix13_En10
  input   signed [12:0] twdlXdin_17_im;  // sfix13_En10
  input   twdlXdin_1_vld;
  output  signed [12:0] dout_1_re;  // sfix13_En10
  output  signed [12:0] dout_1_im;  // sfix13_En10
  output  signed [12:0] dout_2_re;  // sfix13_En10
  output  signed [12:0] dout_2_im;  // sfix13_En10
  output  dout_1_vld;


  reg signed [13:0] Radix22ButterflyG1_NF_btf1_re_reg;  // sfix14
  reg signed [13:0] Radix22ButterflyG1_NF_btf1_im_reg;  // sfix14
  reg signed [13:0] Radix22ButterflyG1_NF_btf2_re_reg;  // sfix14
  reg signed [13:0] Radix22ButterflyG1_NF_btf2_im_reg;  // sfix14
  reg  Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  reg signed [13:0] Radix22ButterflyG1_NF_btf1_re_reg_next;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_btf1_im_reg_next;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_btf2_re_reg_next;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_btf2_im_reg_next;  // sfix14_En10
  reg  Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next;
  reg signed [12:0] dout_1_re_1;  // sfix13_En10
  reg signed [12:0] dout_1_im_1;  // sfix13_En10
  reg signed [12:0] dout_2_re_1;  // sfix13_En10
  reg signed [12:0] dout_2_im_1;  // sfix13_En10
  reg  dout_1_vld_1;
  reg signed [13:0] Radix22ButterflyG1_NF_add_cast;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_add_cast_0;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_sub_cast;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_sub_cast_0;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_add_cast_1;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_add_cast_2;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_sub_cast_1;  // sfix14_En10
  reg signed [13:0] Radix22ButterflyG1_NF_sub_cast_2;  // sfix14_En10


  // Radix22ButterflyG1_NF
  always @(posedge clk or posedge reset)
    begin : Radix22ButterflyG1_NF_process
      if (reset == 1'b1) begin
        Radix22ButterflyG1_NF_btf1_re_reg <= 14'sb00000000000000;
        Radix22ButterflyG1_NF_btf1_im_reg <= 14'sb00000000000000;
        Radix22ButterflyG1_NF_btf2_re_reg <= 14'sb00000000000000;
        Radix22ButterflyG1_NF_btf2_im_reg <= 14'sb00000000000000;
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
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, twdlXdin_17_im, twdlXdin_17_re,
       twdlXdin_1_im, twdlXdin_1_re, twdlXdin_1_vld) begin
    Radix22ButterflyG1_NF_add_cast = 14'sb00000000000000;
    Radix22ButterflyG1_NF_add_cast_0 = 14'sb00000000000000;
    Radix22ButterflyG1_NF_sub_cast = 14'sb00000000000000;
    Radix22ButterflyG1_NF_sub_cast_0 = 14'sb00000000000000;
    Radix22ButterflyG1_NF_add_cast_1 = 14'sb00000000000000;
    Radix22ButterflyG1_NF_add_cast_2 = 14'sb00000000000000;
    Radix22ButterflyG1_NF_sub_cast_1 = 14'sb00000000000000;
    Radix22ButterflyG1_NF_sub_cast_2 = 14'sb00000000000000;
    Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_btf1_re_reg;
    Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_btf1_im_reg;
    Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_btf2_re_reg;
    Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_btf2_im_reg;
    Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next = twdlXdin_1_vld;
    if (twdlXdin_1_vld) begin
      Radix22ButterflyG1_NF_add_cast = {twdlXdin_1_re[12], twdlXdin_1_re};
      Radix22ButterflyG1_NF_add_cast_0 = {twdlXdin_17_re[12], twdlXdin_17_re};
      Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_add_cast + Radix22ButterflyG1_NF_add_cast_0;
      Radix22ButterflyG1_NF_sub_cast = {twdlXdin_1_re[12], twdlXdin_1_re};
      Radix22ButterflyG1_NF_sub_cast_0 = {twdlXdin_17_re[12], twdlXdin_17_re};
      Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_sub_cast - Radix22ButterflyG1_NF_sub_cast_0;
      Radix22ButterflyG1_NF_add_cast_1 = {twdlXdin_1_im[12], twdlXdin_1_im};
      Radix22ButterflyG1_NF_add_cast_2 = {twdlXdin_17_im[12], twdlXdin_17_im};
      Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_add_cast_1 + Radix22ButterflyG1_NF_add_cast_2;
      Radix22ButterflyG1_NF_sub_cast_1 = {twdlXdin_1_im[12], twdlXdin_1_im};
      Radix22ButterflyG1_NF_sub_cast_2 = {twdlXdin_17_im[12], twdlXdin_17_im};
      Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_sub_cast_1 - Radix22ButterflyG1_NF_sub_cast_2;
    end
    dout_1_re_1 = Radix22ButterflyG1_NF_btf1_re_reg[12:0];
    dout_1_im_1 = Radix22ButterflyG1_NF_btf1_im_reg[12:0];
    dout_2_re_1 = Radix22ButterflyG1_NF_btf2_re_reg[12:0];
    dout_2_im_1 = Radix22ButterflyG1_NF_btf2_im_reg[12:0];
    dout_1_vld_1 = Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  end



  assign dout_1_re = dout_1_re_1;

  assign dout_1_im = dout_1_im_1;

  assign dout_2_re = dout_2_re_1;

  assign dout_2_im = dout_2_im_1;

  assign dout_1_vld = dout_1_vld_1;

endmodule  // RADIX22FFT_SDNF1_1
