// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\RADIX22FFT_SDNF1_1_block10.v
// Created: 2023-09-29 06:09:50
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: RADIX22FFT_SDNF1_1_block10
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/RADIX22FFT_SDNF1_1
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module RADIX22FFT_SDNF1_1_block10
          (clk,
           reset,
           enb,
           twdlXdin_12_re,
           twdlXdin_12_im,
           twdlXdin_28_re,
           twdlXdin_28_im,
           twdlXdin_1_vld,
           dout_23_re,
           dout_23_im,
           dout_24_re,
           dout_24_im);


  input   clk;
  input   reset;
  input   enb;
  input   signed [12:0] twdlXdin_12_re;  // sfix13_En10
  input   signed [12:0] twdlXdin_12_im;  // sfix13_En10
  input   signed [12:0] twdlXdin_28_re;  // sfix13_En10
  input   signed [12:0] twdlXdin_28_im;  // sfix13_En10
  input   twdlXdin_1_vld;
  output  signed [12:0] dout_23_re;  // sfix13_En10
  output  signed [12:0] dout_23_im;  // sfix13_En10
  output  signed [12:0] dout_24_re;  // sfix13_En10
  output  signed [12:0] dout_24_im;  // sfix13_En10


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
  reg signed [12:0] dout_23_re_1;  // sfix13_En10
  reg signed [12:0] dout_23_im_1;  // sfix13_En10
  reg signed [12:0] dout_24_re_1;  // sfix13_En10
  reg signed [12:0] dout_24_im_1;  // sfix13_En10
  reg  dout_23_vld;
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
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, twdlXdin_12_im, twdlXdin_12_re,
       twdlXdin_1_vld, twdlXdin_28_im, twdlXdin_28_re) begin
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
      Radix22ButterflyG1_NF_add_cast = {twdlXdin_12_re[12], twdlXdin_12_re};
      Radix22ButterflyG1_NF_add_cast_0 = {twdlXdin_28_re[12], twdlXdin_28_re};
      Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_add_cast + Radix22ButterflyG1_NF_add_cast_0;
      Radix22ButterflyG1_NF_sub_cast = {twdlXdin_12_re[12], twdlXdin_12_re};
      Radix22ButterflyG1_NF_sub_cast_0 = {twdlXdin_28_re[12], twdlXdin_28_re};
      Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_sub_cast - Radix22ButterflyG1_NF_sub_cast_0;
      Radix22ButterflyG1_NF_add_cast_1 = {twdlXdin_12_im[12], twdlXdin_12_im};
      Radix22ButterflyG1_NF_add_cast_2 = {twdlXdin_28_im[12], twdlXdin_28_im};
      Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_add_cast_1 + Radix22ButterflyG1_NF_add_cast_2;
      Radix22ButterflyG1_NF_sub_cast_1 = {twdlXdin_12_im[12], twdlXdin_12_im};
      Radix22ButterflyG1_NF_sub_cast_2 = {twdlXdin_28_im[12], twdlXdin_28_im};
      Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_sub_cast_1 - Radix22ButterflyG1_NF_sub_cast_2;
    end
    dout_23_re_1 = Radix22ButterflyG1_NF_btf1_re_reg[12:0];
    dout_23_im_1 = Radix22ButterflyG1_NF_btf1_im_reg[12:0];
    dout_24_re_1 = Radix22ButterflyG1_NF_btf2_re_reg[12:0];
    dout_24_im_1 = Radix22ButterflyG1_NF_btf2_im_reg[12:0];
    dout_23_vld = Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  end



  assign dout_23_re = dout_23_re_1;

  assign dout_23_im = dout_23_im_1;

  assign dout_24_re = dout_24_re_1;

  assign dout_24_im = dout_24_im_1;

endmodule  // RADIX22FFT_SDNF1_1_block10

