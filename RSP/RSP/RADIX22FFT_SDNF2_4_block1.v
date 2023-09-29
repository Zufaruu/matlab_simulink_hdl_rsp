// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\RADIX22FFT_SDNF2_4_block1.v
// Created: 2023-09-29 06:09:51
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: RADIX22FFT_SDNF2_4_block1
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/RADIX22FFT_SDNF2_4
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module RADIX22FFT_SDNF2_4_block1
          (clk,
           reset,
           enb,
           rotate_5,
           dout_2_re,
           dout_2_im,
           dout_6_re,
           dout_6_im,
           dout_1_vld,
           dout_5_re,
           dout_5_im,
           dout_6_re_1,
           dout_6_im_1);


  input   clk;
  input   reset;
  input   enb;
  input   rotate_5;  // ufix1
  input   signed [14:0] dout_2_re;  // sfix15_En10
  input   signed [14:0] dout_2_im;  // sfix15_En10
  input   signed [14:0] dout_6_re;  // sfix15_En10
  input   signed [14:0] dout_6_im;  // sfix15_En10
  input   dout_1_vld;
  output  signed [15:0] dout_5_re;  // sfix16_En10
  output  signed [15:0] dout_5_im;  // sfix16_En10
  output  signed [15:0] dout_6_re_1;  // sfix16_En10
  output  signed [15:0] dout_6_im_1;  // sfix16_En10


  wire signed [15:0] din1_re;  // sfix16_En10
  wire signed [15:0] din1_im;  // sfix16_En10
  wire signed [15:0] din2_re;  // sfix16_En10
  wire signed [15:0] din2_im;  // sfix16_En10
  reg  Radix22ButterflyG2_NF_din_vld_dly;
  reg signed [16:0] Radix22ButterflyG2_NF_btf1_re_reg;  // sfix17
  reg signed [16:0] Radix22ButterflyG2_NF_btf1_im_reg;  // sfix17
  reg signed [16:0] Radix22ButterflyG2_NF_btf2_re_reg;  // sfix17
  reg signed [16:0] Radix22ButterflyG2_NF_btf2_im_reg;  // sfix17
  reg  Radix22ButterflyG2_NF_din_vld_dly_next;
  reg signed [16:0] Radix22ButterflyG2_NF_btf1_re_reg_next;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_btf1_im_reg_next;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_btf2_re_reg_next;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_btf2_im_reg_next;  // sfix17_En10
  reg signed [15:0] dout_5_re_1;  // sfix16_En10
  reg signed [15:0] dout_5_im_1;  // sfix16_En10
  reg signed [15:0] dout_6_re_2;  // sfix16_En10
  reg signed [15:0] dout_6_im_2;  // sfix16_En10
  reg  dout_4_vld;
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_0;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_1;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_2;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_0;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_1;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_2;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_3;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_4;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_5;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_add_cast_6;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_3;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_4;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_5;  // sfix17_En10
  reg signed [16:0] Radix22ButterflyG2_NF_sub_cast_6;  // sfix17_En10


  assign din1_re = {dout_2_re[14], dout_2_re};



  assign din1_im = {dout_2_im[14], dout_2_im};



  assign din2_re = {dout_6_re[14], dout_6_re};



  assign din2_im = {dout_6_im[14], dout_6_im};



  // Radix22ButterflyG2_NF
  always @(posedge clk or posedge reset)
    begin : Radix22ButterflyG2_NF_process
      if (reset == 1'b1) begin
        Radix22ButterflyG2_NF_din_vld_dly <= 1'b0;
        Radix22ButterflyG2_NF_btf1_re_reg <= 17'sb00000000000000000;
        Radix22ButterflyG2_NF_btf1_im_reg <= 17'sb00000000000000000;
        Radix22ButterflyG2_NF_btf2_re_reg <= 17'sb00000000000000000;
        Radix22ButterflyG2_NF_btf2_im_reg <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          Radix22ButterflyG2_NF_din_vld_dly <= Radix22ButterflyG2_NF_din_vld_dly_next;
          Radix22ButterflyG2_NF_btf1_re_reg <= Radix22ButterflyG2_NF_btf1_re_reg_next;
          Radix22ButterflyG2_NF_btf1_im_reg <= Radix22ButterflyG2_NF_btf1_im_reg_next;
          Radix22ButterflyG2_NF_btf2_re_reg <= Radix22ButterflyG2_NF_btf2_re_reg_next;
          Radix22ButterflyG2_NF_btf2_im_reg <= Radix22ButterflyG2_NF_btf2_im_reg_next;
        end
      end
    end

  always @(Radix22ButterflyG2_NF_btf1_im_reg, Radix22ButterflyG2_NF_btf1_re_reg,
       Radix22ButterflyG2_NF_btf2_im_reg, Radix22ButterflyG2_NF_btf2_re_reg,
       Radix22ButterflyG2_NF_din_vld_dly, din1_im, din1_re, din2_im, din2_re,
       dout_1_vld, rotate_5) begin
    Radix22ButterflyG2_NF_add_cast_1 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast_2 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_1 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_2 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast_5 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast_6 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_5 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_6 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast_0 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_0 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast_3 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_add_cast_4 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_3 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_sub_cast_4 = 17'sb00000000000000000;
    Radix22ButterflyG2_NF_btf1_re_reg_next = Radix22ButterflyG2_NF_btf1_re_reg;
    Radix22ButterflyG2_NF_btf1_im_reg_next = Radix22ButterflyG2_NF_btf1_im_reg;
    Radix22ButterflyG2_NF_btf2_re_reg_next = Radix22ButterflyG2_NF_btf2_re_reg;
    Radix22ButterflyG2_NF_btf2_im_reg_next = Radix22ButterflyG2_NF_btf2_im_reg;
    Radix22ButterflyG2_NF_din_vld_dly_next = dout_1_vld;
    if (rotate_5 != 1'b0) begin
      if (dout_1_vld) begin
        Radix22ButterflyG2_NF_add_cast_1 = {din1_re[15], din1_re};
        Radix22ButterflyG2_NF_add_cast_2 = {din2_im[15], din2_im};
        Radix22ButterflyG2_NF_btf1_re_reg_next = Radix22ButterflyG2_NF_add_cast_1 + Radix22ButterflyG2_NF_add_cast_2;
        Radix22ButterflyG2_NF_sub_cast_1 = {din1_re[15], din1_re};
        Radix22ButterflyG2_NF_sub_cast_2 = {din2_im[15], din2_im};
        Radix22ButterflyG2_NF_btf2_re_reg_next = Radix22ButterflyG2_NF_sub_cast_1 - Radix22ButterflyG2_NF_sub_cast_2;
        Radix22ButterflyG2_NF_add_cast_5 = {din1_im[15], din1_im};
        Radix22ButterflyG2_NF_add_cast_6 = {din2_re[15], din2_re};
        Radix22ButterflyG2_NF_btf2_im_reg_next = Radix22ButterflyG2_NF_add_cast_5 + Radix22ButterflyG2_NF_add_cast_6;
        Radix22ButterflyG2_NF_sub_cast_5 = {din1_im[15], din1_im};
        Radix22ButterflyG2_NF_sub_cast_6 = {din2_re[15], din2_re};
        Radix22ButterflyG2_NF_btf1_im_reg_next = Radix22ButterflyG2_NF_sub_cast_5 - Radix22ButterflyG2_NF_sub_cast_6;
      end
    end
    else if (dout_1_vld) begin
      Radix22ButterflyG2_NF_add_cast = {din1_re[15], din1_re};
      Radix22ButterflyG2_NF_add_cast_0 = {din2_re[15], din2_re};
      Radix22ButterflyG2_NF_btf1_re_reg_next = Radix22ButterflyG2_NF_add_cast + Radix22ButterflyG2_NF_add_cast_0;
      Radix22ButterflyG2_NF_sub_cast = {din1_re[15], din1_re};
      Radix22ButterflyG2_NF_sub_cast_0 = {din2_re[15], din2_re};
      Radix22ButterflyG2_NF_btf2_re_reg_next = Radix22ButterflyG2_NF_sub_cast - Radix22ButterflyG2_NF_sub_cast_0;
      Radix22ButterflyG2_NF_add_cast_3 = {din1_im[15], din1_im};
      Radix22ButterflyG2_NF_add_cast_4 = {din2_im[15], din2_im};
      Radix22ButterflyG2_NF_btf1_im_reg_next = Radix22ButterflyG2_NF_add_cast_3 + Radix22ButterflyG2_NF_add_cast_4;
      Radix22ButterflyG2_NF_sub_cast_3 = {din1_im[15], din1_im};
      Radix22ButterflyG2_NF_sub_cast_4 = {din2_im[15], din2_im};
      Radix22ButterflyG2_NF_btf2_im_reg_next = Radix22ButterflyG2_NF_sub_cast_3 - Radix22ButterflyG2_NF_sub_cast_4;
    end
    dout_5_re_1 = Radix22ButterflyG2_NF_btf1_re_reg[15:0];
    dout_5_im_1 = Radix22ButterflyG2_NF_btf1_im_reg[15:0];
    dout_6_re_2 = Radix22ButterflyG2_NF_btf2_re_reg[15:0];
    dout_6_im_2 = Radix22ButterflyG2_NF_btf2_im_reg[15:0];
    dout_4_vld = Radix22ButterflyG2_NF_din_vld_dly;
  end



  assign dout_5_re = dout_5_re_1;

  assign dout_5_im = dout_5_im_1;

  assign dout_6_re_1 = dout_6_re_2;

  assign dout_6_im_1 = dout_6_im_2;

endmodule  // RADIX22FFT_SDNF2_4_block1

