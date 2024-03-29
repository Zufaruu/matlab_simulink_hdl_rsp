// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\RADIX22FFT_SDNF1_5_block11.v
// Created: 2023-09-29 06:09:51
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: RADIX22FFT_SDNF1_5_block11
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/RADIX22FFT_SDNF1_5
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module RADIX22FFT_SDNF1_5_block11
          (clk,
           reset,
           enb,
           dout_25_re,
           dout_25_im,
           dout_27_re,
           dout_27_im,
           dout_4_vld,
           twdl_5_26_re,
           twdl_5_26_im,
           dout_25_re_1,
           dout_25_im_1,
           dout_26_re,
           dout_26_im);


  input   clk;
  input   reset;
  input   enb;
  input   signed [15:0] dout_25_re;  // sfix16_En10
  input   signed [15:0] dout_25_im;  // sfix16_En10
  input   signed [15:0] dout_27_re;  // sfix16_En10
  input   signed [15:0] dout_27_im;  // sfix16_En10
  input   dout_4_vld;
  input   signed [11:0] twdl_5_26_re;  // sfix12_En10
  input   signed [11:0] twdl_5_26_im;  // sfix12_En10
  output  signed [16:0] dout_25_re_1;  // sfix17_En10
  output  signed [16:0] dout_25_im_1;  // sfix17_En10
  output  signed [16:0] dout_26_re;  // sfix17_En10
  output  signed [16:0] dout_26_im;  // sfix17_En10


  wire signed [16:0] din_re;  // sfix17_En10
  reg signed [16:0] din1_re_dly1;  // sfix17_En10
  reg signed [16:0] din1_re_dly2;  // sfix17_En10
  reg signed [16:0] din1_re_dly3;  // sfix17_En10
  reg signed [16:0] din1_re_dly4;  // sfix17_En10
  reg signed [16:0] din1_re_dly5;  // sfix17_En10
  reg signed [16:0] din1_re_dly6;  // sfix17_En10
  reg signed [16:0] din1_re_dly7;  // sfix17_En10
  reg signed [16:0] din1_re_dly8;  // sfix17_En10
  reg signed [16:0] din1_re_dly9;  // sfix17_En10
  wire signed [16:0] din_im;  // sfix17_En10
  reg signed [16:0] din1_im_dly1;  // sfix17_En10
  reg signed [16:0] din1_im_dly2;  // sfix17_En10
  reg signed [16:0] din1_im_dly3;  // sfix17_En10
  reg signed [16:0] din1_im_dly4;  // sfix17_En10
  reg signed [16:0] din1_im_dly5;  // sfix17_En10
  reg signed [16:0] din1_im_dly6;  // sfix17_En10
  reg signed [16:0] din1_im_dly7;  // sfix17_En10
  reg signed [16:0] din1_im_dly8;  // sfix17_En10
  reg signed [16:0] din1_im_dly9;  // sfix17_En10
  wire signed [16:0] din_re_1;  // sfix17_En10
  reg signed [16:0] din2_re_dly1;  // sfix17_En10
  reg signed [16:0] din2_re_dly2;  // sfix17_En10
  reg signed [16:0] din2_re_dly3;  // sfix17_En10
  wire signed [16:0] din_im_1;  // sfix17_En10
  reg signed [16:0] din2_im_dly1;  // sfix17_En10
  reg signed [16:0] din2_im_dly2;  // sfix17_En10
  reg signed [16:0] din2_im_dly3;  // sfix17_En10
  reg  di2_vld_dly1;
  reg  di2_vld_dly2;
  reg  di2_vld_dly3;
  wire signed [16:0] dinXTwdl2_re;  // sfix17_En10
  wire signed [16:0] dinXTwdl2_im;  // sfix17_En10
  reg  din_vld_dly1;
  reg  din_vld_dly2;
  reg  din_vld_dly3;
  reg  din_vld_dly4;
  reg  din_vld_dly5;
  reg  din_vld_dly6;
  reg  din_vld_dly7;
  reg  din_vld_dly8;
  reg  din_vld_dly9;
  reg signed [17:0] Radix22ButterflyG1_NF_btf1_re_reg;  // sfix18
  reg signed [17:0] Radix22ButterflyG1_NF_btf1_im_reg;  // sfix18
  reg signed [17:0] Radix22ButterflyG1_NF_btf2_re_reg;  // sfix18
  reg signed [17:0] Radix22ButterflyG1_NF_btf2_im_reg;  // sfix18
  reg  Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  reg signed [17:0] Radix22ButterflyG1_NF_btf1_re_reg_next;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_btf1_im_reg_next;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_btf2_re_reg_next;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_btf2_im_reg_next;  // sfix18_En10
  reg  Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next;
  reg signed [16:0] dout_25_re_2;  // sfix17_En10
  reg signed [16:0] dout_25_im_2;  // sfix17_En10
  reg signed [16:0] dout_26_re_1;  // sfix17_En10
  reg signed [16:0] dout_26_im_1;  // sfix17_En10
  reg  dout_26_vld;
  reg signed [17:0] Radix22ButterflyG1_NF_add_cast;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_add_cast_0;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_sub_cast;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_sub_cast_0;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_add_cast_1;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_add_cast_2;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_sub_cast_1;  // sfix18_En10
  reg signed [17:0] Radix22ButterflyG1_NF_sub_cast_2;  // sfix18_En10


  assign din_re = {dout_25_re[15], dout_25_re};



  always @(posedge clk or posedge reset)
    begin : intdelay_process
      if (reset == 1'b1) begin
        din1_re_dly1 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly1 <= din_re;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_1_process
      if (reset == 1'b1) begin
        din1_re_dly2 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly2 <= din1_re_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_2_process
      if (reset == 1'b1) begin
        din1_re_dly3 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly3 <= din1_re_dly2;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_3_process
      if (reset == 1'b1) begin
        din1_re_dly4 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly4 <= din1_re_dly3;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_4_process
      if (reset == 1'b1) begin
        din1_re_dly5 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly5 <= din1_re_dly4;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_5_process
      if (reset == 1'b1) begin
        din1_re_dly6 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly6 <= din1_re_dly5;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_6_process
      if (reset == 1'b1) begin
        din1_re_dly7 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly7 <= din1_re_dly6;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_7_process
      if (reset == 1'b1) begin
        din1_re_dly8 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly8 <= din1_re_dly7;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_8_process
      if (reset == 1'b1) begin
        din1_re_dly9 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly9 <= din1_re_dly8;
        end
      end
    end



  assign din_im = {dout_25_im[15], dout_25_im};



  always @(posedge clk or posedge reset)
    begin : intdelay_9_process
      if (reset == 1'b1) begin
        din1_im_dly1 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly1 <= din_im;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_10_process
      if (reset == 1'b1) begin
        din1_im_dly2 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly2 <= din1_im_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_11_process
      if (reset == 1'b1) begin
        din1_im_dly3 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly3 <= din1_im_dly2;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_12_process
      if (reset == 1'b1) begin
        din1_im_dly4 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly4 <= din1_im_dly3;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_13_process
      if (reset == 1'b1) begin
        din1_im_dly5 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly5 <= din1_im_dly4;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_14_process
      if (reset == 1'b1) begin
        din1_im_dly6 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly6 <= din1_im_dly5;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_15_process
      if (reset == 1'b1) begin
        din1_im_dly7 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly7 <= din1_im_dly6;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_16_process
      if (reset == 1'b1) begin
        din1_im_dly8 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly8 <= din1_im_dly7;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_17_process
      if (reset == 1'b1) begin
        din1_im_dly9 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly9 <= din1_im_dly8;
        end
      end
    end



  assign din_re_1 = {dout_27_re[15], dout_27_re};



  always @(posedge clk or posedge reset)
    begin : intdelay_18_process
      if (reset == 1'b1) begin
        din2_re_dly1 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din2_re_dly1 <= din_re_1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_19_process
      if (reset == 1'b1) begin
        din2_re_dly2 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din2_re_dly2 <= din2_re_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_20_process
      if (reset == 1'b1) begin
        din2_re_dly3 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din2_re_dly3 <= din2_re_dly2;
        end
      end
    end



  assign din_im_1 = {dout_27_im[15], dout_27_im};



  always @(posedge clk or posedge reset)
    begin : intdelay_21_process
      if (reset == 1'b1) begin
        din2_im_dly1 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din2_im_dly1 <= din_im_1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_22_process
      if (reset == 1'b1) begin
        din2_im_dly2 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din2_im_dly2 <= din2_im_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_23_process
      if (reset == 1'b1) begin
        din2_im_dly3 <= 17'sb00000000000000000;
      end
      else begin
        if (enb) begin
          din2_im_dly3 <= din2_im_dly2;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_24_process
      if (reset == 1'b1) begin
        di2_vld_dly1 <= 1'b0;
      end
      else begin
        if (enb) begin
          di2_vld_dly1 <= dout_4_vld;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_25_process
      if (reset == 1'b1) begin
        di2_vld_dly2 <= 1'b0;
      end
      else begin
        if (enb) begin
          di2_vld_dly2 <= di2_vld_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_26_process
      if (reset == 1'b1) begin
        di2_vld_dly3 <= 1'b0;
      end
      else begin
        if (enb) begin
          di2_vld_dly3 <= di2_vld_dly2;
        end
      end
    end



  Complex4Multiply_block39 u_MUL4_2 (.clk(clk),
                                     .reset(reset),
                                     .enb(enb),
                                     .din2_re_dly3(din2_re_dly3),  // sfix17_En10
                                     .din2_im_dly3(din2_im_dly3),  // sfix17_En10
                                     .di2_vld_dly3(di2_vld_dly3),
                                     .twdl_5_26_re(twdl_5_26_re),  // sfix12_En10
                                     .twdl_5_26_im(twdl_5_26_im),  // sfix12_En10
                                     .dinXTwdl2_re(dinXTwdl2_re),  // sfix17_En10
                                     .dinXTwdl2_im(dinXTwdl2_im)  // sfix17_En10
                                     );

  always @(posedge clk or posedge reset)
    begin : intdelay_27_process
      if (reset == 1'b1) begin
        din_vld_dly1 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly1 <= dout_4_vld;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_28_process
      if (reset == 1'b1) begin
        din_vld_dly2 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly2 <= din_vld_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_29_process
      if (reset == 1'b1) begin
        din_vld_dly3 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly3 <= din_vld_dly2;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_30_process
      if (reset == 1'b1) begin
        din_vld_dly4 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly4 <= din_vld_dly3;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_31_process
      if (reset == 1'b1) begin
        din_vld_dly5 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly5 <= din_vld_dly4;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_32_process
      if (reset == 1'b1) begin
        din_vld_dly6 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly6 <= din_vld_dly5;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_33_process
      if (reset == 1'b1) begin
        din_vld_dly7 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly7 <= din_vld_dly6;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_34_process
      if (reset == 1'b1) begin
        din_vld_dly8 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly8 <= din_vld_dly7;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_35_process
      if (reset == 1'b1) begin
        din_vld_dly9 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly9 <= din_vld_dly8;
        end
      end
    end



  // Radix22ButterflyG1_NF
  always @(posedge clk or posedge reset)
    begin : Radix22ButterflyG1_NF_process
      if (reset == 1'b1) begin
        Radix22ButterflyG1_NF_btf1_re_reg <= 18'sb000000000000000000;
        Radix22ButterflyG1_NF_btf1_im_reg <= 18'sb000000000000000000;
        Radix22ButterflyG1_NF_btf2_re_reg <= 18'sb000000000000000000;
        Radix22ButterflyG1_NF_btf2_im_reg <= 18'sb000000000000000000;
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
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, din1_im_dly9, din1_re_dly9,
       dinXTwdl2_im, dinXTwdl2_re, din_vld_dly9) begin
    Radix22ButterflyG1_NF_add_cast = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_add_cast_0 = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_sub_cast = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_sub_cast_0 = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_add_cast_1 = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_add_cast_2 = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_sub_cast_1 = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_sub_cast_2 = 18'sb000000000000000000;
    Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_btf1_re_reg;
    Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_btf1_im_reg;
    Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_btf2_re_reg;
    Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_btf2_im_reg;
    Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next = din_vld_dly9;
    if (din_vld_dly9) begin
      Radix22ButterflyG1_NF_add_cast = {din1_re_dly9[16], din1_re_dly9};
      Radix22ButterflyG1_NF_add_cast_0 = {dinXTwdl2_re[16], dinXTwdl2_re};
      Radix22ButterflyG1_NF_btf1_re_reg_next = Radix22ButterflyG1_NF_add_cast + Radix22ButterflyG1_NF_add_cast_0;
      Radix22ButterflyG1_NF_sub_cast = {din1_re_dly9[16], din1_re_dly9};
      Radix22ButterflyG1_NF_sub_cast_0 = {dinXTwdl2_re[16], dinXTwdl2_re};
      Radix22ButterflyG1_NF_btf2_re_reg_next = Radix22ButterflyG1_NF_sub_cast - Radix22ButterflyG1_NF_sub_cast_0;
      Radix22ButterflyG1_NF_add_cast_1 = {din1_im_dly9[16], din1_im_dly9};
      Radix22ButterflyG1_NF_add_cast_2 = {dinXTwdl2_im[16], dinXTwdl2_im};
      Radix22ButterflyG1_NF_btf1_im_reg_next = Radix22ButterflyG1_NF_add_cast_1 + Radix22ButterflyG1_NF_add_cast_2;
      Radix22ButterflyG1_NF_sub_cast_1 = {din1_im_dly9[16], din1_im_dly9};
      Radix22ButterflyG1_NF_sub_cast_2 = {dinXTwdl2_im[16], dinXTwdl2_im};
      Radix22ButterflyG1_NF_btf2_im_reg_next = Radix22ButterflyG1_NF_sub_cast_1 - Radix22ButterflyG1_NF_sub_cast_2;
    end
    dout_25_re_2 = Radix22ButterflyG1_NF_btf1_re_reg[16:0];
    dout_25_im_2 = Radix22ButterflyG1_NF_btf1_im_reg[16:0];
    dout_26_re_1 = Radix22ButterflyG1_NF_btf2_re_reg[16:0];
    dout_26_im_1 = Radix22ButterflyG1_NF_btf2_im_reg[16:0];
    dout_26_vld = Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  end



  assign dout_25_re_1 = dout_25_re_2;

  assign dout_25_im_1 = dout_25_im_2;

  assign dout_26_re = dout_26_re_1;

  assign dout_26_im = dout_26_im_1;

endmodule  // RADIX22FFT_SDNF1_5_block11

