// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\TWDLMULT_SDNF1_3_block8.v
// Created: 2023-09-29 06:09:50
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: TWDLMULT_SDNF1_3_block8
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/TWDLMULT_SDNF1_3
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module TWDLMULT_SDNF1_3_block8
          (clk,
           reset,
           enb,
           dout_21_re,
           dout_21_im,
           dout_23_re,
           dout_23_im,
           dout_2_vld,
           twdl_3_19_re,
           twdl_3_19_im,
           twdl_3_20_re,
           twdl_3_20_im,
           twdlXdin_19_re,
           twdlXdin_19_im,
           twdlXdin_20_re,
           twdlXdin_20_im);


  input   clk;
  input   reset;
  input   enb;
  input   signed [13:0] dout_21_re;  // sfix14_En10
  input   signed [13:0] dout_21_im;  // sfix14_En10
  input   signed [13:0] dout_23_re;  // sfix14_En10
  input   signed [13:0] dout_23_im;  // sfix14_En10
  input   dout_2_vld;
  input   signed [11:0] twdl_3_19_re;  // sfix12_En10
  input   signed [11:0] twdl_3_19_im;  // sfix12_En10
  input   signed [11:0] twdl_3_20_re;  // sfix12_En10
  input   signed [11:0] twdl_3_20_im;  // sfix12_En10
  output  signed [14:0] twdlXdin_19_re;  // sfix15_En10
  output  signed [14:0] twdlXdin_19_im;  // sfix15_En10
  output  signed [14:0] twdlXdin_20_re;  // sfix15_En10
  output  signed [14:0] twdlXdin_20_im;  // sfix15_En10


  wire signed [14:0] din_re;  // sfix15_En10
  reg signed [14:0] din1_re_dly1;  // sfix15_En10
  reg signed [14:0] din1_re_dly2;  // sfix15_En10
  reg signed [14:0] din1_re_dly3;  // sfix15_En10
  wire signed [14:0] din_im;  // sfix15_En10
  reg signed [14:0] din1_im_dly1;  // sfix15_En10
  reg signed [14:0] din1_im_dly2;  // sfix15_En10
  reg signed [14:0] din1_im_dly3;  // sfix15_En10
  reg  din1_vld_dly1;
  reg  din1_vld_dly2;
  reg  din1_vld_dly3;
  wire signed [14:0] din_re_1;  // sfix15_En10
  reg signed [14:0] din2_re_dly1;  // sfix15_En10
  reg signed [14:0] din2_re_dly2;  // sfix15_En10
  reg signed [14:0] din2_re_dly3;  // sfix15_En10
  wire signed [14:0] din_im_1;  // sfix15_En10
  reg signed [14:0] din2_im_dly1;  // sfix15_En10
  reg signed [14:0] din2_im_dly2;  // sfix15_En10
  reg signed [14:0] din2_im_dly3;  // sfix15_En10
  reg  di2_vld_dly1;
  reg  di2_vld_dly2;
  reg  di2_vld_dly3;


  assign din_re = {dout_21_re[13], dout_21_re};



  always @(posedge clk or posedge reset)
    begin : intdelay_process
      if (reset == 1'b1) begin
        din1_re_dly1 <= 15'sb000000000000000;
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
        din1_re_dly2 <= 15'sb000000000000000;
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
        din1_re_dly3 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din1_re_dly3 <= din1_re_dly2;
        end
      end
    end



  assign din_im = {dout_21_im[13], dout_21_im};



  always @(posedge clk or posedge reset)
    begin : intdelay_3_process
      if (reset == 1'b1) begin
        din1_im_dly1 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly1 <= din_im;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_4_process
      if (reset == 1'b1) begin
        din1_im_dly2 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly2 <= din1_im_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_5_process
      if (reset == 1'b1) begin
        din1_im_dly3 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din1_im_dly3 <= din1_im_dly2;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_6_process
      if (reset == 1'b1) begin
        din1_vld_dly1 <= 1'b0;
      end
      else begin
        if (enb) begin
          din1_vld_dly1 <= dout_2_vld;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_7_process
      if (reset == 1'b1) begin
        din1_vld_dly2 <= 1'b0;
      end
      else begin
        if (enb) begin
          din1_vld_dly2 <= din1_vld_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_8_process
      if (reset == 1'b1) begin
        din1_vld_dly3 <= 1'b0;
      end
      else begin
        if (enb) begin
          din1_vld_dly3 <= din1_vld_dly2;
        end
      end
    end



  Complex4Multiply_block14 u_MUL4_1 (.clk(clk),
                                     .reset(reset),
                                     .enb(enb),
                                     .din1_re_dly3(din1_re_dly3),  // sfix15_En10
                                     .din1_im_dly3(din1_im_dly3),  // sfix15_En10
                                     .din1_vld_dly3(din1_vld_dly3),
                                     .twdl_3_19_re(twdl_3_19_re),  // sfix12_En10
                                     .twdl_3_19_im(twdl_3_19_im),  // sfix12_En10
                                     .twdlXdin_19_re(twdlXdin_19_re),  // sfix15_En10
                                     .twdlXdin_19_im(twdlXdin_19_im)  // sfix15_En10
                                     );

  assign din_re_1 = {dout_23_re[13], dout_23_re};



  always @(posedge clk or posedge reset)
    begin : intdelay_9_process
      if (reset == 1'b1) begin
        din2_re_dly1 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din2_re_dly1 <= din_re_1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_10_process
      if (reset == 1'b1) begin
        din2_re_dly2 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din2_re_dly2 <= din2_re_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_11_process
      if (reset == 1'b1) begin
        din2_re_dly3 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din2_re_dly3 <= din2_re_dly2;
        end
      end
    end



  assign din_im_1 = {dout_23_im[13], dout_23_im};



  always @(posedge clk or posedge reset)
    begin : intdelay_12_process
      if (reset == 1'b1) begin
        din2_im_dly1 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din2_im_dly1 <= din_im_1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_13_process
      if (reset == 1'b1) begin
        din2_im_dly2 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din2_im_dly2 <= din2_im_dly1;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_14_process
      if (reset == 1'b1) begin
        din2_im_dly3 <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din2_im_dly3 <= din2_im_dly2;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_15_process
      if (reset == 1'b1) begin
        di2_vld_dly1 <= 1'b0;
      end
      else begin
        if (enb) begin
          di2_vld_dly1 <= dout_2_vld;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_16_process
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
    begin : intdelay_17_process
      if (reset == 1'b1) begin
        di2_vld_dly3 <= 1'b0;
      end
      else begin
        if (enb) begin
          di2_vld_dly3 <= di2_vld_dly2;
        end
      end
    end



  Complex4Multiply_block15 u_MUL4_2 (.clk(clk),
                                     .reset(reset),
                                     .enb(enb),
                                     .din2_re_dly3(din2_re_dly3),  // sfix15_En10
                                     .din2_im_dly3(din2_im_dly3),  // sfix15_En10
                                     .di2_vld_dly3(di2_vld_dly3),
                                     .twdl_3_20_re(twdl_3_20_re),  // sfix12_En10
                                     .twdl_3_20_im(twdl_3_20_im),  // sfix12_En10
                                     .twdlXdin_20_re(twdlXdin_20_re),  // sfix15_En10
                                     .twdlXdin_20_im(twdlXdin_20_im)  // sfix15_En10
                                     );

endmodule  // TWDLMULT_SDNF1_3_block8

