// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\Complex4Multiply_block.v
// Created: 2023-09-29 06:09:50
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: Complex4Multiply_block
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/TWDLMULT_SDNF1_3/Complex4Multiply
// Hierarchy Level: 6
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module Complex4Multiply_block
          (clk,
           reset,
           enb,
           din1_re_dly3,
           din1_im_dly3,
           din1_vld_dly3,
           twdl_3_3_re,
           twdl_3_3_im,
           twdlXdin_3_re,
           twdlXdin_3_im);


  input   clk;
  input   reset;
  input   enb;
  input   signed [14:0] din1_re_dly3;  // sfix15_En10
  input   signed [14:0] din1_im_dly3;  // sfix15_En10
  input   din1_vld_dly3;
  input   signed [11:0] twdl_3_3_re;  // sfix12_En10
  input   signed [11:0] twdl_3_3_im;  // sfix12_En10
  output  signed [14:0] twdlXdin_3_re;  // sfix15_En10
  output  signed [14:0] twdlXdin_3_im;  // sfix15_En10


  reg signed [14:0] din_re_reg;  // sfix15_En10
  reg signed [14:0] din_im_reg;  // sfix15_En10
  reg signed [11:0] twdl_re_reg;  // sfix12_En10
  reg signed [11:0] twdl_im_reg;  // sfix12_En10
  reg signed [14:0] Complex4Multiply_din1_re_pipe1;  // sfix15
  reg signed [14:0] Complex4Multiply_din1_im_pipe1;  // sfix15
  reg signed [26:0] Complex4Multiply_mult1_re_pipe1;  // sfix27
  reg signed [26:0] Complex4Multiply_mult2_re_pipe1;  // sfix27
  reg signed [26:0] Complex4Multiply_mult1_im_pipe1;  // sfix27
  reg signed [26:0] Complex4Multiply_mult2_im_pipe1;  // sfix27
  reg signed [11:0] Complex4Multiply_twiddle_re_pipe1;  // sfix12
  reg signed [11:0] Complex4Multiply_twiddle_im_pipe1;  // sfix12
  reg signed [26:0] prod1_re;  // sfix27_En20
  reg signed [26:0] prod1_im;  // sfix27_En20
  reg signed [26:0] prod2_re;  // sfix27_En20
  reg signed [26:0] prod2_im;  // sfix27_En20
  reg  din_vld_dly1;
  reg  din_vld_dly2;
  reg  din_vld_dly3;
  reg  prod_vld;
  reg signed [27:0] Complex4Add_multRes_re_reg;  // sfix28
  reg signed [27:0] Complex4Add_multRes_im_reg;  // sfix28
  reg  Complex4Add_prod_vld_reg1;
  reg signed [26:0] Complex4Add_prod1_re_reg;  // sfix27
  reg signed [26:0] Complex4Add_prod1_im_reg;  // sfix27
  reg signed [26:0] Complex4Add_prod2_re_reg;  // sfix27
  reg signed [26:0] Complex4Add_prod2_im_reg;  // sfix27
  wire signed [27:0] Complex4Add_multRes_re_reg_next;  // sfix28_En20
  wire signed [27:0] Complex4Add_multRes_im_reg_next;  // sfix28_En20
  wire signed [27:0] Complex4Add_sub_cast;  // sfix28_En20
  wire signed [27:0] Complex4Add_sub_cast_1;  // sfix28_En20
  wire signed [27:0] Complex4Add_add_cast;  // sfix28_En20
  wire signed [27:0] Complex4Add_add_cast_1;  // sfix28_En20
  wire signed [27:0] mulResFP_re;  // sfix28_En20
  wire signed [27:0] mulResFP_im;  // sfix28_En20
  reg  twdlXdin1_vld;

  initial begin
    Complex4Multiply_din1_re_pipe1 = 15'sb000000000000000;
    Complex4Multiply_din1_im_pipe1 = 15'sb000000000000000;
    Complex4Multiply_twiddle_re_pipe1 = 12'sb000000000000;
    Complex4Multiply_twiddle_im_pipe1 = 12'sb000000000000;
    Complex4Multiply_mult1_re_pipe1 = 27'sb000000000000000000000000000;
    Complex4Multiply_mult2_re_pipe1 = 27'sb000000000000000000000000000;
    Complex4Multiply_mult1_im_pipe1 = 27'sb000000000000000000000000000;
    Complex4Multiply_mult2_im_pipe1 = 27'sb000000000000000000000000000;
    prod1_re = 27'sb000000000000000000000000000;
    prod2_re = 27'sb000000000000000000000000000;
    prod1_im = 27'sb000000000000000000000000000;
    prod2_im = 27'sb000000000000000000000000000;
  end

  always @(posedge clk or posedge reset)
    begin : intdelay_process
      if (reset == 1'b1) begin
        din_re_reg <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din_re_reg <= din1_re_dly3;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_1_process
      if (reset == 1'b1) begin
        din_im_reg <= 15'sb000000000000000;
      end
      else begin
        if (enb) begin
          din_im_reg <= din1_im_dly3;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_2_process
      if (reset == 1'b1) begin
        twdl_re_reg <= 12'sb000000000000;
      end
      else begin
        if (enb) begin
          twdl_re_reg <= twdl_3_3_re;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_3_process
      if (reset == 1'b1) begin
        twdl_im_reg <= 12'sb000000000000;
      end
      else begin
        if (enb) begin
          twdl_im_reg <= twdl_3_3_im;
        end
      end
    end



  // Complex4Multiply
  always @(posedge clk)
    begin : Complex4Multiply_process
      if (enb) begin
        prod1_re <= Complex4Multiply_mult1_re_pipe1;
        prod2_re <= Complex4Multiply_mult2_re_pipe1;
        prod1_im <= Complex4Multiply_mult1_im_pipe1;
        prod2_im <= Complex4Multiply_mult2_im_pipe1;
        Complex4Multiply_mult1_re_pipe1 <= Complex4Multiply_din1_re_pipe1 * Complex4Multiply_twiddle_re_pipe1;
        Complex4Multiply_mult2_re_pipe1 <= Complex4Multiply_din1_im_pipe1 * Complex4Multiply_twiddle_im_pipe1;
        Complex4Multiply_mult1_im_pipe1 <= Complex4Multiply_din1_re_pipe1 * Complex4Multiply_twiddle_im_pipe1;
        Complex4Multiply_mult2_im_pipe1 <= Complex4Multiply_din1_im_pipe1 * Complex4Multiply_twiddle_re_pipe1;
        Complex4Multiply_twiddle_re_pipe1 <= twdl_re_reg;
        Complex4Multiply_twiddle_im_pipe1 <= twdl_im_reg;
        Complex4Multiply_din1_re_pipe1 <= din_re_reg;
        Complex4Multiply_din1_im_pipe1 <= din_im_reg;
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_4_process
      if (reset == 1'b1) begin
        din_vld_dly1 <= 1'b0;
      end
      else begin
        if (enb) begin
          din_vld_dly1 <= din1_vld_dly3;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_5_process
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
    begin : intdelay_6_process
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
    begin : intdelay_7_process
      if (reset == 1'b1) begin
        prod_vld <= 1'b0;
      end
      else begin
        if (enb) begin
          prod_vld <= din_vld_dly3;
        end
      end
    end



  // Complex4Add
  always @(posedge clk or posedge reset)
    begin : Complex4Add_process
      if (reset == 1'b1) begin
        Complex4Add_multRes_re_reg <= 28'sb0000000000000000000000000000;
        Complex4Add_multRes_im_reg <= 28'sb0000000000000000000000000000;
        Complex4Add_prod1_re_reg <= 27'sb000000000000000000000000000;
        Complex4Add_prod1_im_reg <= 27'sb000000000000000000000000000;
        Complex4Add_prod2_re_reg <= 27'sb000000000000000000000000000;
        Complex4Add_prod2_im_reg <= 27'sb000000000000000000000000000;
        Complex4Add_prod_vld_reg1 <= 1'b0;
        twdlXdin1_vld <= 1'b0;
      end
      else begin
        if (enb) begin
          Complex4Add_multRes_re_reg <= Complex4Add_multRes_re_reg_next;
          Complex4Add_multRes_im_reg <= Complex4Add_multRes_im_reg_next;
          Complex4Add_prod1_re_reg <= prod1_re;
          Complex4Add_prod1_im_reg <= prod1_im;
          Complex4Add_prod2_re_reg <= prod2_re;
          Complex4Add_prod2_im_reg <= prod2_im;
          twdlXdin1_vld <= Complex4Add_prod_vld_reg1;
          Complex4Add_prod_vld_reg1 <= prod_vld;
        end
      end
    end

  assign Complex4Add_sub_cast = {Complex4Add_prod1_re_reg[26], Complex4Add_prod1_re_reg};
  assign Complex4Add_sub_cast_1 = {Complex4Add_prod2_re_reg[26], Complex4Add_prod2_re_reg};
  assign Complex4Add_multRes_re_reg_next = Complex4Add_sub_cast - Complex4Add_sub_cast_1;
  assign Complex4Add_add_cast = {Complex4Add_prod1_im_reg[26], Complex4Add_prod1_im_reg};
  assign Complex4Add_add_cast_1 = {Complex4Add_prod2_im_reg[26], Complex4Add_prod2_im_reg};
  assign Complex4Add_multRes_im_reg_next = Complex4Add_add_cast + Complex4Add_add_cast_1;
  assign mulResFP_re = Complex4Add_multRes_re_reg;
  assign mulResFP_im = Complex4Add_multRes_im_reg;



  assign twdlXdin_3_re = mulResFP_re[24:10];



  assign twdlXdin_3_im = mulResFP_im[24:10];



endmodule  // Complex4Multiply_block

