// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\Vel_Calculation.v
// Created: 2023-09-27 04:10:36
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: Vel_Calculation
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Vel Calculation
// Hierarchy Level: 3
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module Vel_Calculation
          (doppler_freq,
           vel);


  input   [31:0] doppler_freq;  // ufix32_En5
  output  [43:0] vel;  // ufix44_En11


  wire [11:0] Lambda_out1;  // ufix12_En6
  wire [43:0] Product_out1;  // ufix44_En11
  wire [7:0] RTT_Constant_out1;  // uint8
  reg [43:0] Divide_out1;  // ufix44_En11
  reg [44:0] Divide_div_temp;  // ufix45_En11
  reg [44:0] Divide_cast;  // ufix45_En11


  assign Lambda_out1 = 12'b000000000110;



  assign Product_out1 = Lambda_out1 * doppler_freq;



  assign RTT_Constant_out1 = 8'b00000010;



  always @(Product_out1, RTT_Constant_out1) begin
    Divide_div_temp = 45'h000000000000;
    Divide_cast = 45'h000000000000;
    if (RTT_Constant_out1 == 8'b00000000) begin
      Divide_out1 = 44'hFFFFFFFFFFF;
    end
    else begin
      Divide_cast = {1'b0, Product_out1};
      Divide_div_temp = Divide_cast / RTT_Constant_out1;
      if (Divide_div_temp[44] != 1'b0) begin
        Divide_out1 = 44'hFFFFFFFFFFF;
      end
      else begin
        Divide_out1 = Divide_div_temp[43:0];
      end
    end
  end



  assign vel = Divide_out1;

endmodule  // Vel_Calculation
