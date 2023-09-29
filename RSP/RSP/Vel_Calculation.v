// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\Vel_Calculation.v
// Created: 2023-09-29 06:09:51
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


  assign Lambda_out1 = 12'b000000000110;



  assign Product_out1 = Lambda_out1 * doppler_freq;



  assign vel = Product_out1;

endmodule  // Vel_Calculation

