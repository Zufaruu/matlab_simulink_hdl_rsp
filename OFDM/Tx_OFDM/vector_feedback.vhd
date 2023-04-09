-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\vector_feedback.vhd
-- Created: 2023-04-09 12:21:03
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: vector_feedback
-- Source Path: Tx_OFDM/Transmitter/custom buffer/vector feedback
-- Hierarchy Level: 2
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;
USE work.Transmitter_pkg.ALL;

ENTITY vector_feedback IS
  PORT( vector_before_re                  :   IN    vector_of_std_logic_vector16(0 TO 63);  -- sfix16_En14 [64]
        vector_before_im                  :   IN    vector_of_std_logic_vector16(0 TO 63);  -- sfix16_En14 [64]
        y_re                              :   OUT   vector_of_std_logic_vector16(0 TO 63);  -- sfix16_En14 [64]
        y_im                              :   OUT   vector_of_std_logic_vector16(0 TO 63)  -- sfix16_En14 [64]
        );
END vector_feedback;


ARCHITECTURE rtl OF vector_feedback IS

  -- Signals
  SIGNAL vector_before_re_signed          : vector_of_signed16(0 TO 63);  -- sfix16_En14 [64]
  SIGNAL vector_before_im_signed          : vector_of_signed16(0 TO 63);  -- sfix16_En14 [64]
  SIGNAL y_re_tmp                         : vector_of_signed16(0 TO 63);  -- sfix16_En14 [64]
  SIGNAL y_im_tmp                         : vector_of_signed16(0 TO 63);  -- sfix16_En14 [64]

BEGIN
  outputgen3: FOR k IN 0 TO 63 GENERATE
    vector_before_re_signed(k) <= signed(vector_before_re(k));
  END GENERATE;

  outputgen2: FOR k IN 0 TO 63 GENERATE
    vector_before_im_signed(k) <= signed(vector_before_im(k));
  END GENERATE;

  y_re_tmp <= vector_before_re_signed;
  y_im_tmp <= vector_before_im_signed;

  outputgen1: FOR k IN 0 TO 63 GENERATE
    y_re(k) <= std_logic_vector(y_re_tmp(k));
  END GENERATE;

  outputgen: FOR k IN 0 TO 63 GENERATE
    y_im(k) <= std_logic_vector(y_im_tmp(k));
  END GENERATE;

END rtl;

