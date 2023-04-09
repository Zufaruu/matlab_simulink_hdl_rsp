-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\RADIX22FFT_SDNF1_7_block14.vhd
-- Created: 2023-04-09 19:16:04
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: RADIX22FFT_SDNF1_7_block14
-- Source Path: dsphdl.IFFT/RADIX22FFT_SDNF1_7
-- Hierarchy Level: 3
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY RADIX22FFT_SDNF1_7_block14 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_960_0                       :   IN    std_logic;
        dout_30_re                        :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_30_im                        :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_32_re                        :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_32_im                        :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_6_vld                        :   IN    std_logic;
        twdl_7_32_re                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdl_7_32_im                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_31_re                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_31_im                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_32_re_1                      :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_32_im_1                      :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END RADIX22FFT_SDNF1_7_block14;


ARCHITECTURE rtl OF RADIX22FFT_SDNF1_7_block14 IS

  -- Component Declarations
  COMPONENT Complex4Multiply_block132
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb_1_960_0                     :   IN    std_logic;
          din2_re_dly3                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          din2_im_dly3                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          di2_vld_dly3                    :   IN    std_logic;
          twdl_7_32_re                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_7_32_im                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          dinXTwdl2_re                    :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          dinXTwdl2_im                    :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
          );
  END COMPONENT;

  -- Component Configuration Statements
  FOR ALL : Complex4Multiply_block132
    USE ENTITY work.Complex4Multiply_block132(rtl);

  -- Signals
  SIGNAL dout_30_re_signed                : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly4                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly5                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly6                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly7                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly8                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly9                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_30_im_signed                : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly4                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly5                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly6                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly7                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly8                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly9                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_32_re_signed                : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_re_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_re_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_32_im_signed                : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_im_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_im_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_re_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_im_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL di2_vld_dly1                     : std_logic;
  SIGNAL di2_vld_dly2                     : std_logic;
  SIGNAL di2_vld_dly3                     : std_logic;
  SIGNAL dinXTwdl2_re                     : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dinXTwdl2_im                     : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dinXTwdl2_re_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dinXTwdl2_im_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din_vld_dly1                     : std_logic;
  SIGNAL din_vld_dly2                     : std_logic;
  SIGNAL din_vld_dly3                     : std_logic;
  SIGNAL din_vld_dly4                     : std_logic;
  SIGNAL din_vld_dly5                     : std_logic;
  SIGNAL din_vld_dly6                     : std_logic;
  SIGNAL din_vld_dly7                     : std_logic;
  SIGNAL din_vld_dly8                     : std_logic;
  SIGNAL din_vld_dly9                     : std_logic;
  SIGNAL Radix22ButterflyG1_NF_btf1_re_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_NF_btf1_im_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_NF_btf2_re_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_NF_btf2_im_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_NF_dinXtwdl_vld_dly1 : std_logic;
  SIGNAL Radix22ButterflyG1_NF_btf1_re_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_NF_btf1_im_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_NF_btf2_re_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_NF_btf2_im_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next : std_logic;
  SIGNAL dout_31_re_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_31_im_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_32_re_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_32_im_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_32_vld                      : std_logic;

BEGIN
  u_MUL4_2 : Complex4Multiply_block132
    PORT MAP( clk => clk,
              reset => reset,
              enb_1_960_0 => enb_1_960_0,
              din2_re_dly3 => std_logic_vector(din2_re_dly3),  -- sfix16_En14
              din2_im_dly3 => std_logic_vector(din2_im_dly3),  -- sfix16_En14
              di2_vld_dly3 => di2_vld_dly3,
              twdl_7_32_re => twdl_7_32_re,  -- sfix16_En14
              twdl_7_32_im => twdl_7_32_im,  -- sfix16_En14
              dinXTwdl2_re => dinXTwdl2_re,  -- sfix16_En14
              dinXTwdl2_im => dinXTwdl2_im  -- sfix16_En14
              );

  dout_30_re_signed <= signed(dout_30_re);

  intdelay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly1 <= dout_30_re_signed;
      END IF;
    END IF;
  END PROCESS intdelay_process;


  intdelay_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly2 <= din1_re_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_1_process;


  intdelay_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly3 <= din1_re_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_2_process;


  intdelay_3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly4 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly4 <= din1_re_dly3;
      END IF;
    END IF;
  END PROCESS intdelay_3_process;


  intdelay_4_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly5 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly5 <= din1_re_dly4;
      END IF;
    END IF;
  END PROCESS intdelay_4_process;


  intdelay_5_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly6 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly6 <= din1_re_dly5;
      END IF;
    END IF;
  END PROCESS intdelay_5_process;


  intdelay_6_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly7 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly7 <= din1_re_dly6;
      END IF;
    END IF;
  END PROCESS intdelay_6_process;


  intdelay_7_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly8 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly8 <= din1_re_dly7;
      END IF;
    END IF;
  END PROCESS intdelay_7_process;


  intdelay_8_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly9 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly9 <= din1_re_dly8;
      END IF;
    END IF;
  END PROCESS intdelay_8_process;


  dout_30_im_signed <= signed(dout_30_im);

  intdelay_9_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly1 <= dout_30_im_signed;
      END IF;
    END IF;
  END PROCESS intdelay_9_process;


  intdelay_10_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly2 <= din1_im_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_10_process;


  intdelay_11_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly3 <= din1_im_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_11_process;


  intdelay_12_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly4 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly4 <= din1_im_dly3;
      END IF;
    END IF;
  END PROCESS intdelay_12_process;


  intdelay_13_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly5 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly5 <= din1_im_dly4;
      END IF;
    END IF;
  END PROCESS intdelay_13_process;


  intdelay_14_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly6 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly6 <= din1_im_dly5;
      END IF;
    END IF;
  END PROCESS intdelay_14_process;


  intdelay_15_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly7 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly7 <= din1_im_dly6;
      END IF;
    END IF;
  END PROCESS intdelay_15_process;


  intdelay_16_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly8 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly8 <= din1_im_dly7;
      END IF;
    END IF;
  END PROCESS intdelay_16_process;


  intdelay_17_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly9 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly9 <= din1_im_dly8;
      END IF;
    END IF;
  END PROCESS intdelay_17_process;


  dout_32_re_signed <= signed(dout_32_re);

  intdelay_18_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_re_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_re_dly1 <= dout_32_re_signed;
      END IF;
    END IF;
  END PROCESS intdelay_18_process;


  intdelay_19_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_re_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_re_dly2 <= din2_re_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_19_process;


  dout_32_im_signed <= signed(dout_32_im);

  intdelay_20_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_im_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_im_dly1 <= dout_32_im_signed;
      END IF;
    END IF;
  END PROCESS intdelay_20_process;


  intdelay_21_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_im_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_im_dly2 <= din2_im_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_21_process;


  intdelay_22_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_re_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_re_dly3 <= din2_re_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_22_process;


  intdelay_23_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_im_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_im_dly3 <= din2_im_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_23_process;


  intdelay_24_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      di2_vld_dly1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        di2_vld_dly1 <= dout_6_vld;
      END IF;
    END IF;
  END PROCESS intdelay_24_process;


  intdelay_25_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      di2_vld_dly2 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        di2_vld_dly2 <= di2_vld_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_25_process;


  intdelay_26_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      di2_vld_dly3 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        di2_vld_dly3 <= di2_vld_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_26_process;


  dinXTwdl2_re_signed <= signed(dinXTwdl2_re);

  dinXTwdl2_im_signed <= signed(dinXTwdl2_im);

  intdelay_27_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly1 <= dout_6_vld;
      END IF;
    END IF;
  END PROCESS intdelay_27_process;


  intdelay_28_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly2 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly2 <= din_vld_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_28_process;


  intdelay_29_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly3 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly3 <= din_vld_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_29_process;


  intdelay_30_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly4 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly4 <= din_vld_dly3;
      END IF;
    END IF;
  END PROCESS intdelay_30_process;


  intdelay_31_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly5 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly5 <= din_vld_dly4;
      END IF;
    END IF;
  END PROCESS intdelay_31_process;


  intdelay_32_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly6 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly6 <= din_vld_dly5;
      END IF;
    END IF;
  END PROCESS intdelay_32_process;


  intdelay_33_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly7 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly7 <= din_vld_dly6;
      END IF;
    END IF;
  END PROCESS intdelay_33_process;


  intdelay_34_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly8 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly8 <= din_vld_dly7;
      END IF;
    END IF;
  END PROCESS intdelay_34_process;


  intdelay_35_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly9 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly9 <= din_vld_dly8;
      END IF;
    END IF;
  END PROCESS intdelay_35_process;


  -- Radix22ButterflyG1_NF
  Radix22ButterflyG1_NF_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Radix22ButterflyG1_NF_btf1_re_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_NF_btf1_im_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_NF_btf2_re_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_NF_btf2_im_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_NF_dinXtwdl_vld_dly1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        Radix22ButterflyG1_NF_btf1_re_reg <= Radix22ButterflyG1_NF_btf1_re_reg_next;
        Radix22ButterflyG1_NF_btf1_im_reg <= Radix22ButterflyG1_NF_btf1_im_reg_next;
        Radix22ButterflyG1_NF_btf2_re_reg <= Radix22ButterflyG1_NF_btf2_re_reg_next;
        Radix22ButterflyG1_NF_btf2_im_reg <= Radix22ButterflyG1_NF_btf2_im_reg_next;
        Radix22ButterflyG1_NF_dinXtwdl_vld_dly1 <= Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next;
      END IF;
    END IF;
  END PROCESS Radix22ButterflyG1_NF_process;

  Radix22ButterflyG1_NF_output : PROCESS (Radix22ButterflyG1_NF_btf1_im_reg, Radix22ButterflyG1_NF_btf1_re_reg,
       Radix22ButterflyG1_NF_btf2_im_reg, Radix22ButterflyG1_NF_btf2_re_reg,
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, din1_im_dly9, din1_re_dly9,
       dinXTwdl2_im_signed, dinXTwdl2_re_signed, din_vld_dly9)
    VARIABLE add_cast : signed(16 DOWNTO 0);
    VARIABLE add_cast_0 : signed(16 DOWNTO 0);
    VARIABLE sra_temp : signed(16 DOWNTO 0);
    VARIABLE sub_cast : signed(16 DOWNTO 0);
    VARIABLE sub_cast_0 : signed(16 DOWNTO 0);
    VARIABLE sra_temp_0 : signed(16 DOWNTO 0);
    VARIABLE add_cast_1 : signed(16 DOWNTO 0);
    VARIABLE add_cast_2 : signed(16 DOWNTO 0);
    VARIABLE sra_temp_1 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_1 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_2 : signed(16 DOWNTO 0);
    VARIABLE sra_temp_2 : signed(16 DOWNTO 0);
  BEGIN
    add_cast := to_signed(16#00000#, 17);
    add_cast_0 := to_signed(16#00000#, 17);
    sub_cast := to_signed(16#00000#, 17);
    sub_cast_0 := to_signed(16#00000#, 17);
    add_cast_1 := to_signed(16#00000#, 17);
    add_cast_2 := to_signed(16#00000#, 17);
    sub_cast_1 := to_signed(16#00000#, 17);
    sub_cast_2 := to_signed(16#00000#, 17);
    Radix22ButterflyG1_NF_btf1_re_reg_next <= Radix22ButterflyG1_NF_btf1_re_reg;
    Radix22ButterflyG1_NF_btf1_im_reg_next <= Radix22ButterflyG1_NF_btf1_im_reg;
    Radix22ButterflyG1_NF_btf2_re_reg_next <= Radix22ButterflyG1_NF_btf2_re_reg;
    Radix22ButterflyG1_NF_btf2_im_reg_next <= Radix22ButterflyG1_NF_btf2_im_reg;
    Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next <= din_vld_dly9;
    IF din_vld_dly9 = '1' THEN 
      add_cast := resize(din1_re_dly9, 17);
      add_cast_0 := resize(dinXTwdl2_re_signed, 17);
      Radix22ButterflyG1_NF_btf1_re_reg_next <= add_cast + add_cast_0;
      sub_cast := resize(din1_re_dly9, 17);
      sub_cast_0 := resize(dinXTwdl2_re_signed, 17);
      Radix22ButterflyG1_NF_btf2_re_reg_next <= sub_cast - sub_cast_0;
      add_cast_1 := resize(din1_im_dly9, 17);
      add_cast_2 := resize(dinXTwdl2_im_signed, 17);
      Radix22ButterflyG1_NF_btf1_im_reg_next <= add_cast_1 + add_cast_2;
      sub_cast_1 := resize(din1_im_dly9, 17);
      sub_cast_2 := resize(dinXTwdl2_im_signed, 17);
      Radix22ButterflyG1_NF_btf2_im_reg_next <= sub_cast_1 - sub_cast_2;
    END IF;
    sra_temp := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf1_re_reg, 1);
    dout_31_re_tmp <= sra_temp(15 DOWNTO 0);
    sra_temp_0 := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf1_im_reg, 1);
    dout_31_im_tmp <= sra_temp_0(15 DOWNTO 0);
    sra_temp_1 := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf2_re_reg, 1);
    dout_32_re_tmp <= sra_temp_1(15 DOWNTO 0);
    sra_temp_2 := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf2_im_reg, 1);
    dout_32_im_tmp <= sra_temp_2(15 DOWNTO 0);
    dout_32_vld <= Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  END PROCESS Radix22ButterflyG1_NF_output;


  dout_31_re <= std_logic_vector(dout_31_re_tmp);

  dout_31_im <= std_logic_vector(dout_31_im_tmp);

  dout_32_re_1 <= std_logic_vector(dout_32_re_tmp);

  dout_32_im_1 <= std_logic_vector(dout_32_im_tmp);

END rtl;

