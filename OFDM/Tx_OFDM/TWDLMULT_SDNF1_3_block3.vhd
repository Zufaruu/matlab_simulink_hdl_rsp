-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\TWDLMULT_SDNF1_3_block3.vhd
-- Created: 2023-04-09 15:00:11
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: TWDLMULT_SDNF1_3_block3
-- Source Path: Tx_OFDM/Transmitter/IFFT/TWDLMULT_SDNF1_3
-- Hierarchy Level: 2
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY TWDLMULT_SDNF1_3_block3 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_960_0                       :   IN    std_logic;
        dout_2_9_re                       :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_2_9_im                       :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_2_10_re                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_2_10_im                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_2_1_vld                      :   IN    std_logic;
        twdl_3_9_re                       :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdl_3_9_im                       :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdl_3_10_re                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdl_3_10_im                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_9_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_9_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_10_re                    :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_10_im                    :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END TWDLMULT_SDNF1_3_block3;


ARCHITECTURE rtl OF TWDLMULT_SDNF1_3_block3 IS

  -- Component Declarations
  COMPONENT Complex4Multiply_block6
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb_1_960_0                     :   IN    std_logic;
          din1_re_dly3                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          din1_im_dly3                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          din1_vld_dly3                   :   IN    std_logic;
          twdl_3_9_re                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_3_9_im                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdlXdin_9_re                   :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdlXdin_9_im                   :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
          );
  END COMPONENT;

  COMPONENT Complex4Multiply_block7
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb_1_960_0                     :   IN    std_logic;
          din2_re_dly3                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          din2_im_dly3                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          di2_vld_dly3                    :   IN    std_logic;
          twdl_3_10_re                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_3_10_im                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdlXdin_10_re                  :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdlXdin_10_im                  :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
          );
  END COMPONENT;

  -- Component Configuration Statements
  FOR ALL : Complex4Multiply_block6
    USE ENTITY work.Complex4Multiply_block6(rtl);

  FOR ALL : Complex4Multiply_block7
    USE ENTITY work.Complex4Multiply_block7(rtl);

  -- Signals
  SIGNAL dout_2_9_re_signed               : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_2_9_im_signed               : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_re_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_vld_dly1                    : std_logic;
  SIGNAL din1_vld_dly2                    : std_logic;
  SIGNAL din1_vld_dly3                    : std_logic;
  SIGNAL twdlXdin_9_re_tmp                : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL twdlXdin_9_im_tmp                : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_2_10_re_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_re_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_re_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_2_10_im_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_im_dly1                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_im_dly2                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_re_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din2_im_dly3                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL di2_vld_dly1                     : std_logic;
  SIGNAL di2_vld_dly2                     : std_logic;
  SIGNAL di2_vld_dly3                     : std_logic;
  SIGNAL twdlXdin_10_re_tmp               : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL twdlXdin_10_im_tmp               : std_logic_vector(15 DOWNTO 0);  -- ufix16

BEGIN
  u_MUL4_1 : Complex4Multiply_block6
    PORT MAP( clk => clk,
              reset => reset,
              enb_1_960_0 => enb_1_960_0,
              din1_re_dly3 => std_logic_vector(din1_re_dly3),  -- sfix16_En14
              din1_im_dly3 => std_logic_vector(din1_im_dly3),  -- sfix16_En14
              din1_vld_dly3 => din1_vld_dly3,
              twdl_3_9_re => twdl_3_9_re,  -- sfix16_En14
              twdl_3_9_im => twdl_3_9_im,  -- sfix16_En14
              twdlXdin_9_re => twdlXdin_9_re_tmp,  -- sfix16_En14
              twdlXdin_9_im => twdlXdin_9_im_tmp  -- sfix16_En14
              );

  u_MUL4_2 : Complex4Multiply_block7
    PORT MAP( clk => clk,
              reset => reset,
              enb_1_960_0 => enb_1_960_0,
              din2_re_dly3 => std_logic_vector(din2_re_dly3),  -- sfix16_En14
              din2_im_dly3 => std_logic_vector(din2_im_dly3),  -- sfix16_En14
              di2_vld_dly3 => di2_vld_dly3,
              twdl_3_10_re => twdl_3_10_re,  -- sfix16_En14
              twdl_3_10_im => twdl_3_10_im,  -- sfix16_En14
              twdlXdin_10_re => twdlXdin_10_re_tmp,  -- sfix16_En14
              twdlXdin_10_im => twdlXdin_10_im_tmp  -- sfix16_En14
              );

  dout_2_9_re_signed <= signed(dout_2_9_re);

  intdelay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly1 <= dout_2_9_re_signed;
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


  dout_2_9_im_signed <= signed(dout_2_9_im);

  intdelay_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly1 <= dout_2_9_im_signed;
      END IF;
    END IF;
  END PROCESS intdelay_2_process;


  intdelay_3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly2 <= din1_im_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_3_process;


  intdelay_4_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_re_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_re_dly3 <= din1_re_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_4_process;


  intdelay_5_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_im_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_im_dly3 <= din1_im_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_5_process;


  intdelay_6_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_vld_dly1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_vld_dly1 <= dout_2_1_vld;
      END IF;
    END IF;
  END PROCESS intdelay_6_process;


  intdelay_7_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_vld_dly2 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_vld_dly2 <= din1_vld_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_7_process;


  intdelay_8_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din1_vld_dly3 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din1_vld_dly3 <= din1_vld_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_8_process;


  dout_2_10_re_signed <= signed(dout_2_10_re);

  intdelay_9_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_re_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_re_dly1 <= dout_2_10_re_signed;
      END IF;
    END IF;
  END PROCESS intdelay_9_process;


  intdelay_10_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_re_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_re_dly2 <= din2_re_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_10_process;


  dout_2_10_im_signed <= signed(dout_2_10_im);

  intdelay_11_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_im_dly1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_im_dly1 <= dout_2_10_im_signed;
      END IF;
    END IF;
  END PROCESS intdelay_11_process;


  intdelay_12_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_im_dly2 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_im_dly2 <= din2_im_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_12_process;


  intdelay_13_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_re_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_re_dly3 <= din2_re_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_13_process;


  intdelay_14_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din2_im_dly3 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din2_im_dly3 <= din2_im_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_14_process;


  intdelay_15_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      di2_vld_dly1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        di2_vld_dly1 <= dout_2_1_vld;
      END IF;
    END IF;
  END PROCESS intdelay_15_process;


  intdelay_16_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      di2_vld_dly2 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        di2_vld_dly2 <= di2_vld_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_16_process;


  intdelay_17_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      di2_vld_dly3 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        di2_vld_dly3 <= di2_vld_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_17_process;


  twdlXdin_9_re <= twdlXdin_9_re_tmp;

  twdlXdin_9_im <= twdlXdin_9_im_tmp;

  twdlXdin_10_re <= twdlXdin_10_re_tmp;

  twdlXdin_10_im <= twdlXdin_10_im_tmp;

END rtl;

