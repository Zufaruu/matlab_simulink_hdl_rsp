-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\Complex4Multiply_block12.vhd
-- Created: 2023-04-09 19:16:02
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: Complex4Multiply_block12
-- Source Path: dsphdl.IFFT/TWDLMULT_SDNF1_3/Complex4Multiply
-- Hierarchy Level: 4
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY Complex4Multiply_block12 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_960_0                       :   IN    std_logic;
        din1_re_dly3                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        din1_im_dly3                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        din1_vld_dly3                     :   IN    std_logic;
        twdl_3_15_re                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdl_3_15_im                      :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_15_re                    :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_15_im                    :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END Complex4Multiply_block12;


ARCHITECTURE rtl OF Complex4Multiply_block12 IS

  -- Signals
  SIGNAL din1_re_dly3_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din_re_reg                       : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din1_im_dly3_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din_im_reg                       : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdl_3_15_re_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdl_re_reg                      : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdl_3_15_im_signed              : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdl_im_reg                      : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL Complex4Multiply_din1_re_pipe1   : signed(15 DOWNTO 0) := to_signed(16#0000#, 16);  -- sfix16
  SIGNAL Complex4Multiply_din1_im_pipe1   : signed(15 DOWNTO 0) := to_signed(16#0000#, 16);  -- sfix16
  SIGNAL Complex4Multiply_mult1_re_pipe1  : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32
  SIGNAL Complex4Multiply_mult2_re_pipe1  : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32
  SIGNAL Complex4Multiply_mult1_im_pipe1  : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32
  SIGNAL Complex4Multiply_mult2_im_pipe1  : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32
  SIGNAL Complex4Multiply_twiddle_re_pipe1 : signed(15 DOWNTO 0) := to_signed(16#0000#, 16);  -- sfix16
  SIGNAL Complex4Multiply_twiddle_im_pipe1 : signed(15 DOWNTO 0) := to_signed(16#0000#, 16);  -- sfix16
  SIGNAL prod1_re                         : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32_En28
  SIGNAL prod1_im                         : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32_En28
  SIGNAL prod2_re                         : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32_En28
  SIGNAL prod2_im                         : signed(31 DOWNTO 0) := to_signed(0, 32);  -- sfix32_En28
  SIGNAL din_vld_dly1                     : std_logic;
  SIGNAL din_vld_dly2                     : std_logic;
  SIGNAL din_vld_dly3                     : std_logic;
  SIGNAL prod_vld                         : std_logic;
  SIGNAL Complex4Add_multRes_re_reg       : signed(32 DOWNTO 0);  -- sfix33
  SIGNAL Complex4Add_multRes_im_reg       : signed(32 DOWNTO 0);  -- sfix33
  SIGNAL Complex4Add_prod_vld_reg1        : std_logic;
  SIGNAL Complex4Add_prod1_re_reg         : signed(31 DOWNTO 0);  -- sfix32
  SIGNAL Complex4Add_prod1_im_reg         : signed(31 DOWNTO 0);  -- sfix32
  SIGNAL Complex4Add_prod2_re_reg         : signed(31 DOWNTO 0);  -- sfix32
  SIGNAL Complex4Add_prod2_im_reg         : signed(31 DOWNTO 0);  -- sfix32
  SIGNAL Complex4Add_multRes_re_reg_next  : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL Complex4Add_multRes_im_reg_next  : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL Complex4Add_sub_cast             : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL Complex4Add_sub_cast_1           : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL Complex4Add_add_cast             : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL Complex4Add_add_cast_1           : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL mulResFP_re                      : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL mulResFP_im                      : signed(32 DOWNTO 0);  -- sfix33_En28
  SIGNAL twdlXdin1_vld                    : std_logic;
  SIGNAL twdlXdin_15_re_tmp               : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdlXdin_15_im_tmp               : signed(15 DOWNTO 0);  -- sfix16_En14

BEGIN
  din1_re_dly3_signed <= signed(din1_re_dly3);

  intdelay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_re_reg <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_re_reg <= din1_re_dly3_signed;
      END IF;
    END IF;
  END PROCESS intdelay_process;


  din1_im_dly3_signed <= signed(din1_im_dly3);

  intdelay_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_im_reg <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_im_reg <= din1_im_dly3_signed;
      END IF;
    END IF;
  END PROCESS intdelay_1_process;


  twdl_3_15_re_signed <= signed(twdl_3_15_re);

  intdelay_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twdl_re_reg <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twdl_re_reg <= twdl_3_15_re_signed;
      END IF;
    END IF;
  END PROCESS intdelay_2_process;


  twdl_3_15_im_signed <= signed(twdl_3_15_im);

  intdelay_3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twdl_im_reg <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twdl_im_reg <= twdl_3_15_im_signed;
      END IF;
    END IF;
  END PROCESS intdelay_3_process;


  -- Complex4Multiply
  Complex4Multiply_process : PROCESS (clk)
  BEGIN
    IF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
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
      END IF;
    END IF;
  END PROCESS Complex4Multiply_process;


  intdelay_4_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly1 <= din1_vld_dly3;
      END IF;
    END IF;
  END PROCESS intdelay_4_process;


  intdelay_5_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly2 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly2 <= din_vld_dly1;
      END IF;
    END IF;
  END PROCESS intdelay_5_process;


  intdelay_6_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      din_vld_dly3 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        din_vld_dly3 <= din_vld_dly2;
      END IF;
    END IF;
  END PROCESS intdelay_6_process;


  intdelay_7_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      prod_vld <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        prod_vld <= din_vld_dly3;
      END IF;
    END IF;
  END PROCESS intdelay_7_process;


  -- Complex4Add
  Complex4Add_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Complex4Add_multRes_re_reg <= to_signed(0, 33);
      Complex4Add_multRes_im_reg <= to_signed(0, 33);
      Complex4Add_prod1_re_reg <= to_signed(0, 32);
      Complex4Add_prod1_im_reg <= to_signed(0, 32);
      Complex4Add_prod2_re_reg <= to_signed(0, 32);
      Complex4Add_prod2_im_reg <= to_signed(0, 32);
      Complex4Add_prod_vld_reg1 <= '0';
      twdlXdin1_vld <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        Complex4Add_multRes_re_reg <= Complex4Add_multRes_re_reg_next;
        Complex4Add_multRes_im_reg <= Complex4Add_multRes_im_reg_next;
        Complex4Add_prod1_re_reg <= prod1_re;
        Complex4Add_prod1_im_reg <= prod1_im;
        Complex4Add_prod2_re_reg <= prod2_re;
        Complex4Add_prod2_im_reg <= prod2_im;
        twdlXdin1_vld <= Complex4Add_prod_vld_reg1;
        Complex4Add_prod_vld_reg1 <= prod_vld;
      END IF;
    END IF;
  END PROCESS Complex4Add_process;

  Complex4Add_sub_cast <= resize(Complex4Add_prod1_re_reg, 33);
  Complex4Add_sub_cast_1 <= resize(Complex4Add_prod2_re_reg, 33);
  Complex4Add_multRes_re_reg_next <= Complex4Add_sub_cast - Complex4Add_sub_cast_1;
  Complex4Add_add_cast <= resize(Complex4Add_prod1_im_reg, 33);
  Complex4Add_add_cast_1 <= resize(Complex4Add_prod2_im_reg, 33);
  Complex4Add_multRes_im_reg_next <= Complex4Add_add_cast + Complex4Add_add_cast_1;
  mulResFP_re <= Complex4Add_multRes_re_reg;
  mulResFP_im <= Complex4Add_multRes_im_reg;

  twdlXdin_15_re_tmp <= mulResFP_re(29 DOWNTO 14);

  twdlXdin_15_re <= std_logic_vector(twdlXdin_15_re_tmp);

  twdlXdin_15_im_tmp <= mulResFP_im(29 DOWNTO 14);

  twdlXdin_15_im <= std_logic_vector(twdlXdin_15_im_tmp);

END rtl;

