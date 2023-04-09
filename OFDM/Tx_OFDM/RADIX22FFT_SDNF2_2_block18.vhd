-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\RADIX22FFT_SDNF2_2_block18.vhd
-- Created: 2023-04-09 12:21:04
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: RADIX22FFT_SDNF2_2_block18
-- Source Path: Tx_OFDM/Transmitter/IFFT/RADIX22FFT_SDNF2_2
-- Hierarchy Level: 2
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY RADIX22FFT_SDNF2_2_block18 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_60_0                        :   IN    std_logic;
        rotate_39                         :   IN    std_logic;  -- ufix1
        dout_8_re                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_8_im                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_40_re                        :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_40_im                        :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_1_vld                        :   IN    std_logic;
        dout_39_re                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_39_im                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_40_re_1                      :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_40_im_1                      :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END RADIX22FFT_SDNF2_2_block18;


ARCHITECTURE rtl OF RADIX22FFT_SDNF2_2_block18 IS

  -- Signals
  SIGNAL dout_8_re_signed                 : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_8_im_signed                 : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_40_re_signed                : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_40_im_signed                : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL Radix22ButterflyG2_NF_din_vld_dly : std_logic;
  SIGNAL Radix22ButterflyG2_NF_btf1_re_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG2_NF_btf1_im_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG2_NF_btf2_re_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG2_NF_btf2_im_reg : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG2_NF_din_vld_dly_next : std_logic;
  SIGNAL Radix22ButterflyG2_NF_btf1_re_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG2_NF_btf1_im_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG2_NF_btf2_re_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG2_NF_btf2_im_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL dout_39_re_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_39_im_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_40_re_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_40_im_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_2_vld                       : std_logic;

BEGIN
  dout_8_re_signed <= signed(dout_8_re);

  dout_8_im_signed <= signed(dout_8_im);

  dout_40_re_signed <= signed(dout_40_re);

  dout_40_im_signed <= signed(dout_40_im);

  -- Radix22ButterflyG2_NF
  Radix22ButterflyG2_NF_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Radix22ButterflyG2_NF_din_vld_dly <= '0';
      Radix22ButterflyG2_NF_btf1_re_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG2_NF_btf1_im_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG2_NF_btf2_re_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG2_NF_btf2_im_reg <= to_signed(16#00000#, 17);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_60_0 = '1' THEN
        Radix22ButterflyG2_NF_din_vld_dly <= Radix22ButterflyG2_NF_din_vld_dly_next;
        Radix22ButterflyG2_NF_btf1_re_reg <= Radix22ButterflyG2_NF_btf1_re_reg_next;
        Radix22ButterflyG2_NF_btf1_im_reg <= Radix22ButterflyG2_NF_btf1_im_reg_next;
        Radix22ButterflyG2_NF_btf2_re_reg <= Radix22ButterflyG2_NF_btf2_re_reg_next;
        Radix22ButterflyG2_NF_btf2_im_reg <= Radix22ButterflyG2_NF_btf2_im_reg_next;
      END IF;
    END IF;
  END PROCESS Radix22ButterflyG2_NF_process;

  Radix22ButterflyG2_NF_output : PROCESS (Radix22ButterflyG2_NF_btf1_im_reg, Radix22ButterflyG2_NF_btf1_re_reg,
       Radix22ButterflyG2_NF_btf2_im_reg, Radix22ButterflyG2_NF_btf2_re_reg,
       Radix22ButterflyG2_NF_din_vld_dly, dout_1_vld, dout_40_im_signed,
       dout_40_re_signed, dout_8_im_signed, dout_8_re_signed, rotate_39)
    VARIABLE add_cast : signed(16 DOWNTO 0);
    VARIABLE add_cast_0 : signed(16 DOWNTO 0);
    VARIABLE add_cast_1 : signed(16 DOWNTO 0);
    VARIABLE add_cast_2 : signed(16 DOWNTO 0);
    VARIABLE sub_cast : signed(16 DOWNTO 0);
    VARIABLE sub_cast_0 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_1 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_2 : signed(16 DOWNTO 0);
    VARIABLE sra_temp : signed(16 DOWNTO 0);
    VARIABLE add_cast_3 : signed(16 DOWNTO 0);
    VARIABLE add_cast_4 : signed(16 DOWNTO 0);
    VARIABLE add_cast_5 : signed(16 DOWNTO 0);
    VARIABLE add_cast_6 : signed(16 DOWNTO 0);
    VARIABLE sra_temp_0 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_3 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_4 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_5 : signed(16 DOWNTO 0);
    VARIABLE sub_cast_6 : signed(16 DOWNTO 0);
    VARIABLE sra_temp_1 : signed(16 DOWNTO 0);
    VARIABLE sra_temp_2 : signed(16 DOWNTO 0);
  BEGIN
    add_cast_1 := to_signed(16#00000#, 17);
    add_cast_2 := to_signed(16#00000#, 17);
    sub_cast_1 := to_signed(16#00000#, 17);
    sub_cast_2 := to_signed(16#00000#, 17);
    add_cast_5 := to_signed(16#00000#, 17);
    add_cast_6 := to_signed(16#00000#, 17);
    sub_cast_5 := to_signed(16#00000#, 17);
    sub_cast_6 := to_signed(16#00000#, 17);
    add_cast := to_signed(16#00000#, 17);
    add_cast_0 := to_signed(16#00000#, 17);
    sub_cast := to_signed(16#00000#, 17);
    sub_cast_0 := to_signed(16#00000#, 17);
    add_cast_3 := to_signed(16#00000#, 17);
    add_cast_4 := to_signed(16#00000#, 17);
    sub_cast_3 := to_signed(16#00000#, 17);
    sub_cast_4 := to_signed(16#00000#, 17);
    Radix22ButterflyG2_NF_btf1_re_reg_next <= Radix22ButterflyG2_NF_btf1_re_reg;
    Radix22ButterflyG2_NF_btf1_im_reg_next <= Radix22ButterflyG2_NF_btf1_im_reg;
    Radix22ButterflyG2_NF_btf2_re_reg_next <= Radix22ButterflyG2_NF_btf2_re_reg;
    Radix22ButterflyG2_NF_btf2_im_reg_next <= Radix22ButterflyG2_NF_btf2_im_reg;
    Radix22ButterflyG2_NF_din_vld_dly_next <= dout_1_vld;
    IF rotate_39 /= '0' THEN 
      IF dout_1_vld = '1' THEN 
        add_cast_1 := resize(dout_8_re_signed, 17);
        add_cast_2 := resize(dout_40_im_signed, 17);
        Radix22ButterflyG2_NF_btf1_re_reg_next <= add_cast_1 + add_cast_2;
        sub_cast_1 := resize(dout_8_re_signed, 17);
        sub_cast_2 := resize(dout_40_im_signed, 17);
        Radix22ButterflyG2_NF_btf2_re_reg_next <= sub_cast_1 - sub_cast_2;
        add_cast_5 := resize(dout_8_im_signed, 17);
        add_cast_6 := resize(dout_40_re_signed, 17);
        Radix22ButterflyG2_NF_btf2_im_reg_next <= add_cast_5 + add_cast_6;
        sub_cast_5 := resize(dout_8_im_signed, 17);
        sub_cast_6 := resize(dout_40_re_signed, 17);
        Radix22ButterflyG2_NF_btf1_im_reg_next <= sub_cast_5 - sub_cast_6;
      END IF;
    ELSIF dout_1_vld = '1' THEN 
      add_cast := resize(dout_8_re_signed, 17);
      add_cast_0 := resize(dout_40_re_signed, 17);
      Radix22ButterflyG2_NF_btf1_re_reg_next <= add_cast + add_cast_0;
      sub_cast := resize(dout_8_re_signed, 17);
      sub_cast_0 := resize(dout_40_re_signed, 17);
      Radix22ButterflyG2_NF_btf2_re_reg_next <= sub_cast - sub_cast_0;
      add_cast_3 := resize(dout_8_im_signed, 17);
      add_cast_4 := resize(dout_40_im_signed, 17);
      Radix22ButterflyG2_NF_btf1_im_reg_next <= add_cast_3 + add_cast_4;
      sub_cast_3 := resize(dout_8_im_signed, 17);
      sub_cast_4 := resize(dout_40_im_signed, 17);
      Radix22ButterflyG2_NF_btf2_im_reg_next <= sub_cast_3 - sub_cast_4;
    END IF;
    sra_temp := SHIFT_RIGHT(Radix22ButterflyG2_NF_btf1_re_reg, 1);
    dout_39_re_tmp <= sra_temp(15 DOWNTO 0);
    sra_temp_0 := SHIFT_RIGHT(Radix22ButterflyG2_NF_btf1_im_reg, 1);
    dout_39_im_tmp <= sra_temp_0(15 DOWNTO 0);
    sra_temp_1 := SHIFT_RIGHT(Radix22ButterflyG2_NF_btf2_re_reg, 1);
    dout_40_re_tmp <= sra_temp_1(15 DOWNTO 0);
    sra_temp_2 := SHIFT_RIGHT(Radix22ButterflyG2_NF_btf2_im_reg, 1);
    dout_40_im_tmp <= sra_temp_2(15 DOWNTO 0);
    dout_2_vld <= Radix22ButterflyG2_NF_din_vld_dly;
  END PROCESS Radix22ButterflyG2_NF_output;


  dout_39_re <= std_logic_vector(dout_39_re_tmp);

  dout_39_im <= std_logic_vector(dout_39_im_tmp);

  dout_40_re_1 <= std_logic_vector(dout_40_re_tmp);

  dout_40_im_1 <= std_logic_vector(dout_40_im_tmp);

END rtl;

