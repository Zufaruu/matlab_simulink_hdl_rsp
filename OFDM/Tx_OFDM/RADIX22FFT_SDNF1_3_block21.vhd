-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\RADIX22FFT_SDNF1_3_block21.vhd
-- Created: 2023-04-09 19:16:02
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: RADIX22FFT_SDNF1_3_block21
-- Source Path: dsphdl.IFFT/RADIX22FFT_SDNF1_3
-- Hierarchy Level: 3
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY RADIX22FFT_SDNF1_3_block21 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_960_0                       :   IN    std_logic;
        twdlXdin_39_re                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_39_im                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_55_re                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_55_im                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdlXdin_1_vld                    :   IN    std_logic;
        dout_45_re                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_45_im                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_46_re                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_46_im                        :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END RADIX22FFT_SDNF1_3_block21;


ARCHITECTURE rtl OF RADIX22FFT_SDNF1_3_block21 IS

  -- Signals
  SIGNAL twdlXdin_39_re_signed            : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdlXdin_39_im_signed            : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdlXdin_55_re_signed            : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdlXdin_55_im_signed            : signed(15 DOWNTO 0);  -- sfix16_En14
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
  SIGNAL dout_45_re_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_45_im_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_46_re_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_46_im_tmp                   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dout_45_vld                      : std_logic;

BEGIN
  twdlXdin_39_re_signed <= signed(twdlXdin_39_re);

  twdlXdin_39_im_signed <= signed(twdlXdin_39_im);

  twdlXdin_55_re_signed <= signed(twdlXdin_55_re);

  twdlXdin_55_im_signed <= signed(twdlXdin_55_im);

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
       Radix22ButterflyG1_NF_dinXtwdl_vld_dly1, twdlXdin_1_vld,
       twdlXdin_39_im_signed, twdlXdin_39_re_signed, twdlXdin_55_im_signed,
       twdlXdin_55_re_signed)
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
    Radix22ButterflyG1_NF_dinXtwdl_vld_dly1_next <= twdlXdin_1_vld;
    IF twdlXdin_1_vld = '1' THEN 
      add_cast := resize(twdlXdin_39_re_signed, 17);
      add_cast_0 := resize(twdlXdin_55_re_signed, 17);
      Radix22ButterflyG1_NF_btf1_re_reg_next <= add_cast + add_cast_0;
      sub_cast := resize(twdlXdin_39_re_signed, 17);
      sub_cast_0 := resize(twdlXdin_55_re_signed, 17);
      Radix22ButterflyG1_NF_btf2_re_reg_next <= sub_cast - sub_cast_0;
      add_cast_1 := resize(twdlXdin_39_im_signed, 17);
      add_cast_2 := resize(twdlXdin_55_im_signed, 17);
      Radix22ButterflyG1_NF_btf1_im_reg_next <= add_cast_1 + add_cast_2;
      sub_cast_1 := resize(twdlXdin_39_im_signed, 17);
      sub_cast_2 := resize(twdlXdin_55_im_signed, 17);
      Radix22ButterflyG1_NF_btf2_im_reg_next <= sub_cast_1 - sub_cast_2;
    END IF;
    sra_temp := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf1_re_reg, 1);
    dout_45_re_tmp <= sra_temp(15 DOWNTO 0);
    sra_temp_0 := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf1_im_reg, 1);
    dout_45_im_tmp <= sra_temp_0(15 DOWNTO 0);
    sra_temp_1 := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf2_re_reg, 1);
    dout_46_re_tmp <= sra_temp_1(15 DOWNTO 0);
    sra_temp_2 := SHIFT_RIGHT(Radix22ButterflyG1_NF_btf2_im_reg, 1);
    dout_46_im_tmp <= sra_temp_2(15 DOWNTO 0);
    dout_45_vld <= Radix22ButterflyG1_NF_dinXtwdl_vld_dly1;
  END PROCESS Radix22ButterflyG1_NF_output;


  dout_45_re <= std_logic_vector(dout_45_re_tmp);

  dout_45_im <= std_logic_vector(dout_45_im_tmp);

  dout_46_re <= std_logic_vector(dout_46_re_tmp);

  dout_46_im <= std_logic_vector(dout_46_im_tmp);

END rtl;

