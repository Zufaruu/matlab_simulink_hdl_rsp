-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\RADIX22FFT_SDF1_1_block62.vhd
-- Created: 2023-04-09 15:00:10
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: RADIX22FFT_SDF1_1_block62
-- Source Path: Tx_OFDM/Transmitter/IFFT/RADIX22FFT_SDF1_1
-- Hierarchy Level: 2
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY RADIX22FFT_SDF1_1_block62 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_960_0                       :   IN    std_logic;
        din_1_64_re_dly                   :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        din_1_64_im_dly                   :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        din_1_vld_dly                     :   IN    std_logic;
        rd_1_Addr                         :   IN    std_logic;  -- ufix1
        rd_1_Enb                          :   IN    std_logic;
        proc_1_enb                        :   IN    std_logic;
        dout_1_64_re                      :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        dout_1_64_im                      :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END RADIX22FFT_SDF1_1_block62;


ARCHITECTURE rtl OF RADIX22FFT_SDF1_1_block62 IS

  -- Component Declarations
  COMPONENT SDFCommutator1_block62
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb_1_960_0                     :   IN    std_logic;
          din_1_vld_dly                   :   IN    std_logic;
          xf_re                           :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          xf_im                           :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          xf_vld                          :   IN    std_logic;
          dinXTwdlf_re                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          dinXTwdlf_im                    :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          dinxTwdlf_vld                   :   IN    std_logic;
          btf1_re                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          btf1_im                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          btf2_re                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          btf2_im                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          btf_vld                         :   IN    std_logic;
          wrData_re                       :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          wrData_im                       :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          wrAddr                          :   OUT   std_logic;  -- ufix1
          wrEnb                           :   OUT   std_logic;
          dout_1_64_re                    :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          dout_1_64_im                    :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
          );
  END COMPONENT;

  -- Component Configuration Statements
  FOR ALL : SDFCommutator1_block62
    USE ENTITY work.SDFCommutator1_block62(rtl);

  -- Signals
  SIGNAL din_1_64_re_dly_signed           : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dinXTwdl_re                      : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL din_1_64_im_dly_signed           : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dinXTwdl_im                      : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dinXTwdl_1_64_vld                : std_logic;
  SIGNAL x_vld                            : std_logic;
  SIGNAL btf2_im                          : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL btf2_re                          : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL btf1_im                          : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL btf1_re                          : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dinXTwdlf_im                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL dinXTwdlf_re                     : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL xf_im                            : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL wrData_im                        : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL wrData_im_signed                 : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL wrData_re                        : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL wrData_re_signed                 : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL wrAddr                           : std_logic;  -- ufix1
  SIGNAL wrEnb                            : std_logic;
  SIGNAL twoLocationReg_0_MEM_re_0        : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL twoLocationReg_0_MEM_im_0        : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL twoLocationReg_0_MEM_re_1        : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL twoLocationReg_0_MEM_im_1        : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL twoLocationReg_0_dout_re_reg     : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL twoLocationReg_0_dout_im_reg     : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL twoLocationReg_0_MEM_re_0_next   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twoLocationReg_0_MEM_im_0_next   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twoLocationReg_0_MEM_re_1_next   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twoLocationReg_0_MEM_im_1_next   : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twoLocationReg_0_dout_re_reg_next : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twoLocationReg_0_dout_im_reg_next : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL x_re                             : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL x_im                             : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL Radix22ButterflyG1_btf1_re_reg   : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_btf1_im_reg   : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_btf2_re_reg   : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_btf2_im_reg   : signed(16 DOWNTO 0);  -- sfix17
  SIGNAL Radix22ButterflyG1_x_re_dly1     : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL Radix22ButterflyG1_x_im_dly1     : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL Radix22ButterflyG1_x_vld_dly1    : std_logic;
  SIGNAL Radix22ButterflyG1_dinXtwdl_re_dly1 : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL Radix22ButterflyG1_dinXtwdl_im_dly1 : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL Radix22ButterflyG1_dinXtwdl_re_dly2 : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL Radix22ButterflyG1_dinXtwdl_im_dly2 : signed(15 DOWNTO 0);  -- sfix16
  SIGNAL Radix22ButterflyG1_dinXtwdl_vld_dly1 : std_logic;
  SIGNAL Radix22ButterflyG1_dinXtwdl_vld_dly2 : std_logic;
  SIGNAL Radix22ButterflyG1_btf1_re_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_btf1_im_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_btf2_re_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_btf2_im_reg_next : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_add_cast      : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_add_cast_1    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sub_cast      : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sub_cast_1    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_add_cast_2    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_add_cast_3    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sub_cast_2    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sub_cast_3    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sra_temp      : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sra_temp_1    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sra_temp_2    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL Radix22ButterflyG1_sra_temp_3    : signed(16 DOWNTO 0);  -- sfix17_En14
  SIGNAL xf_re                            : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL xf_vld                           : std_logic;
  SIGNAL dinxTwdlf_vld                    : std_logic;
  SIGNAL btf_vld                          : std_logic;
  SIGNAL dout_1_64_re_tmp                 : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_1_64_im_tmp                 : std_logic_vector(15 DOWNTO 0);  -- ufix16

BEGIN
  u_SDFCOMMUTATOR_1 : SDFCommutator1_block62
    PORT MAP( clk => clk,
              reset => reset,
              enb_1_960_0 => enb_1_960_0,
              din_1_vld_dly => din_1_vld_dly,
              xf_re => std_logic_vector(xf_re),  -- sfix16_En14
              xf_im => std_logic_vector(xf_im),  -- sfix16_En14
              xf_vld => xf_vld,
              dinXTwdlf_re => std_logic_vector(dinXTwdlf_re),  -- sfix16_En14
              dinXTwdlf_im => std_logic_vector(dinXTwdlf_im),  -- sfix16_En14
              dinxTwdlf_vld => dinxTwdlf_vld,
              btf1_re => std_logic_vector(btf1_re),  -- sfix16_En14
              btf1_im => std_logic_vector(btf1_im),  -- sfix16_En14
              btf2_re => std_logic_vector(btf2_re),  -- sfix16_En14
              btf2_im => std_logic_vector(btf2_im),  -- sfix16_En14
              btf_vld => btf_vld,
              wrData_re => wrData_re,  -- sfix16_En14
              wrData_im => wrData_im,  -- sfix16_En14
              wrAddr => wrAddr,  -- ufix1
              wrEnb => wrEnb,
              dout_1_64_re => dout_1_64_re_tmp,  -- sfix16_En14
              dout_1_64_im => dout_1_64_im_tmp  -- sfix16_En14
              );

  din_1_64_re_dly_signed <= signed(din_1_64_re_dly);

  intdelay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      dinXTwdl_re <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        dinXTwdl_re <= din_1_64_re_dly_signed;
      END IF;
    END IF;
  END PROCESS intdelay_process;


  din_1_64_im_dly_signed <= signed(din_1_64_im_dly);

  intdelay_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      dinXTwdl_im <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        dinXTwdl_im <= din_1_64_im_dly_signed;
      END IF;
    END IF;
  END PROCESS intdelay_1_process;


  intdelay_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      dinXTwdl_1_64_vld <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        dinXTwdl_1_64_vld <= din_1_vld_dly;
      END IF;
    END IF;
  END PROCESS intdelay_2_process;


  intdelay_3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      x_vld <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        x_vld <= rd_1_Enb;
      END IF;
    END IF;
  END PROCESS intdelay_3_process;


  wrData_im_signed <= signed(wrData_im);

  wrData_re_signed <= signed(wrData_re);

  -- twoLocationReg_0
  twoLocationReg_0_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twoLocationReg_0_MEM_re_0 <= to_signed(16#0000#, 16);
      twoLocationReg_0_MEM_im_0 <= to_signed(16#0000#, 16);
      twoLocationReg_0_MEM_re_1 <= to_signed(16#0000#, 16);
      twoLocationReg_0_MEM_im_1 <= to_signed(16#0000#, 16);
      twoLocationReg_0_dout_re_reg <= to_signed(16#0000#, 16);
      twoLocationReg_0_dout_im_reg <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twoLocationReg_0_MEM_re_0 <= twoLocationReg_0_MEM_re_0_next;
        twoLocationReg_0_MEM_im_0 <= twoLocationReg_0_MEM_im_0_next;
        twoLocationReg_0_MEM_re_1 <= twoLocationReg_0_MEM_re_1_next;
        twoLocationReg_0_MEM_im_1 <= twoLocationReg_0_MEM_im_1_next;
        twoLocationReg_0_dout_re_reg <= twoLocationReg_0_dout_re_reg_next;
        twoLocationReg_0_dout_im_reg <= twoLocationReg_0_dout_im_reg_next;
      END IF;
    END IF;
  END PROCESS twoLocationReg_0_process;

  twoLocationReg_0_output : PROCESS (rd_1_Addr, twoLocationReg_0_MEM_im_0, twoLocationReg_0_MEM_im_1,
       twoLocationReg_0_MEM_re_0, twoLocationReg_0_MEM_re_1,
       twoLocationReg_0_dout_im_reg, twoLocationReg_0_dout_re_reg, wrAddr,
       wrData_im_signed, wrData_re_signed, wrEnb)
  BEGIN
    twoLocationReg_0_MEM_re_0_next <= twoLocationReg_0_MEM_re_0;
    twoLocationReg_0_MEM_im_0_next <= twoLocationReg_0_MEM_im_0;
    twoLocationReg_0_MEM_re_1_next <= twoLocationReg_0_MEM_re_1;
    twoLocationReg_0_MEM_im_1_next <= twoLocationReg_0_MEM_im_1;
    IF rd_1_Addr = '1' THEN 
      twoLocationReg_0_dout_re_reg_next <= twoLocationReg_0_MEM_re_1;
      twoLocationReg_0_dout_im_reg_next <= twoLocationReg_0_MEM_im_1;
    ELSE 
      twoLocationReg_0_dout_re_reg_next <= twoLocationReg_0_MEM_re_0;
      twoLocationReg_0_dout_im_reg_next <= twoLocationReg_0_MEM_im_0;
    END IF;
    IF wrEnb = '1' THEN 
      IF wrAddr = '1' THEN 
        twoLocationReg_0_MEM_re_1_next <= wrData_re_signed;
        twoLocationReg_0_MEM_im_1_next <= wrData_im_signed;
      ELSE 
        twoLocationReg_0_MEM_re_0_next <= wrData_re_signed;
        twoLocationReg_0_MEM_im_0_next <= wrData_im_signed;
      END IF;
    END IF;
    x_re <= twoLocationReg_0_dout_re_reg;
    x_im <= twoLocationReg_0_dout_im_reg;
  END PROCESS twoLocationReg_0_output;


  -- Radix22ButterflyG1
  Radix22ButterflyG1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Radix22ButterflyG1_btf1_re_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_btf1_im_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_btf2_re_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_btf2_im_reg <= to_signed(16#00000#, 17);
      Radix22ButterflyG1_x_re_dly1 <= to_signed(16#0000#, 16);
      Radix22ButterflyG1_x_im_dly1 <= to_signed(16#0000#, 16);
      Radix22ButterflyG1_x_vld_dly1 <= '0';
      xf_re <= to_signed(16#0000#, 16);
      xf_im <= to_signed(16#0000#, 16);
      xf_vld <= '0';
      Radix22ButterflyG1_dinXtwdl_re_dly1 <= to_signed(16#0000#, 16);
      Radix22ButterflyG1_dinXtwdl_im_dly1 <= to_signed(16#0000#, 16);
      Radix22ButterflyG1_dinXtwdl_re_dly2 <= to_signed(16#0000#, 16);
      Radix22ButterflyG1_dinXtwdl_im_dly2 <= to_signed(16#0000#, 16);
      Radix22ButterflyG1_dinXtwdl_vld_dly1 <= '0';
      Radix22ButterflyG1_dinXtwdl_vld_dly2 <= '0';
      btf_vld <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        Radix22ButterflyG1_btf1_re_reg <= Radix22ButterflyG1_btf1_re_reg_next;
        Radix22ButterflyG1_btf1_im_reg <= Radix22ButterflyG1_btf1_im_reg_next;
        Radix22ButterflyG1_btf2_re_reg <= Radix22ButterflyG1_btf2_re_reg_next;
        Radix22ButterflyG1_btf2_im_reg <= Radix22ButterflyG1_btf2_im_reg_next;
        xf_re <= Radix22ButterflyG1_x_re_dly1;
        xf_im <= Radix22ButterflyG1_x_im_dly1;
        xf_vld <= Radix22ButterflyG1_x_vld_dly1;
        btf_vld <= Radix22ButterflyG1_dinXtwdl_vld_dly2;
        Radix22ButterflyG1_dinXtwdl_vld_dly2 <= Radix22ButterflyG1_dinXtwdl_vld_dly1;
        Radix22ButterflyG1_dinXtwdl_re_dly2 <= Radix22ButterflyG1_dinXtwdl_re_dly1;
        Radix22ButterflyG1_dinXtwdl_im_dly2 <= Radix22ButterflyG1_dinXtwdl_im_dly1;
        Radix22ButterflyG1_dinXtwdl_re_dly1 <= dinXTwdl_re;
        Radix22ButterflyG1_dinXtwdl_im_dly1 <= dinXTwdl_im;
        Radix22ButterflyG1_x_re_dly1 <= x_re;
        Radix22ButterflyG1_x_im_dly1 <= x_im;
        Radix22ButterflyG1_x_vld_dly1 <= x_vld;
        Radix22ButterflyG1_dinXtwdl_vld_dly1 <= proc_1_enb AND dinXTwdl_1_64_vld;
      END IF;
    END IF;
  END PROCESS Radix22ButterflyG1_process;

  dinxTwdlf_vld <= ( NOT proc_1_enb) AND dinXTwdl_1_64_vld;
  Radix22ButterflyG1_add_cast <= resize(Radix22ButterflyG1_x_re_dly1, 17);
  Radix22ButterflyG1_add_cast_1 <= resize(Radix22ButterflyG1_dinXtwdl_re_dly2, 17);
  Radix22ButterflyG1_btf1_re_reg_next <= Radix22ButterflyG1_add_cast + Radix22ButterflyG1_add_cast_1;
  Radix22ButterflyG1_sub_cast <= resize(Radix22ButterflyG1_x_re_dly1, 17);
  Radix22ButterflyG1_sub_cast_1 <= resize(Radix22ButterflyG1_dinXtwdl_re_dly2, 17);
  Radix22ButterflyG1_btf2_re_reg_next <= Radix22ButterflyG1_sub_cast - Radix22ButterflyG1_sub_cast_1;
  Radix22ButterflyG1_add_cast_2 <= resize(Radix22ButterflyG1_x_im_dly1, 17);
  Radix22ButterflyG1_add_cast_3 <= resize(Radix22ButterflyG1_dinXtwdl_im_dly2, 17);
  Radix22ButterflyG1_btf1_im_reg_next <= Radix22ButterflyG1_add_cast_2 + Radix22ButterflyG1_add_cast_3;
  Radix22ButterflyG1_sub_cast_2 <= resize(Radix22ButterflyG1_x_im_dly1, 17);
  Radix22ButterflyG1_sub_cast_3 <= resize(Radix22ButterflyG1_dinXtwdl_im_dly2, 17);
  Radix22ButterflyG1_btf2_im_reg_next <= Radix22ButterflyG1_sub_cast_2 - Radix22ButterflyG1_sub_cast_3;
  dinXTwdlf_re <= dinXTwdl_re;
  dinXTwdlf_im <= dinXTwdl_im;
  Radix22ButterflyG1_sra_temp <= SHIFT_RIGHT(Radix22ButterflyG1_btf1_re_reg, 1);
  btf1_re <= Radix22ButterflyG1_sra_temp(15 DOWNTO 0);
  Radix22ButterflyG1_sra_temp_1 <= SHIFT_RIGHT(Radix22ButterflyG1_btf1_im_reg, 1);
  btf1_im <= Radix22ButterflyG1_sra_temp_1(15 DOWNTO 0);
  Radix22ButterflyG1_sra_temp_2 <= SHIFT_RIGHT(Radix22ButterflyG1_btf2_re_reg, 1);
  btf2_re <= Radix22ButterflyG1_sra_temp_2(15 DOWNTO 0);
  Radix22ButterflyG1_sra_temp_3 <= SHIFT_RIGHT(Radix22ButterflyG1_btf2_im_reg, 1);
  btf2_im <= Radix22ButterflyG1_sra_temp_3(15 DOWNTO 0);

  dout_1_64_re <= dout_1_64_re_tmp;

  dout_1_64_im <= dout_1_64_im_tmp;

END rtl;

