-- -------------------------------------------------------------
-- 
-- File Name: hdlsrc\Tx_OFDM\IFFT1.vhd
-- Created: 2023-04-08 13:30:41
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: IFFT1
-- Source Path: Tx_OFDM/Transmitter/IFFT1
-- Hierarchy Level: 1
-- 
-- FFT
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;
USE work.Transmitter_pkg.ALL;

ENTITY IFFT1 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb                               :   IN    std_logic;
        dataIn_re                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
        dataIn_im                         :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
        validIn                           :   IN    std_logic;
        dataOut_re                        :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
        dataOut_im                        :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En15
        );
END IFFT1;


ARCHITECTURE rtl OF IFFT1 IS

  -- Component Declarations
  COMPONENT RADIX22FFT_CTRL1_1
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dinXTwdl_1_1_vld                :   IN    std_logic;
          dinXTwdl_1_1_vld_1              :   IN    std_logic;
          rd_1_Addr                       :   OUT   std_logic_vector(4 DOWNTO 0);  -- ufix5
          rd_1_Enb                        :   OUT   std_logic;
          proc_1_enb                      :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_SDF1_1
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          din_1_1_re_dly                  :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          din_1_1_im_dly                  :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          din_1_vld_dly                   :   IN    std_logic;
          rd_1_Addr                       :   IN    std_logic_vector(4 DOWNTO 0);  -- ufix5
          rd_1_Enb                        :   IN    std_logic;
          proc_1_enb                      :   IN    std_logic;
          dout_1_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_1_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_1_1_vld                    :   OUT   std_logic;
          dinXTwdl_1_1_vld                :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_CTRL1_2
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_1_1_vld                    :   IN    std_logic;
          dinXTwdl_2_1_vld                :   IN    std_logic;
          rd_2_Addr                       :   OUT   std_logic_vector(3 DOWNTO 0);  -- ufix4
          rd_2_Enb                        :   OUT   std_logic;
          proc_2_enb                      :   OUT   std_logic;
          multiply_2_J                    :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_SDF2_2
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_1_1_re                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_1_1_im                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_1_1_vld                    :   IN    std_logic;
          rd_2_Addr                       :   IN    std_logic_vector(3 DOWNTO 0);  -- ufix4
          rd_2_Enb                        :   IN    std_logic;
          proc_2_enb                      :   IN    std_logic;
          multiply_2_J                    :   IN    std_logic;
          dout_2_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_2_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_2_1_vld                    :   OUT   std_logic;
          dinXTwdl_2_1_vld                :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT TWDLROM_3_1
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_2_1_vld                    :   IN    std_logic;
          twdl_3_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_3_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_CTRL1_3
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dinXTwdl_3_1_vld                :   IN    std_logic;
          dinXTwdl_3_1_vld_1              :   IN    std_logic;
          rd_3_Addr                       :   OUT   std_logic_vector(2 DOWNTO 0);  -- ufix3
          rd_3_Enb                        :   OUT   std_logic;
          proc_3_enb                      :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_SDF1_3
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          din_3_1_re_dly                  :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          din_3_1_im_dly                  :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          din_3_vld_dly                   :   IN    std_logic;
          rd_3_Addr                       :   IN    std_logic_vector(2 DOWNTO 0);  -- ufix3
          rd_3_Enb                        :   IN    std_logic;
          twdl_3_1_re                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_3_1_im                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          proc_3_enb                      :   IN    std_logic;
          dout_3_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_3_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_3_1_vld                    :   OUT   std_logic;
          dinXTwdl_3_1_vld                :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_CTRL1_4
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_3_1_vld                    :   IN    std_logic;
          dinXTwdl_4_1_vld                :   IN    std_logic;
          rd_4_Addr                       :   OUT   std_logic_vector(1 DOWNTO 0);  -- ufix2
          rd_4_Enb                        :   OUT   std_logic;
          proc_4_enb                      :   OUT   std_logic;
          multiply_4_J                    :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_SDF2_4
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_3_1_re                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_3_1_im                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_3_1_vld                    :   IN    std_logic;
          rd_4_Addr                       :   IN    std_logic_vector(1 DOWNTO 0);  -- ufix2
          rd_4_Enb                        :   IN    std_logic;
          proc_4_enb                      :   IN    std_logic;
          multiply_4_J                    :   IN    std_logic;
          dout_4_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_4_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_4_1_vld                    :   OUT   std_logic;
          dinXTwdl_4_1_vld                :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT TWDLROM_5_1
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_4_1_vld                    :   IN    std_logic;
          twdl_5_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_5_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_CTRL1_5
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dinXTwdl_5_1_vld                :   IN    std_logic;
          dinXTwdl_5_1_vld_1              :   IN    std_logic;
          rd_5_Addr                       :   OUT   std_logic;  -- ufix1
          rd_5_Enb                        :   OUT   std_logic;
          proc_5_enb                      :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_SDF1_5
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          din_5_1_re_dly                  :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          din_5_1_im_dly                  :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          din_5_vld_dly                   :   IN    std_logic;
          rd_5_Addr                       :   IN    std_logic;  -- ufix1
          rd_5_Enb                        :   IN    std_logic;
          twdl_5_1_re                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          twdl_5_1_im                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
          proc_5_enb                      :   IN    std_logic;
          dout_5_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_5_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_5_1_vld                    :   OUT   std_logic;
          dinXTwdl_5_1_vld                :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_CTRL1_6
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_5_1_vld                    :   IN    std_logic;
          dinXTwdl_6_1_vld                :   IN    std_logic;
          rd_6_Addr                       :   OUT   std_logic;
          rd_6_Enb                        :   OUT   std_logic;
          proc_6_enb                      :   OUT   std_logic;
          multiply_6_J                    :   OUT   std_logic
          );
  END COMPONENT;

  COMPONENT RADIX22FFT_SDF2_6
    PORT( clk                             :   IN    std_logic;
          reset                           :   IN    std_logic;
          enb                             :   IN    std_logic;
          dout_5_1_re                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_5_1_im                     :   IN    std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_5_1_vld                    :   IN    std_logic;
          rd_6_Addr                       :   IN    std_logic;
          rd_6_Enb                        :   IN    std_logic;
          proc_6_enb                      :   IN    std_logic;
          multiply_6_J                    :   IN    std_logic;
          dout_6_1_re                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dout_6_1_im                     :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En15
          dinXTwdl_6_1_vld                :   OUT   std_logic
          );
  END COMPONENT;

  -- Component Configuration Statements
  FOR ALL : RADIX22FFT_CTRL1_1
    USE ENTITY work.RADIX22FFT_CTRL1_1(rtl);

  FOR ALL : RADIX22FFT_SDF1_1
    USE ENTITY work.RADIX22FFT_SDF1_1(rtl);

  FOR ALL : RADIX22FFT_CTRL1_2
    USE ENTITY work.RADIX22FFT_CTRL1_2(rtl);

  FOR ALL : RADIX22FFT_SDF2_2
    USE ENTITY work.RADIX22FFT_SDF2_2(rtl);

  FOR ALL : TWDLROM_3_1
    USE ENTITY work.TWDLROM_3_1(rtl);

  FOR ALL : RADIX22FFT_CTRL1_3
    USE ENTITY work.RADIX22FFT_CTRL1_3(rtl);

  FOR ALL : RADIX22FFT_SDF1_3
    USE ENTITY work.RADIX22FFT_SDF1_3(rtl);

  FOR ALL : RADIX22FFT_CTRL1_4
    USE ENTITY work.RADIX22FFT_CTRL1_4(rtl);

  FOR ALL : RADIX22FFT_SDF2_4
    USE ENTITY work.RADIX22FFT_SDF2_4(rtl);

  FOR ALL : TWDLROM_5_1
    USE ENTITY work.TWDLROM_5_1(rtl);

  FOR ALL : RADIX22FFT_CTRL1_5
    USE ENTITY work.RADIX22FFT_CTRL1_5(rtl);

  FOR ALL : RADIX22FFT_SDF1_5
    USE ENTITY work.RADIX22FFT_SDF1_5(rtl);

  FOR ALL : RADIX22FFT_CTRL1_6
    USE ENTITY work.RADIX22FFT_CTRL1_6(rtl);

  FOR ALL : RADIX22FFT_SDF2_6
    USE ENTITY work.RADIX22FFT_SDF2_6(rtl);

  -- Signals
  SIGNAL softReset                        : std_logic;
  SIGNAL dataIn_im_signed                 : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL dataIn_re_signed                 : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg                     : vector_of_signed16(0 TO 2);  -- sfix16 [3]
  SIGNAL intdelay_reg_next                : vector_of_signed16(0 TO 2);  -- sfix16_En15 [3]
  SIGNAL din_1_1_re_dly                   : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_1                   : vector_of_signed16(0 TO 2);  -- sfix16 [3]
  SIGNAL intdelay_reg_next_1              : vector_of_signed16(0 TO 2);  -- sfix16_En15 [3]
  SIGNAL din_1_1_im_dly                   : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_2                   : std_logic_vector(0 TO 2);  -- ufix1 [3]
  SIGNAL intdelay_reg_next_2              : std_logic_vector(0 TO 2);  -- ufix1 [3]
  SIGNAL din_1_vld_dly                    : std_logic;
  SIGNAL dinXTwdl_1_1_vld                 : std_logic;
  SIGNAL rd_1_Addr                        : std_logic_vector(4 DOWNTO 0);  -- ufix5
  SIGNAL rd_1_Enb                         : std_logic;
  SIGNAL proc_1_enb                       : std_logic;
  SIGNAL dout_1_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_1_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_1_1_vld                     : std_logic;
  SIGNAL dinXTwdl_2_1_vld                 : std_logic;
  SIGNAL rd_2_Addr                        : std_logic_vector(3 DOWNTO 0);  -- ufix4
  SIGNAL rd_2_Enb                         : std_logic;
  SIGNAL proc_2_enb                       : std_logic;
  SIGNAL multiply_2_J                     : std_logic;
  SIGNAL dout_2_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_2_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_2_1_vld                     : std_logic;
  SIGNAL dout_2_1_re_signed               : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL dout_2_1_im_signed               : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_3                   : vector_of_signed16(0 TO 2);  -- sfix16 [3]
  SIGNAL intdelay_reg_next_3              : vector_of_signed16(0 TO 2);  -- sfix16_En15 [3]
  SIGNAL din_3_1_re_dly                   : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_4                   : vector_of_signed16(0 TO 2);  -- sfix16 [3]
  SIGNAL intdelay_reg_next_4              : vector_of_signed16(0 TO 2);  -- sfix16_En15 [3]
  SIGNAL din_3_1_im_dly                   : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_5                   : std_logic_vector(0 TO 2);  -- ufix1 [3]
  SIGNAL intdelay_reg_next_5              : std_logic_vector(0 TO 2);  -- ufix1 [3]
  SIGNAL din_3_vld_dly                    : std_logic;
  SIGNAL twdl_3_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL twdl_3_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dinXTwdl_3_1_vld                 : std_logic;
  SIGNAL rd_3_Addr                        : std_logic_vector(2 DOWNTO 0);  -- ufix3
  SIGNAL rd_3_Enb                         : std_logic;
  SIGNAL proc_3_enb                       : std_logic;
  SIGNAL dout_3_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_3_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_3_1_vld                     : std_logic;
  SIGNAL dinXTwdl_4_1_vld                 : std_logic;
  SIGNAL rd_4_Addr                        : std_logic_vector(1 DOWNTO 0);  -- ufix2
  SIGNAL rd_4_Enb                         : std_logic;
  SIGNAL proc_4_enb                       : std_logic;
  SIGNAL multiply_4_J                     : std_logic;
  SIGNAL dout_4_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_4_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_4_1_vld                     : std_logic;
  SIGNAL dout_4_1_re_signed               : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL dout_4_1_im_signed               : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_6                   : vector_of_signed16(0 TO 2);  -- sfix16 [3]
  SIGNAL intdelay_reg_next_6              : vector_of_signed16(0 TO 2);  -- sfix16_En15 [3]
  SIGNAL din_5_1_re_dly                   : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_7                   : vector_of_signed16(0 TO 2);  -- sfix16 [3]
  SIGNAL intdelay_reg_next_7              : vector_of_signed16(0 TO 2);  -- sfix16_En15 [3]
  SIGNAL din_5_1_im_dly                   : signed(15 DOWNTO 0);  -- sfix16_En15
  SIGNAL intdelay_reg_8                   : std_logic_vector(0 TO 2);  -- ufix1 [3]
  SIGNAL intdelay_reg_next_8              : std_logic_vector(0 TO 2);  -- ufix1 [3]
  SIGNAL din_5_vld_dly                    : std_logic;
  SIGNAL twdl_5_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL twdl_5_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dinXTwdl_5_1_vld                 : std_logic;
  SIGNAL rd_5_Addr                        : std_logic;  -- ufix1
  SIGNAL rd_5_Enb                         : std_logic;
  SIGNAL proc_5_enb                       : std_logic;
  SIGNAL dout_5_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_5_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_5_1_vld                     : std_logic;
  SIGNAL dinXTwdl_6_1_vld                 : std_logic;
  SIGNAL rd_6_Addr                        : std_logic;
  SIGNAL rd_6_Enb                         : std_logic;
  SIGNAL proc_6_enb                       : std_logic;
  SIGNAL multiply_6_J                     : std_logic;
  SIGNAL dout_6_1_re                      : std_logic_vector(15 DOWNTO 0);  -- ufix16
  SIGNAL dout_6_1_im                      : std_logic_vector(15 DOWNTO 0);  -- ufix16

BEGIN
  u_CTRL1_1_1 : RADIX22FFT_CTRL1_1
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dinXTwdl_1_1_vld => dinXTwdl_1_1_vld,
              dinXTwdl_1_1_vld_1 => dinXTwdl_1_1_vld,
              rd_1_Addr => rd_1_Addr,  -- ufix5
              rd_1_Enb => rd_1_Enb,
              proc_1_enb => proc_1_enb
              );

  u_SDF1_1_1 : RADIX22FFT_SDF1_1
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              din_1_1_re_dly => std_logic_vector(din_1_1_re_dly),  -- sfix16_En15
              din_1_1_im_dly => std_logic_vector(din_1_1_im_dly),  -- sfix16_En15
              din_1_vld_dly => din_1_vld_dly,
              rd_1_Addr => rd_1_Addr,  -- ufix5
              rd_1_Enb => rd_1_Enb,
              proc_1_enb => proc_1_enb,
              dout_1_1_re => dout_1_1_re,  -- sfix16_En15
              dout_1_1_im => dout_1_1_im,  -- sfix16_En15
              dout_1_1_vld => dout_1_1_vld,
              dinXTwdl_1_1_vld => dinXTwdl_1_1_vld
              );

  u_CTRL2_2_1 : RADIX22FFT_CTRL1_2
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_1_1_vld => dout_1_1_vld,
              dinXTwdl_2_1_vld => dinXTwdl_2_1_vld,
              rd_2_Addr => rd_2_Addr,  -- ufix4
              rd_2_Enb => rd_2_Enb,
              proc_2_enb => proc_2_enb,
              multiply_2_J => multiply_2_J
              );

  u_SDF2_2_1 : RADIX22FFT_SDF2_2
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_1_1_re => dout_1_1_re,  -- sfix16_En15
              dout_1_1_im => dout_1_1_im,  -- sfix16_En15
              dout_1_1_vld => dout_1_1_vld,
              rd_2_Addr => rd_2_Addr,  -- ufix4
              rd_2_Enb => rd_2_Enb,
              proc_2_enb => proc_2_enb,
              multiply_2_J => multiply_2_J,
              dout_2_1_re => dout_2_1_re,  -- sfix16_En15
              dout_2_1_im => dout_2_1_im,  -- sfix16_En15
              dout_2_1_vld => dout_2_1_vld,
              dinXTwdl_2_1_vld => dinXTwdl_2_1_vld
              );

  u_twdlROM_3_1 : TWDLROM_3_1
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_2_1_vld => dout_2_1_vld,
              twdl_3_1_re => twdl_3_1_re,  -- sfix16_En14
              twdl_3_1_im => twdl_3_1_im  -- sfix16_En14
              );

  u_CTRL1_3_1 : RADIX22FFT_CTRL1_3
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dinXTwdl_3_1_vld => dinXTwdl_3_1_vld,
              dinXTwdl_3_1_vld_1 => dinXTwdl_3_1_vld,
              rd_3_Addr => rd_3_Addr,  -- ufix3
              rd_3_Enb => rd_3_Enb,
              proc_3_enb => proc_3_enb
              );

  u_SDF1_3_1 : RADIX22FFT_SDF1_3
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              din_3_1_re_dly => std_logic_vector(din_3_1_re_dly),  -- sfix16_En15
              din_3_1_im_dly => std_logic_vector(din_3_1_im_dly),  -- sfix16_En15
              din_3_vld_dly => din_3_vld_dly,
              rd_3_Addr => rd_3_Addr,  -- ufix3
              rd_3_Enb => rd_3_Enb,
              twdl_3_1_re => twdl_3_1_re,  -- sfix16_En14
              twdl_3_1_im => twdl_3_1_im,  -- sfix16_En14
              proc_3_enb => proc_3_enb,
              dout_3_1_re => dout_3_1_re,  -- sfix16_En15
              dout_3_1_im => dout_3_1_im,  -- sfix16_En15
              dout_3_1_vld => dout_3_1_vld,
              dinXTwdl_3_1_vld => dinXTwdl_3_1_vld
              );

  u_CTRL2_4_1 : RADIX22FFT_CTRL1_4
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_3_1_vld => dout_3_1_vld,
              dinXTwdl_4_1_vld => dinXTwdl_4_1_vld,
              rd_4_Addr => rd_4_Addr,  -- ufix2
              rd_4_Enb => rd_4_Enb,
              proc_4_enb => proc_4_enb,
              multiply_4_J => multiply_4_J
              );

  u_SDF2_4_1 : RADIX22FFT_SDF2_4
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_3_1_re => dout_3_1_re,  -- sfix16_En15
              dout_3_1_im => dout_3_1_im,  -- sfix16_En15
              dout_3_1_vld => dout_3_1_vld,
              rd_4_Addr => rd_4_Addr,  -- ufix2
              rd_4_Enb => rd_4_Enb,
              proc_4_enb => proc_4_enb,
              multiply_4_J => multiply_4_J,
              dout_4_1_re => dout_4_1_re,  -- sfix16_En15
              dout_4_1_im => dout_4_1_im,  -- sfix16_En15
              dout_4_1_vld => dout_4_1_vld,
              dinXTwdl_4_1_vld => dinXTwdl_4_1_vld
              );

  u_twdlROM_5_1 : TWDLROM_5_1
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_4_1_vld => dout_4_1_vld,
              twdl_5_1_re => twdl_5_1_re,  -- sfix16_En14
              twdl_5_1_im => twdl_5_1_im  -- sfix16_En14
              );

  u_CTRL1_5_1 : RADIX22FFT_CTRL1_5
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dinXTwdl_5_1_vld => dinXTwdl_5_1_vld,
              dinXTwdl_5_1_vld_1 => dinXTwdl_5_1_vld,
              rd_5_Addr => rd_5_Addr,  -- ufix1
              rd_5_Enb => rd_5_Enb,
              proc_5_enb => proc_5_enb
              );

  u_SDF1_5_1 : RADIX22FFT_SDF1_5
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              din_5_1_re_dly => std_logic_vector(din_5_1_re_dly),  -- sfix16_En15
              din_5_1_im_dly => std_logic_vector(din_5_1_im_dly),  -- sfix16_En15
              din_5_vld_dly => din_5_vld_dly,
              rd_5_Addr => rd_5_Addr,  -- ufix1
              rd_5_Enb => rd_5_Enb,
              twdl_5_1_re => twdl_5_1_re,  -- sfix16_En14
              twdl_5_1_im => twdl_5_1_im,  -- sfix16_En14
              proc_5_enb => proc_5_enb,
              dout_5_1_re => dout_5_1_re,  -- sfix16_En15
              dout_5_1_im => dout_5_1_im,  -- sfix16_En15
              dout_5_1_vld => dout_5_1_vld,
              dinXTwdl_5_1_vld => dinXTwdl_5_1_vld
              );

  u_CTRL2_6_1 : RADIX22FFT_CTRL1_6
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_5_1_vld => dout_5_1_vld,
              dinXTwdl_6_1_vld => dinXTwdl_6_1_vld,
              rd_6_Addr => rd_6_Addr,
              rd_6_Enb => rd_6_Enb,
              proc_6_enb => proc_6_enb,
              multiply_6_J => multiply_6_J
              );

  u_SDF2_6_1 : RADIX22FFT_SDF2_6
    PORT MAP( clk => clk,
              reset => reset,
              enb => enb,
              dout_5_1_re => dout_5_1_re,  -- sfix16_En15
              dout_5_1_im => dout_5_1_im,  -- sfix16_En15
              dout_5_1_vld => dout_5_1_vld,
              rd_6_Addr => rd_6_Addr,
              rd_6_Enb => rd_6_Enb,
              proc_6_enb => proc_6_enb,
              multiply_6_J => multiply_6_J,
              dout_6_1_re => dout_6_1_re,  -- sfix16_En15
              dout_6_1_im => dout_6_1_im,  -- sfix16_En15
              dinXTwdl_6_1_vld => dinXTwdl_6_1_vld
              );

  softReset <= '0';

  dataIn_im_signed <= signed(dataIn_im);

  dataIn_re_signed <= signed(dataIn_re);

  intdelay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg(0) <= to_signed(16#0000#, 16);
      intdelay_reg(1) <= to_signed(16#0000#, 16);
      intdelay_reg(2) <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg(0) <= to_signed(16#0000#, 16);
          intdelay_reg(1) <= to_signed(16#0000#, 16);
          intdelay_reg(2) <= to_signed(16#0000#, 16);
        ELSE 
          intdelay_reg(0) <= intdelay_reg_next(0);
          intdelay_reg(1) <= intdelay_reg_next(1);
          intdelay_reg(2) <= intdelay_reg_next(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_process;

  din_1_1_re_dly <= intdelay_reg(2);
  intdelay_reg_next(0) <= dataIn_im_signed;
  intdelay_reg_next(1) <= intdelay_reg(0);
  intdelay_reg_next(2) <= intdelay_reg(1);

  intdelay_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_1(0) <= to_signed(16#0000#, 16);
      intdelay_reg_1(1) <= to_signed(16#0000#, 16);
      intdelay_reg_1(2) <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_1(0) <= to_signed(16#0000#, 16);
          intdelay_reg_1(1) <= to_signed(16#0000#, 16);
          intdelay_reg_1(2) <= to_signed(16#0000#, 16);
        ELSE 
          intdelay_reg_1(0) <= intdelay_reg_next_1(0);
          intdelay_reg_1(1) <= intdelay_reg_next_1(1);
          intdelay_reg_1(2) <= intdelay_reg_next_1(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_1_process;

  din_1_1_im_dly <= intdelay_reg_1(2);
  intdelay_reg_next_1(0) <= dataIn_re_signed;
  intdelay_reg_next_1(1) <= intdelay_reg_1(0);
  intdelay_reg_next_1(2) <= intdelay_reg_1(1);

  intdelay_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_2(0) <= '0';
      intdelay_reg_2(1) <= '0';
      intdelay_reg_2(2) <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_2(0) <= '0';
          intdelay_reg_2(1) <= '0';
          intdelay_reg_2(2) <= '0';
        ELSE 
          intdelay_reg_2(0) <= intdelay_reg_next_2(0);
          intdelay_reg_2(1) <= intdelay_reg_next_2(1);
          intdelay_reg_2(2) <= intdelay_reg_next_2(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_2_process;

  din_1_vld_dly <= intdelay_reg_2(2);
  intdelay_reg_next_2(0) <= validIn;
  intdelay_reg_next_2(1) <= intdelay_reg_2(0);
  intdelay_reg_next_2(2) <= intdelay_reg_2(1);

  dout_2_1_re_signed <= signed(dout_2_1_re);

  dout_2_1_im_signed <= signed(dout_2_1_im);

  intdelay_3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_3(0) <= to_signed(16#0000#, 16);
      intdelay_reg_3(1) <= to_signed(16#0000#, 16);
      intdelay_reg_3(2) <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_3(0) <= to_signed(16#0000#, 16);
          intdelay_reg_3(1) <= to_signed(16#0000#, 16);
          intdelay_reg_3(2) <= to_signed(16#0000#, 16);
        ELSE 
          intdelay_reg_3(0) <= intdelay_reg_next_3(0);
          intdelay_reg_3(1) <= intdelay_reg_next_3(1);
          intdelay_reg_3(2) <= intdelay_reg_next_3(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_3_process;

  din_3_1_re_dly <= intdelay_reg_3(2);
  intdelay_reg_next_3(0) <= dout_2_1_re_signed;
  intdelay_reg_next_3(1) <= intdelay_reg_3(0);
  intdelay_reg_next_3(2) <= intdelay_reg_3(1);

  intdelay_4_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_4(0) <= to_signed(16#0000#, 16);
      intdelay_reg_4(1) <= to_signed(16#0000#, 16);
      intdelay_reg_4(2) <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_4(0) <= to_signed(16#0000#, 16);
          intdelay_reg_4(1) <= to_signed(16#0000#, 16);
          intdelay_reg_4(2) <= to_signed(16#0000#, 16);
        ELSE 
          intdelay_reg_4(0) <= intdelay_reg_next_4(0);
          intdelay_reg_4(1) <= intdelay_reg_next_4(1);
          intdelay_reg_4(2) <= intdelay_reg_next_4(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_4_process;

  din_3_1_im_dly <= intdelay_reg_4(2);
  intdelay_reg_next_4(0) <= dout_2_1_im_signed;
  intdelay_reg_next_4(1) <= intdelay_reg_4(0);
  intdelay_reg_next_4(2) <= intdelay_reg_4(1);

  intdelay_5_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_5(0) <= '0';
      intdelay_reg_5(1) <= '0';
      intdelay_reg_5(2) <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_5(0) <= '0';
          intdelay_reg_5(1) <= '0';
          intdelay_reg_5(2) <= '0';
        ELSE 
          intdelay_reg_5(0) <= intdelay_reg_next_5(0);
          intdelay_reg_5(1) <= intdelay_reg_next_5(1);
          intdelay_reg_5(2) <= intdelay_reg_next_5(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_5_process;

  din_3_vld_dly <= intdelay_reg_5(2);
  intdelay_reg_next_5(0) <= dout_2_1_vld;
  intdelay_reg_next_5(1) <= intdelay_reg_5(0);
  intdelay_reg_next_5(2) <= intdelay_reg_5(1);

  dout_4_1_re_signed <= signed(dout_4_1_re);

  dout_4_1_im_signed <= signed(dout_4_1_im);

  intdelay_6_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_6(0) <= to_signed(16#0000#, 16);
      intdelay_reg_6(1) <= to_signed(16#0000#, 16);
      intdelay_reg_6(2) <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_6(0) <= to_signed(16#0000#, 16);
          intdelay_reg_6(1) <= to_signed(16#0000#, 16);
          intdelay_reg_6(2) <= to_signed(16#0000#, 16);
        ELSE 
          intdelay_reg_6(0) <= intdelay_reg_next_6(0);
          intdelay_reg_6(1) <= intdelay_reg_next_6(1);
          intdelay_reg_6(2) <= intdelay_reg_next_6(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_6_process;

  din_5_1_re_dly <= intdelay_reg_6(2);
  intdelay_reg_next_6(0) <= dout_4_1_re_signed;
  intdelay_reg_next_6(1) <= intdelay_reg_6(0);
  intdelay_reg_next_6(2) <= intdelay_reg_6(1);

  intdelay_7_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_7(0) <= to_signed(16#0000#, 16);
      intdelay_reg_7(1) <= to_signed(16#0000#, 16);
      intdelay_reg_7(2) <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_7(0) <= to_signed(16#0000#, 16);
          intdelay_reg_7(1) <= to_signed(16#0000#, 16);
          intdelay_reg_7(2) <= to_signed(16#0000#, 16);
        ELSE 
          intdelay_reg_7(0) <= intdelay_reg_next_7(0);
          intdelay_reg_7(1) <= intdelay_reg_next_7(1);
          intdelay_reg_7(2) <= intdelay_reg_next_7(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_7_process;

  din_5_1_im_dly <= intdelay_reg_7(2);
  intdelay_reg_next_7(0) <= dout_4_1_im_signed;
  intdelay_reg_next_7(1) <= intdelay_reg_7(0);
  intdelay_reg_next_7(2) <= intdelay_reg_7(1);

  intdelay_8_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      intdelay_reg_8(0) <= '0';
      intdelay_reg_8(1) <= '0';
      intdelay_reg_8(2) <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        IF softReset = '1' THEN
          intdelay_reg_8(0) <= '0';
          intdelay_reg_8(1) <= '0';
          intdelay_reg_8(2) <= '0';
        ELSE 
          intdelay_reg_8(0) <= intdelay_reg_next_8(0);
          intdelay_reg_8(1) <= intdelay_reg_next_8(1);
          intdelay_reg_8(2) <= intdelay_reg_next_8(2);
        END IF;
      END IF;
    END IF;
  END PROCESS intdelay_8_process;

  din_5_vld_dly <= intdelay_reg_8(2);
  intdelay_reg_next_8(0) <= dout_4_1_vld;
  intdelay_reg_next_8(1) <= intdelay_reg_8(0);
  intdelay_reg_next_8(2) <= intdelay_reg_8(1);

  dataOut_re <= dout_6_1_im;

  dataOut_im <= dout_6_1_re;

END rtl;

