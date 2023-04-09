-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\TWDLROM_5_15.vhd
-- Created: 2023-04-09 15:00:12
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: TWDLROM_5_15
-- Source Path: Tx_OFDM/Transmitter/IFFT/TWDLROM_5_15
-- Hierarchy Level: 2
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;
USE work.Transmitter_pkg.ALL;

ENTITY TWDLROM_5_15 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb_1_960_0                       :   IN    std_logic;
        dout_4_vld                        :   IN    std_logic;
        twdl_5_15_re                      :   OUT   std_logic_vector(15 DOWNTO 0);  -- sfix16_En14
        twdl_5_15_im                      :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END TWDLROM_5_15;


ARCHITECTURE rtl OF TWDLROM_5_15 IS

  -- Constants
  CONSTANT Twiddle_re_table_data          : vector_of_signed16(0 TO 31) := 
    (to_signed(16#4000#, 16), to_signed(16#3FFB#, 16), to_signed(16#3FEC#, 16), to_signed(16#3FD4#, 16),
     to_signed(16#3FB1#, 16), to_signed(16#3F85#, 16), to_signed(16#3F4F#, 16), to_signed(16#3F0F#, 16),
     to_signed(16#3EC5#, 16), to_signed(16#3E72#, 16), to_signed(16#3E15#, 16), to_signed(16#3DAF#, 16),
     to_signed(16#3D3F#, 16), to_signed(16#3CC5#, 16), to_signed(16#3C42#, 16), to_signed(16#3BB6#, 16),
     to_signed(16#3B21#, 16), to_signed(16#3A82#, 16), to_signed(16#39DB#, 16), to_signed(16#392B#, 16),
     to_signed(16#3871#, 16), to_signed(16#37B0#, 16), to_signed(16#36E5#, 16), to_signed(16#3612#, 16),
     to_signed(16#3537#, 16), to_signed(16#3453#, 16), to_signed(16#3368#, 16), to_signed(16#3274#, 16),
     to_signed(16#3179#, 16), to_signed(16#3076#, 16), to_signed(16#2F6C#, 16), to_signed(16#2E5A#, 16));  -- sfix16 [32]
  CONSTANT Twiddle_im_table_data          : vector_of_signed16(0 TO 31) := 
    (to_signed(16#0000#, 16), to_signed(-16#0192#, 16), to_signed(-16#0324#, 16), to_signed(-16#04B5#, 16),
     to_signed(-16#0646#, 16), to_signed(-16#07D6#, 16), to_signed(-16#0964#, 16), to_signed(-16#0AF1#, 16),
     to_signed(-16#0C7C#, 16), to_signed(-16#0E06#, 16), to_signed(-16#0F8D#, 16), to_signed(-16#1112#, 16),
     to_signed(-16#1294#, 16), to_signed(-16#1413#, 16), to_signed(-16#1590#, 16), to_signed(-16#1709#, 16),
     to_signed(-16#187E#, 16), to_signed(-16#19EF#, 16), to_signed(-16#1B5D#, 16), to_signed(-16#1CC6#, 16),
     to_signed(-16#1E2B#, 16), to_signed(-16#1F8C#, 16), to_signed(-16#20E7#, 16), to_signed(-16#223D#, 16),
     to_signed(-16#238E#, 16), to_signed(-16#24DA#, 16), to_signed(-16#2620#, 16), to_signed(-16#2760#, 16),
     to_signed(-16#289A#, 16), to_signed(-16#29CE#, 16), to_signed(-16#2AFB#, 16), to_signed(-16#2C21#, 16));  -- sfix16 [32]

  -- Signals
  SIGNAL Radix22TwdlMapping_cnt           : unsigned(5 DOWNTO 0);  -- ufix6
  SIGNAL Radix22TwdlMapping_phase         : unsigned(1 DOWNTO 0);  -- ufix2
  SIGNAL Radix22TwdlMapping_octantReg1    : unsigned(2 DOWNTO 0);  -- ufix3
  SIGNAL Radix22TwdlMapping_twdlAddr_raw  : unsigned(7 DOWNTO 0);  -- ufix8
  SIGNAL Radix22TwdlMapping_twdlAddrMap   : unsigned(4 DOWNTO 0);  -- ufix5
  SIGNAL Radix22TwdlMapping_twdl45Reg     : std_logic;
  SIGNAL Radix22TwdlMapping_dvldReg1      : std_logic;
  SIGNAL Radix22TwdlMapping_dvldReg2      : std_logic;
  SIGNAL Radix22TwdlMapping_cnt_next      : unsigned(5 DOWNTO 0);  -- ufix6
  SIGNAL Radix22TwdlMapping_phase_next    : unsigned(1 DOWNTO 0);  -- ufix2
  SIGNAL Radix22TwdlMapping_octantReg1_next : unsigned(2 DOWNTO 0);  -- ufix3
  SIGNAL Radix22TwdlMapping_twdlAddr_raw_next : unsigned(7 DOWNTO 0);  -- ufix8
  SIGNAL Radix22TwdlMapping_twdlAddrMap_next : unsigned(4 DOWNTO 0);  -- ufix5
  SIGNAL Radix22TwdlMapping_twdl45Reg_next : std_logic;
  SIGNAL Radix22TwdlMapping_dvldReg1_next : std_logic;
  SIGNAL Radix22TwdlMapping_dvldReg2_next : std_logic;
  SIGNAL twdlAddr                         : unsigned(4 DOWNTO 0);  -- ufix5
  SIGNAL twdlAddrVld                      : std_logic;
  SIGNAL twdlOctant                       : unsigned(2 DOWNTO 0);  -- ufix3
  SIGNAL twdl45                           : std_logic;
  SIGNAL twiddleS_re                      : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twiddleReg_re                    : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twiddleS_im                      : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twiddleReg_im                    : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdlOctantReg                    : unsigned(2 DOWNTO 0);  -- ufix3
  SIGNAL twdl45Reg                        : std_logic;
  SIGNAL twdl_5_15_re_tmp                 : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL twdl_5_15_im_tmp                 : signed(15 DOWNTO 0);  -- sfix16_En14

BEGIN
  -- Radix22TwdlMapping
  Radix22TwdlMapping_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Radix22TwdlMapping_octantReg1 <= to_unsigned(16#0#, 3);
      Radix22TwdlMapping_twdlAddr_raw <= to_unsigned(16#00#, 8);
      Radix22TwdlMapping_twdlAddrMap <= to_unsigned(16#00#, 5);
      Radix22TwdlMapping_twdl45Reg <= '0';
      Radix22TwdlMapping_dvldReg1 <= '0';
      Radix22TwdlMapping_dvldReg2 <= '0';
      Radix22TwdlMapping_cnt <= to_unsigned(16#38#, 6);
      Radix22TwdlMapping_phase <= to_unsigned(16#0#, 2);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        Radix22TwdlMapping_cnt <= Radix22TwdlMapping_cnt_next;
        Radix22TwdlMapping_phase <= Radix22TwdlMapping_phase_next;
        Radix22TwdlMapping_octantReg1 <= Radix22TwdlMapping_octantReg1_next;
        Radix22TwdlMapping_twdlAddr_raw <= Radix22TwdlMapping_twdlAddr_raw_next;
        Radix22TwdlMapping_twdlAddrMap <= Radix22TwdlMapping_twdlAddrMap_next;
        Radix22TwdlMapping_twdl45Reg <= Radix22TwdlMapping_twdl45Reg_next;
        Radix22TwdlMapping_dvldReg1 <= Radix22TwdlMapping_dvldReg1_next;
        Radix22TwdlMapping_dvldReg2 <= Radix22TwdlMapping_dvldReg2_next;
      END IF;
    END IF;
  END PROCESS Radix22TwdlMapping_process;

  Radix22TwdlMapping_output : PROCESS (Radix22TwdlMapping_cnt, Radix22TwdlMapping_dvldReg1,
       Radix22TwdlMapping_dvldReg2, Radix22TwdlMapping_octantReg1,
       Radix22TwdlMapping_phase, Radix22TwdlMapping_twdl45Reg,
       Radix22TwdlMapping_twdlAddrMap, Radix22TwdlMapping_twdlAddr_raw,
       dout_4_vld)
    VARIABLE octant : unsigned(2 DOWNTO 0);
    VARIABLE cnt_cast : unsigned(7 DOWNTO 0);
    VARIABLE sub_cast : signed(17 DOWNTO 0);
    VARIABLE sub_temp : signed(17 DOWNTO 0);
    VARIABLE sub_cast_0 : signed(9 DOWNTO 0);
    VARIABLE sub_temp_0 : signed(9 DOWNTO 0);
    VARIABLE sub_cast_1 : signed(9 DOWNTO 0);
    VARIABLE sub_temp_1 : signed(9 DOWNTO 0);
    VARIABLE sub_cast_2 : signed(17 DOWNTO 0);
    VARIABLE sub_temp_2 : signed(17 DOWNTO 0);
    VARIABLE sub_cast_3 : signed(17 DOWNTO 0);
    VARIABLE sub_temp_3 : signed(17 DOWNTO 0);
  BEGIN
    sub_temp := to_signed(16#00000#, 18);
    sub_temp_0 := to_signed(16#000#, 10);
    sub_temp_1 := to_signed(16#000#, 10);
    sub_temp_2 := to_signed(16#00000#, 18);
    sub_temp_3 := to_signed(16#00000#, 18);
    sub_cast_3 := to_signed(16#00000#, 18);
    cnt_cast := to_unsigned(16#00#, 8);
    sub_cast_2 := to_signed(16#00000#, 18);
    sub_cast_1 := to_signed(16#000#, 10);
    sub_cast_0 := to_signed(16#000#, 10);
    sub_cast := to_signed(16#00000#, 18);
    Radix22TwdlMapping_dvldReg2_next <= Radix22TwdlMapping_dvldReg1;
    Radix22TwdlMapping_dvldReg1_next <= dout_4_vld;
    CASE Radix22TwdlMapping_twdlAddr_raw IS
      WHEN "00100000" =>
        octant := to_unsigned(16#0#, 3);
        Radix22TwdlMapping_twdl45Reg_next <= '1';
      WHEN "01000000" =>
        octant := to_unsigned(16#1#, 3);
        Radix22TwdlMapping_twdl45Reg_next <= '0';
      WHEN "01100000" =>
        octant := to_unsigned(16#2#, 3);
        Radix22TwdlMapping_twdl45Reg_next <= '1';
      WHEN "10000000" =>
        octant := to_unsigned(16#3#, 3);
        Radix22TwdlMapping_twdl45Reg_next <= '0';
      WHEN "10100000" =>
        octant := to_unsigned(16#4#, 3);
        Radix22TwdlMapping_twdl45Reg_next <= '1';
      WHEN OTHERS => 
        octant := Radix22TwdlMapping_twdlAddr_raw(7 DOWNTO 5);
        Radix22TwdlMapping_twdl45Reg_next <= '0';
    END CASE;
    Radix22TwdlMapping_octantReg1_next <= octant;
    CASE octant IS
      WHEN "000" =>
        Radix22TwdlMapping_twdlAddrMap_next <= Radix22TwdlMapping_twdlAddr_raw(4 DOWNTO 0);
      WHEN "001" =>
        sub_cast_0 := signed(resize(Radix22TwdlMapping_twdlAddr_raw, 10));
        sub_temp_0 := to_signed(16#040#, 10) - sub_cast_0;
        Radix22TwdlMapping_twdlAddrMap_next <= unsigned(sub_temp_0(4 DOWNTO 0));
      WHEN "010" =>
        sub_cast_1 := signed(resize(Radix22TwdlMapping_twdlAddr_raw, 10));
        sub_temp_1 := sub_cast_1 - to_signed(16#040#, 10);
        Radix22TwdlMapping_twdlAddrMap_next <= unsigned(sub_temp_1(4 DOWNTO 0));
      WHEN "011" =>
        sub_cast_2 := signed(resize(Radix22TwdlMapping_twdlAddr_raw & '0' & '0' & '0' & '0' & '0', 18));
        sub_temp_2 := to_signed(16#01000#, 18) - sub_cast_2;
        Radix22TwdlMapping_twdlAddrMap_next <= unsigned(sub_temp_2(9 DOWNTO 5));
      WHEN "100" =>
        sub_cast_3 := signed(resize(Radix22TwdlMapping_twdlAddr_raw & '0' & '0' & '0' & '0' & '0', 18));
        sub_temp_3 := sub_cast_3 - to_signed(16#01000#, 18);
        Radix22TwdlMapping_twdlAddrMap_next <= unsigned(sub_temp_3(9 DOWNTO 5));
      WHEN OTHERS => 
        sub_cast := signed(resize(Radix22TwdlMapping_twdlAddr_raw & '0' & '0' & '0' & '0' & '0', 18));
        sub_temp := to_signed(16#01800#, 18) - sub_cast;
        Radix22TwdlMapping_twdlAddrMap_next <= unsigned(sub_temp(9 DOWNTO 5));
    END CASE;
    IF Radix22TwdlMapping_phase = to_unsigned(16#0#, 2) THEN 
      Radix22TwdlMapping_twdlAddr_raw_next <= to_unsigned(16#00#, 8);
    ELSIF Radix22TwdlMapping_phase = to_unsigned(16#1#, 2) THEN 
      Radix22TwdlMapping_twdlAddr_raw_next <= resize(Radix22TwdlMapping_cnt, 8) sll 1;
    ELSIF Radix22TwdlMapping_phase = to_unsigned(16#2#, 2) THEN 
      Radix22TwdlMapping_twdlAddr_raw_next <= resize(Radix22TwdlMapping_cnt, 8);
    ELSE 
      cnt_cast := resize(Radix22TwdlMapping_cnt, 8);
      Radix22TwdlMapping_twdlAddr_raw_next <= (cnt_cast sll 1) + cnt_cast;
    END IF;
    Radix22TwdlMapping_phase_next <= to_unsigned(16#0#, 2);
    Radix22TwdlMapping_cnt_next <= to_unsigned(16#38#, 6);
    twdlAddr <= Radix22TwdlMapping_twdlAddrMap;
    twdlAddrVld <= Radix22TwdlMapping_dvldReg2;
    twdlOctant <= Radix22TwdlMapping_octantReg1;
    twdl45 <= Radix22TwdlMapping_twdl45Reg;
  END PROCESS Radix22TwdlMapping_output;


  -- Twiddle ROM1
  twiddleS_re <= Twiddle_re_table_data(to_integer(twdlAddr));

  TWIDDLEROM_RE_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twiddleReg_re <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twiddleReg_re <= twiddleS_re;
      END IF;
    END IF;
  END PROCESS TWIDDLEROM_RE_process;


  -- Twiddle ROM2
  twiddleS_im <= Twiddle_im_table_data(to_integer(twdlAddr));

  TWIDDLEROM_IM_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twiddleReg_im <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twiddleReg_im <= twiddleS_im;
      END IF;
    END IF;
  END PROCESS TWIDDLEROM_IM_process;


  intdelay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twdlOctantReg <= to_unsigned(16#0#, 3);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twdlOctantReg <= twdlOctant;
      END IF;
    END IF;
  END PROCESS intdelay_process;


  intdelay_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      twdl45Reg <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb_1_960_0 = '1' THEN
        twdl45Reg <= twdl45;
      END IF;
    END IF;
  END PROCESS intdelay_1_process;


  -- Radix22TwdlOctCorr
  Radix22TwdlOctCorr_output : PROCESS (twdl45Reg, twdlOctantReg, twiddleReg_im, twiddleReg_re)
    VARIABLE twdlIn_re : signed(15 DOWNTO 0);
    VARIABLE twdlIn_im : signed(15 DOWNTO 0);
    VARIABLE cast : signed(16 DOWNTO 0);
    VARIABLE cast_0 : signed(16 DOWNTO 0);
    VARIABLE cast_1 : signed(16 DOWNTO 0);
    VARIABLE cast_2 : signed(16 DOWNTO 0);
    VARIABLE cast_3 : signed(16 DOWNTO 0);
    VARIABLE cast_4 : signed(16 DOWNTO 0);
    VARIABLE cast_5 : signed(16 DOWNTO 0);
    VARIABLE cast_6 : signed(16 DOWNTO 0);
    VARIABLE cast_7 : signed(16 DOWNTO 0);
    VARIABLE cast_8 : signed(16 DOWNTO 0);
    VARIABLE cast_9 : signed(16 DOWNTO 0);
    VARIABLE cast_10 : signed(16 DOWNTO 0);
  BEGIN
    cast_0 := to_signed(16#00000#, 17);
    cast_2 := to_signed(16#00000#, 17);
    cast_4 := to_signed(16#00000#, 17);
    cast_6 := to_signed(16#00000#, 17);
    cast_8 := to_signed(16#00000#, 17);
    cast_10 := to_signed(16#00000#, 17);
    cast_3 := to_signed(16#00000#, 17);
    cast_9 := to_signed(16#00000#, 17);
    cast_1 := to_signed(16#00000#, 17);
    cast_7 := to_signed(16#00000#, 17);
    cast := to_signed(16#00000#, 17);
    cast_5 := to_signed(16#00000#, 17);
    twdlIn_re := twiddleReg_re;
    twdlIn_im := twiddleReg_im;
    IF twdl45Reg = '1' THEN 
      CASE twdlOctantReg IS
        WHEN "000" =>
          twdlIn_re := to_signed(16#2D41#, 16);
          twdlIn_im := to_signed(-16#2D41#, 16);
        WHEN "010" =>
          twdlIn_re := to_signed(-16#2D41#, 16);
          twdlIn_im := to_signed(-16#2D41#, 16);
        WHEN "100" =>
          twdlIn_re := to_signed(-16#2D41#, 16);
          twdlIn_im := to_signed(16#2D41#, 16);
        WHEN OTHERS => 
          twdlIn_re := to_signed(16#2D41#, 16);
          twdlIn_im := to_signed(-16#2D41#, 16);
      END CASE;
    ELSE 
      CASE twdlOctantReg IS
        WHEN "000" =>
          NULL;
        WHEN "001" =>
          cast := resize(twiddleReg_im, 17);
          cast_0 :=  - (cast);
          twdlIn_re := cast_0(15 DOWNTO 0);
          cast_5 := resize(twiddleReg_re, 17);
          cast_6 :=  - (cast_5);
          twdlIn_im := cast_6(15 DOWNTO 0);
        WHEN "010" =>
          twdlIn_re := twiddleReg_im;
          cast_7 := resize(twiddleReg_re, 17);
          cast_8 :=  - (cast_7);
          twdlIn_im := cast_8(15 DOWNTO 0);
        WHEN "011" =>
          cast_1 := resize(twiddleReg_re, 17);
          cast_2 :=  - (cast_1);
          twdlIn_re := cast_2(15 DOWNTO 0);
          twdlIn_im := twiddleReg_im;
        WHEN "100" =>
          cast_3 := resize(twiddleReg_re, 17);
          cast_4 :=  - (cast_3);
          twdlIn_re := cast_4(15 DOWNTO 0);
          cast_9 := resize(twiddleReg_im, 17);
          cast_10 :=  - (cast_9);
          twdlIn_im := cast_10(15 DOWNTO 0);
        WHEN OTHERS => 
          twdlIn_re := twiddleReg_im;
          twdlIn_im := twiddleReg_re;
      END CASE;
    END IF;
    twdl_5_15_re_tmp <= twdlIn_re;
    twdl_5_15_im_tmp <= twdlIn_im;
  END PROCESS Radix22TwdlOctCorr_output;


  twdl_5_15_re <= std_logic_vector(twdl_5_15_re_tmp);

  twdl_5_15_im <= std_logic_vector(twdl_5_15_im_tmp);

END rtl;

