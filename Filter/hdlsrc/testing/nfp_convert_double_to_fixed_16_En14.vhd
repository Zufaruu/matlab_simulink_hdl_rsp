-- -------------------------------------------------------------
-- 
-- File Name: hdlsrc\testing\nfp_convert_double_to_fixed_16_En14.vhd
-- Created: 2023-04-09 11:40:40
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: nfp_convert_double_to_fixed_16_En14
-- Source Path: testing/custom buffer/nfp_convert_double_to_fixed_16_En14
-- Hierarchy Level: 1
-- 
-- {Latency Strategy = "Max"}
-- 
-- {Rounding Mode = Floor}
-- {Overflow Mode = Wrap}
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY nfp_convert_double_to_fixed_16_En14 IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        enb                               :   IN    std_logic;
        nfp_in                            :   IN    std_logic_vector(63 DOWNTO 0);  -- ufix64
        nfp_out                           :   OUT   std_logic_vector(15 DOWNTO 0)  -- sfix16_En14
        );
END nfp_convert_double_to_fixed_16_En14;


ARCHITECTURE rtl OF nfp_convert_double_to_fixed_16_En14 IS

  -- Signals
  SIGNAL nfp_in_unsigned                  : unsigned(63 DOWNTO 0);  -- ufix64
  SIGNAL In1                              : std_logic;  -- ufix1
  SIGNAL In2                              : unsigned(10 DOWNTO 0);  -- ufix11
  SIGNAL In3                              : unsigned(51 DOWNTO 0);  -- ufix52
  SIGNAL Delay8_reg                       : std_logic_vector(0 TO 3);  -- ufix1 [4]
  SIGNAL Delay8_reg_next                  : std_logic_vector(0 TO 3);  -- ufix1 [4]
  SIGNAL Delay8_out1                      : std_logic;  -- ufix1
  SIGNAL Delay3_out1                      : unsigned(10 DOWNTO 0);  -- ufix11
  SIGNAL Constant3_out1                   : unsigned(15 DOWNTO 0);  -- uint16
  SIGNAL Constant2_out1                   : unsigned(15 DOWNTO 0);  -- uint16
  SIGNAL Switch2_out1                     : unsigned(15 DOWNTO 0);  -- uint16
  SIGNAL Add_sub_cast                     : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL Add_sub_cast_1                   : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL alphave                          : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL alphave_1                        : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL alphave_2                        : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL switch_compare_1                 : std_logic;
  SIGNAL Compare_To_Zero_out1             : std_logic;  -- ufix1
  SIGNAL Delay4_out1                      : unsigned(51 DOWNTO 0);  -- ufix52
  SIGNAL Bit_Concat_out1                  : unsigned(52 DOWNTO 0);  -- ufix53
  SIGNAL Data_Type_Conversion1_out1       : unsigned(52 DOWNTO 0);  -- ufix53_En52
  SIGNAL Data_Type_Conversion1_out1_1     : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Data_Type_Conversion1_out1_2     : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Unary_Minus_in0                  : signed(12 DOWNTO 0);  -- sfix13
  SIGNAL alphave_3                        : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL Delay10_out1                     : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic1_zerosig        : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic1_selsig         : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic1_zerosig_1      : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic1_selsig_1       : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic1_cast           : signed(15 DOWNTO 0);  -- int16
  SIGNAL Shift_Arithmetic1_out1           : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Compare_To_Zero1_out1            : std_logic;  -- ufix1
  SIGNAL Compare_To_Zero_out1_1           : std_logic;  -- ufix1
  SIGNAL Logical_Operator_out1            : std_logic;  -- ufix1
  SIGNAL Bit_Set_out1                     : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Switch2_out1_1                   : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Delay4_out1_1                    : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL shift_arithmetic2_zerosig        : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic2_selsig         : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic2_zerosig_1      : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic2_selsig_1       : signed(11 DOWNTO 0);  -- sfix12
  SIGNAL shift_arithmetic2_cast           : signed(15 DOWNTO 0);  -- int16
  SIGNAL Shift_Arithmetic2_out1           : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Delay6_out1                      : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Switch1_out1                     : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Delay2_out1                      : unsigned(54 DOWNTO 0);  -- ufix55_En53
  SIGNAL Data_Type_Conversion3_out1       : signed(55 DOWNTO 0);  -- sfix56_En53
  SIGNAL Unary_Minus_cast                 : signed(56 DOWNTO 0);  -- sfix57_En53
  SIGNAL Unary_Minus_cast_1               : signed(56 DOWNTO 0);  -- sfix57_En53
  SIGNAL Unary_Minus_out1                 : signed(55 DOWNTO 0);  -- sfix56_En53
  SIGNAL Switch1_out1_1                   : signed(55 DOWNTO 0);  -- sfix56_En53
  SIGNAL Delay1_out1                      : signed(55 DOWNTO 0);  -- sfix56_En53
  SIGNAL Data_Type_Conversion2_out1       : signed(15 DOWNTO 0);  -- sfix16_En14
  SIGNAL Delay8_out1_1                    : signed(15 DOWNTO 0);  -- sfix16_En14

BEGIN
  nfp_in_unsigned <= unsigned(nfp_in);

  -- Split 64 bit word into FP sign, exponent, mantissa
  In1 <= nfp_in_unsigned(63);
  In2 <= nfp_in_unsigned(62 DOWNTO 52);
  In3 <= nfp_in_unsigned(51 DOWNTO 0);

  Delay8_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay8_reg(0) <= '0';
      Delay8_reg(1) <= '0';
      Delay8_reg(2) <= '0';
      Delay8_reg(3) <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay8_reg(0) <= Delay8_reg_next(0);
        Delay8_reg(1) <= Delay8_reg_next(1);
        Delay8_reg(2) <= Delay8_reg_next(2);
        Delay8_reg(3) <= Delay8_reg_next(3);
      END IF;
    END IF;
  END PROCESS Delay8_process;

  Delay8_out1 <= Delay8_reg(3);
  Delay8_reg_next(0) <= In1;
  Delay8_reg_next(1) <= Delay8_reg(0);
  Delay8_reg_next(2) <= Delay8_reg(1);
  Delay8_reg_next(3) <= Delay8_reg(2);

  Delay3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay3_out1 <= to_unsigned(16#000#, 11);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay3_out1 <= In2;
      END IF;
    END IF;
  END PROCESS Delay3_process;


  Constant3_out1 <= to_unsigned(16#03FE#, 16);

  Constant2_out1 <= to_unsigned(16#03FF#, 16);

  
  Switch2_out1 <= Constant3_out1 WHEN Delay3_out1 = to_unsigned(16#000#, 11) ELSE
      Constant2_out1;

  Add_sub_cast <= signed(resize(Delay3_out1, 12));
  Add_sub_cast_1 <= signed(Switch2_out1(11 DOWNTO 0));
  alphave <= Add_sub_cast - Add_sub_cast_1;

  reduced_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      alphave_1 <= to_signed(16#000#, 12);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        alphave_1 <= alphave;
      END IF;
    END IF;
  END PROCESS reduced_process;


  reduced_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      alphave_2 <= to_signed(16#000#, 12);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        alphave_2 <= alphave_1;
      END IF;
    END IF;
  END PROCESS reduced_1_process;


  
  switch_compare_1 <= '1' WHEN alphave_2 >= to_signed(16#000#, 12) ELSE
      '0';

  
  Compare_To_Zero_out1 <= '1' WHEN Delay3_out1 /= to_unsigned(16#000#, 11) ELSE
      '0';

  Delay4_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay4_out1 <= to_unsigned(0, 52);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay4_out1 <= In3;
      END IF;
    END IF;
  END PROCESS Delay4_process;


  Bit_Concat_out1 <= Compare_To_Zero_out1 & Delay4_out1;

  Data_Type_Conversion1_out1 <= Bit_Concat_out1;

  Data_Type_Conversion1_out1_1 <= resize(Data_Type_Conversion1_out1 & '0', 55);

  reduced_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Data_Type_Conversion1_out1_2 <= to_unsigned(0, 55);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Data_Type_Conversion1_out1_2 <= Data_Type_Conversion1_out1_1;
      END IF;
    END IF;
  END PROCESS reduced_2_process;


  Unary_Minus_in0 <=  - (resize(alphave, 13));
  alphave_3 <= Unary_Minus_in0(11 DOWNTO 0);

  Delay10_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay10_out1 <= to_signed(16#000#, 12);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay10_out1 <= alphave_3;
      END IF;
    END IF;
  END PROCESS Delay10_process;


  shift_arithmetic1_zerosig <= to_signed(16#000#, 12);

  
  shift_arithmetic1_selsig <= Delay10_out1 WHEN Delay10_out1 >= shift_arithmetic1_zerosig ELSE
      shift_arithmetic1_zerosig;

  shift_arithmetic1_zerosig_1 <= to_signed(16#000#, 12);

  
  shift_arithmetic1_selsig_1 <= shift_arithmetic1_selsig WHEN shift_arithmetic1_selsig >= shift_arithmetic1_zerosig_1 ELSE
      shift_arithmetic1_zerosig_1;

  shift_arithmetic1_cast <= resize(shift_arithmetic1_selsig_1, 16);
  Shift_Arithmetic1_out1 <= SHIFT_RIGHT(Data_Type_Conversion1_out1_2, to_integer(shift_arithmetic1_cast));

  
  Compare_To_Zero1_out1 <= '1' WHEN Shift_Arithmetic1_out1 = to_unsigned(0, 55) ELSE
      '0';

  
  Compare_To_Zero_out1_1 <= '1' WHEN Data_Type_Conversion1_out1_2 /= to_unsigned(0, 55) ELSE
      '0';

  Logical_Operator_out1 <= Compare_To_Zero1_out1 AND Compare_To_Zero_out1_1;

  Bit_Set_out1 <= Shift_Arithmetic1_out1 OR to_unsigned(1, 55);

  
  Switch2_out1_1 <= Shift_Arithmetic1_out1 WHEN Logical_Operator_out1 = '0' ELSE
      Bit_Set_out1;

  Delay4_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay4_out1_1 <= to_unsigned(0, 55);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay4_out1_1 <= Switch2_out1_1;
      END IF;
    END IF;
  END PROCESS Delay4_1_process;


  shift_arithmetic2_zerosig <= to_signed(16#000#, 12);

  
  shift_arithmetic2_selsig <= alphave_1 WHEN alphave_1 >= shift_arithmetic2_zerosig ELSE
      shift_arithmetic2_zerosig;

  shift_arithmetic2_zerosig_1 <= to_signed(16#000#, 12);

  
  shift_arithmetic2_selsig_1 <= shift_arithmetic2_selsig WHEN shift_arithmetic2_selsig >= shift_arithmetic2_zerosig_1 ELSE
      shift_arithmetic2_zerosig_1;

  shift_arithmetic2_cast <= resize(shift_arithmetic2_selsig_1, 16);
  Shift_Arithmetic2_out1 <= Data_Type_Conversion1_out1_2 sll to_integer(shift_arithmetic2_cast);

  Delay6_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay6_out1 <= to_unsigned(0, 55);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay6_out1 <= Shift_Arithmetic2_out1;
      END IF;
    END IF;
  END PROCESS Delay6_process;


  
  Switch1_out1 <= Delay4_out1_1 WHEN switch_compare_1 = '0' ELSE
      Delay6_out1;

  Delay2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay2_out1 <= to_unsigned(0, 55);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay2_out1 <= Switch1_out1;
      END IF;
    END IF;
  END PROCESS Delay2_process;


  Data_Type_Conversion3_out1 <= signed(resize(Delay2_out1, 56));

  Unary_Minus_cast <= resize(Data_Type_Conversion3_out1, 57);
  Unary_Minus_cast_1 <=  - (Unary_Minus_cast);
  Unary_Minus_out1 <= Unary_Minus_cast_1(55 DOWNTO 0);

  
  Switch1_out1_1 <= Data_Type_Conversion3_out1 WHEN Delay8_out1 = '0' ELSE
      Unary_Minus_out1;

  Delay1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay1_out1 <= to_signed(0, 56);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay1_out1 <= Switch1_out1_1;
      END IF;
    END IF;
  END PROCESS Delay1_process;


  Data_Type_Conversion2_out1 <= Delay1_out1(54 DOWNTO 39);

  Delay8_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      Delay8_out1_1 <= to_signed(16#0000#, 16);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF enb = '1' THEN
        Delay8_out1_1 <= Data_Type_Conversion2_out1;
      END IF;
    END IF;
  END PROCESS Delay8_1_process;


  nfp_out <= std_logic_vector(Delay8_out1_1);

END rtl;

