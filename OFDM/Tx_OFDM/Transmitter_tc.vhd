-- -------------------------------------------------------------
-- 
-- File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\OFDM\Tx_OFDM\Transmitter_tc.vhd
-- Created: 2023-04-09 15:00:09
-- 
-- Generated by MATLAB 9.14 and HDL Coder 4.1
-- 
-- -------------------------------------------------------------


-- -------------------------------------------------------------
-- 
-- Module: Transmitter_tc
-- Source Path: Transmitter_tc
-- Hierarchy Level: 1
-- 
-- Master clock enable input: clk_enable
-- 
-- enb         : identical to clk_enable
-- enb_1_1_1   : identical to clk_enable
-- enb_1_15_0  : 15x slower than clk with last phase
-- enb_1_15_1  : 15x slower than clk with phase 1
-- enb_1_960_0 : 960x slower than clk with last phase
-- enb_1_960_1 : 960x slower than clk with phase 1
-- 
-- -------------------------------------------------------------
LIBRARY IEEE;
USE IEEE.std_logic_1164.ALL;
USE IEEE.numeric_std.ALL;

ENTITY Transmitter_tc IS
  PORT( clk                               :   IN    std_logic;
        reset                             :   IN    std_logic;
        clk_enable                        :   IN    std_logic;
        enb                               :   OUT   std_logic;
        enb_1_1_1                         :   OUT   std_logic;
        enb_1_15_0                        :   OUT   std_logic;
        enb_1_15_1                        :   OUT   std_logic;
        enb_1_960_0                       :   OUT   std_logic;
        enb_1_960_1                       :   OUT   std_logic
        );
END Transmitter_tc;


ARCHITECTURE rtl OF Transmitter_tc IS

  -- Signals
  SIGNAL count15                          : unsigned(3 DOWNTO 0);  -- ufix4
  SIGNAL comp_0_tmp                       : std_logic;
  SIGNAL phase_0_tmp                      : std_logic;
  SIGNAL phase_0                          : std_logic;
  SIGNAL enb_1_15_0_1                     : std_logic;
  SIGNAL comp_1_tmp                       : std_logic;
  SIGNAL phase_1_tmp                      : std_logic;
  SIGNAL phase_1                          : std_logic;
  SIGNAL enb_1_15_1_1                     : std_logic;
  SIGNAL count960                         : unsigned(9 DOWNTO 0);  -- ufix10
  SIGNAL comp_0_tmp_1                     : std_logic;
  SIGNAL phase_0_tmp_1                    : std_logic;
  SIGNAL phase_0_1                        : std_logic;
  SIGNAL enb_1_960_0_1                    : std_logic;
  SIGNAL comp_1_tmp_1                     : std_logic;
  SIGNAL phase_1_tmp_1                    : std_logic;
  SIGNAL phase_1_1                        : std_logic;
  SIGNAL enb_1_960_1_1                    : std_logic;

BEGIN
  enb <= clk_enable;

  enb_1_1_1 <= clk_enable;

  -- Count limited, Unsigned Counter
  --  initial value   = 1
  --  step value      = 1
  --  count to value  = 14
  counter_15_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      count15 <= to_unsigned(16#1#, 4);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF clk_enable = '1' THEN
        IF count15 >= to_unsigned(16#E#, 4) THEN 
          count15 <= to_unsigned(16#0#, 4);
        ELSE 
          count15 <= count15 + to_unsigned(16#1#, 4);
        END IF;
      END IF;
    END IF;
  END PROCESS counter_15_process;


  
  comp_0_tmp <= '1' WHEN count15 = to_unsigned(16#E#, 4) ELSE
      '0';

  phase_0_tmp <= comp_0_tmp AND clk_enable;

  phase_delay_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      phase_0 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF clk_enable = '1' THEN
        phase_0 <= phase_0_tmp;
      END IF;
    END IF;
  END PROCESS phase_delay_process;


  enb_1_15_0_1 <= phase_0 AND clk_enable;

  enb_1_15_0 <= enb_1_15_0_1;

  
  comp_1_tmp <= '1' WHEN count15 = to_unsigned(16#0#, 4) ELSE
      '0';

  phase_1_tmp <= comp_1_tmp AND clk_enable;

  phase_delay_1_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      phase_1 <= '1';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF clk_enable = '1' THEN
        phase_1 <= phase_1_tmp;
      END IF;
    END IF;
  END PROCESS phase_delay_1_process;


  enb_1_15_1_1 <= phase_1 AND clk_enable;

  enb_1_15_1 <= enb_1_15_1_1;

  -- Count limited, Unsigned Counter
  --  initial value   = 1
  --  step value      = 1
  --  count to value  = 959
  counter_960_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      count960 <= to_unsigned(16#001#, 10);
    ELSIF clk'EVENT AND clk = '1' THEN
      IF clk_enable = '1' THEN
        IF count960 >= to_unsigned(16#3BF#, 10) THEN 
          count960 <= to_unsigned(16#000#, 10);
        ELSE 
          count960 <= count960 + to_unsigned(16#001#, 10);
        END IF;
      END IF;
    END IF;
  END PROCESS counter_960_process;


  
  comp_0_tmp_1 <= '1' WHEN count960 = to_unsigned(16#3BF#, 10) ELSE
      '0';

  phase_0_tmp_1 <= comp_0_tmp_1 AND clk_enable;

  phase_delay_2_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      phase_0_1 <= '0';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF clk_enable = '1' THEN
        phase_0_1 <= phase_0_tmp_1;
      END IF;
    END IF;
  END PROCESS phase_delay_2_process;


  enb_1_960_0_1 <= phase_0_1 AND clk_enable;

  enb_1_960_0 <= enb_1_960_0_1;

  
  comp_1_tmp_1 <= '1' WHEN count960 = to_unsigned(16#000#, 10) ELSE
      '0';

  phase_1_tmp_1 <= comp_1_tmp_1 AND clk_enable;

  phase_delay_3_process : PROCESS (clk, reset)
  BEGIN
    IF reset = '1' THEN
      phase_1_1 <= '1';
    ELSIF clk'EVENT AND clk = '1' THEN
      IF clk_enable = '1' THEN
        phase_1_1 <= phase_1_tmp_1;
      END IF;
    END IF;
  END PROCESS phase_delay_3_process;


  enb_1_960_1_1 <= phase_1_1 AND clk_enable;

  enb_1_960_1 <= enb_1_960_1_1;

END rtl;

