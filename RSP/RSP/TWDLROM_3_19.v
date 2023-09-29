// -------------------------------------------------------------
// 
// File Name: C:\Users\acer\OneDrive\Documents\ITS CAK V2\Magang\BRIN\Progress\Matlab Simulink\HDL Coder\proyek\RSP\RSP\RSP\TWDLROM_3_19.v
// Created: 2023-09-29 06:09:50
// 
// Generated by MATLAB 9.14 and HDL Coder 4.1
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: TWDLROM_3_19
// Source Path: RSP/RSP/Velocity Processor/Doppler Processing/Doppler Processor/FFT/TWDLROM_3_19
// Hierarchy Level: 5
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module TWDLROM_3_19
          (clk,
           reset,
           enb,
           dout_2_vld,
           twdl_3_19_re,
           twdl_3_19_im);


  input   clk;
  input   reset;
  input   enb;
  input   dout_2_vld;
  output  signed [11:0] twdl_3_19_re;  // sfix12_En10
  output  signed [11:0] twdl_3_19_im;  // sfix12_En10


  reg [2:0] Radix22TwdlMapping_cnt;  // ufix3
  reg [1:0] Radix22TwdlMapping_phase;  // ufix2
  reg [2:0] Radix22TwdlMapping_octantReg1;  // ufix3
  reg [4:0] Radix22TwdlMapping_twdlAddr_raw;  // ufix5
  reg [1:0] Radix22TwdlMapping_twdlAddrMap;  // ufix2
  reg  Radix22TwdlMapping_twdl45Reg;
  reg  Radix22TwdlMapping_dvldReg1;
  reg  Radix22TwdlMapping_dvldReg2;
  reg [2:0] Radix22TwdlMapping_cnt_next;  // ufix3
  reg [1:0] Radix22TwdlMapping_phase_next;  // ufix2
  reg [2:0] Radix22TwdlMapping_octantReg1_next;  // ufix3
  reg [4:0] Radix22TwdlMapping_twdlAddr_raw_next;  // ufix5
  reg [1:0] Radix22TwdlMapping_twdlAddrMap_next;  // ufix2
  reg  Radix22TwdlMapping_twdl45Reg_next;
  reg  Radix22TwdlMapping_dvldReg1_next;
  reg  Radix22TwdlMapping_dvldReg2_next;
  reg [1:0] twdlAddr;  // ufix2
  reg  twdlAddrVld;
  reg [2:0] twdlOctant;  // ufix3
  reg  twdl45;
  wire signed [11:0] Twiddle_re_table_data [0:3];  // sfix12_En10 [4]
  wire signed [11:0] twiddleS_re;  // sfix12_En10
  reg signed [11:0] twiddleReg_re;  // sfix12_En10
  wire signed [11:0] Twiddle_im_table_data [0:3];  // sfix12_En10 [4]
  wire signed [11:0] twiddleS_im;  // sfix12_En10
  reg signed [11:0] twiddleReg_im;  // sfix12_En10
  reg [2:0] twdlOctantReg;  // ufix3
  reg  twdl45Reg;
  reg signed [11:0] twdl_3_19_re_1;  // sfix12_En10
  reg signed [11:0] twdl_3_19_im_1;  // sfix12_En10
  reg [2:0] Radix22TwdlMapping_octant;  // ufix3
  reg [4:0] Radix22TwdlMapping_cnt_cast;  // ufix5
  reg signed [11:0] Radix22TwdlMapping_sub_cast;  // sfix12_En2
  reg signed [11:0] Radix22TwdlMapping_sub_temp;  // sfix12_En2
  reg signed [6:0] Radix22TwdlMapping_sub_temp_0;  // sfix7
  reg signed [6:0] Radix22TwdlMapping_sub_temp_1;  // sfix7
  reg signed [11:0] Radix22TwdlMapping_sub_cast_0;  // sfix12_En2
  reg signed [11:0] Radix22TwdlMapping_sub_temp_2;  // sfix12_En2
  reg signed [11:0] Radix22TwdlMapping_sub_cast_1;  // sfix12_En2
  reg signed [11:0] Radix22TwdlMapping_sub_temp_3;  // sfix12_En2
  reg [4:0] Radix22TwdlMapping_t_0_0;  // ufix5
  reg signed [6:0] Radix22TwdlMapping_t_1;  // sfix7
  reg signed [6:0] Radix22TwdlMapping_t_2_0;  // sfix7
  reg signed [11:0] Radix22TwdlOctCorr_twdlIn_re;  // sfix12_En10
  reg signed [11:0] Radix22TwdlOctCorr_twdlIn_im;  // sfix12_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_0;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_1;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_2;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_3;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_4;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_5;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_6;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_7;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_8;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_9;  // sfix13_En10
  reg signed [12:0] Radix22TwdlOctCorr_cast_10;  // sfix13_En10


  // Radix22TwdlMapping
  always @(posedge clk or posedge reset)
    begin : Radix22TwdlMapping_process
      if (reset == 1'b1) begin
        Radix22TwdlMapping_octantReg1 <= 3'b000;
        Radix22TwdlMapping_twdlAddr_raw <= 5'b00000;
        Radix22TwdlMapping_twdlAddrMap <= 2'b00;
        Radix22TwdlMapping_twdl45Reg <= 1'b0;
        Radix22TwdlMapping_dvldReg1 <= 1'b0;
        Radix22TwdlMapping_dvldReg2 <= 1'b0;
        Radix22TwdlMapping_cnt <= 3'b010;
        Radix22TwdlMapping_phase <= 2'b10;
      end
      else begin
        if (enb) begin
          Radix22TwdlMapping_cnt <= Radix22TwdlMapping_cnt_next;
          Radix22TwdlMapping_phase <= Radix22TwdlMapping_phase_next;
          Radix22TwdlMapping_octantReg1 <= Radix22TwdlMapping_octantReg1_next;
          Radix22TwdlMapping_twdlAddr_raw <= Radix22TwdlMapping_twdlAddr_raw_next;
          Radix22TwdlMapping_twdlAddrMap <= Radix22TwdlMapping_twdlAddrMap_next;
          Radix22TwdlMapping_twdl45Reg <= Radix22TwdlMapping_twdl45Reg_next;
          Radix22TwdlMapping_dvldReg1 <= Radix22TwdlMapping_dvldReg1_next;
          Radix22TwdlMapping_dvldReg2 <= Radix22TwdlMapping_dvldReg2_next;
        end
      end
    end

  always @(Radix22TwdlMapping_cnt, Radix22TwdlMapping_dvldReg1,
       Radix22TwdlMapping_dvldReg2, Radix22TwdlMapping_octantReg1,
       Radix22TwdlMapping_phase, Radix22TwdlMapping_twdl45Reg,
       Radix22TwdlMapping_twdlAddrMap, Radix22TwdlMapping_twdlAddr_raw,
       dout_2_vld) begin
    Radix22TwdlMapping_sub_temp = 12'sb000000000000;
    Radix22TwdlMapping_sub_temp_0 = 7'sb0000000;
    Radix22TwdlMapping_sub_temp_1 = 7'sb0000000;
    Radix22TwdlMapping_sub_temp_2 = 12'sb000000000000;
    Radix22TwdlMapping_sub_temp_3 = 12'sb000000000000;
    Radix22TwdlMapping_sub_cast_1 = 12'sb000000000000;
    Radix22TwdlMapping_t_0_0 = 5'b00000;
    Radix22TwdlMapping_cnt_cast = 5'b00000;
    Radix22TwdlMapping_sub_cast_0 = 12'sb000000000000;
    Radix22TwdlMapping_t_2_0 = 7'sb0000000;
    Radix22TwdlMapping_t_1 = 7'sb0000000;
    Radix22TwdlMapping_sub_cast = 12'sb000000000000;
    Radix22TwdlMapping_dvldReg2_next = Radix22TwdlMapping_dvldReg1;
    Radix22TwdlMapping_dvldReg1_next = dout_2_vld;
    case ( Radix22TwdlMapping_twdlAddr_raw)
      5'b00100 :
        begin
          Radix22TwdlMapping_octant = 3'b000;
          Radix22TwdlMapping_twdl45Reg_next = 1'b1;
        end
      5'b01000 :
        begin
          Radix22TwdlMapping_octant = 3'b001;
          Radix22TwdlMapping_twdl45Reg_next = 1'b0;
        end
      5'b01100 :
        begin
          Radix22TwdlMapping_octant = 3'b010;
          Radix22TwdlMapping_twdl45Reg_next = 1'b1;
        end
      5'b10000 :
        begin
          Radix22TwdlMapping_octant = 3'b011;
          Radix22TwdlMapping_twdl45Reg_next = 1'b0;
        end
      5'b10100 :
        begin
          Radix22TwdlMapping_octant = 3'b100;
          Radix22TwdlMapping_twdl45Reg_next = 1'b1;
        end
      default :
        begin
          Radix22TwdlMapping_octant = Radix22TwdlMapping_twdlAddr_raw[4:2];
          Radix22TwdlMapping_twdl45Reg_next = 1'b0;
        end
    endcase
    Radix22TwdlMapping_octantReg1_next = Radix22TwdlMapping_octant;
    case ( Radix22TwdlMapping_octant)
      3'b000 :
        begin
          Radix22TwdlMapping_twdlAddrMap_next = Radix22TwdlMapping_twdlAddr_raw[1:0];
        end
      3'b001 :
        begin
          Radix22TwdlMapping_t_1 = {2'b0, Radix22TwdlMapping_twdlAddr_raw};
          Radix22TwdlMapping_sub_temp_0 = 7'sb0001000 - Radix22TwdlMapping_t_1;
          Radix22TwdlMapping_twdlAddrMap_next = Radix22TwdlMapping_sub_temp_0[1:0];
        end
      3'b010 :
        begin
          Radix22TwdlMapping_t_2_0 = {2'b0, Radix22TwdlMapping_twdlAddr_raw};
          Radix22TwdlMapping_sub_temp_1 = Radix22TwdlMapping_t_2_0 - 7'sb0001000;
          Radix22TwdlMapping_twdlAddrMap_next = Radix22TwdlMapping_sub_temp_1[1:0];
        end
      3'b011 :
        begin
          Radix22TwdlMapping_sub_cast_0 = {5'b0, {Radix22TwdlMapping_twdlAddr_raw, 2'b00}};
          Radix22TwdlMapping_sub_temp_2 = 12'sb000001000000 - Radix22TwdlMapping_sub_cast_0;
          Radix22TwdlMapping_twdlAddrMap_next = Radix22TwdlMapping_sub_temp_2[3:2];
        end
      3'b100 :
        begin
          Radix22TwdlMapping_sub_cast_1 = {5'b0, {Radix22TwdlMapping_twdlAddr_raw, 2'b00}};
          Radix22TwdlMapping_sub_temp_3 = Radix22TwdlMapping_sub_cast_1 - 12'sb000001000000;
          Radix22TwdlMapping_twdlAddrMap_next = Radix22TwdlMapping_sub_temp_3[3:2];
        end
      default :
        begin
          Radix22TwdlMapping_sub_cast = {5'b0, {Radix22TwdlMapping_twdlAddr_raw, 2'b00}};
          Radix22TwdlMapping_sub_temp = 12'sb000001100000 - Radix22TwdlMapping_sub_cast;
          Radix22TwdlMapping_twdlAddrMap_next = Radix22TwdlMapping_sub_temp[3:2];
        end
    endcase
    if (Radix22TwdlMapping_phase == 2'b00) begin
      Radix22TwdlMapping_twdlAddr_raw_next = 5'b00000;
    end
    else if (Radix22TwdlMapping_phase == 2'b01) begin
      Radix22TwdlMapping_t_0_0 = {2'b0, Radix22TwdlMapping_cnt};
      Radix22TwdlMapping_twdlAddr_raw_next = Radix22TwdlMapping_t_0_0 <<< 8'd1;
    end
    else if (Radix22TwdlMapping_phase == 2'b10) begin
      Radix22TwdlMapping_twdlAddr_raw_next = {2'b0, Radix22TwdlMapping_cnt};
    end
    else begin
      Radix22TwdlMapping_cnt_cast = {2'b0, Radix22TwdlMapping_cnt};
      Radix22TwdlMapping_twdlAddr_raw_next = (Radix22TwdlMapping_cnt_cast <<< 8'd1) + Radix22TwdlMapping_cnt_cast;
    end
    Radix22TwdlMapping_phase_next = 2'b10;
    Radix22TwdlMapping_cnt_next = 3'b010;
    twdlAddr = Radix22TwdlMapping_twdlAddrMap;
    twdlAddrVld = Radix22TwdlMapping_dvldReg2;
    twdlOctant = Radix22TwdlMapping_octantReg1;
    twdl45 = Radix22TwdlMapping_twdl45Reg;
  end



  // Twiddle ROM1
  assign Twiddle_re_table_data[0] = 12'sb010000000000;
  assign Twiddle_re_table_data[1] = 12'sb001111101100;
  assign Twiddle_re_table_data[2] = 12'sb001110110010;
  assign Twiddle_re_table_data[3] = 12'sb001101010011;
  assign twiddleS_re = Twiddle_re_table_data[twdlAddr];



  always @(posedge clk or posedge reset)
    begin : TWIDDLEROM_RE_process
      if (reset == 1'b1) begin
        twiddleReg_re <= 12'sb000000000000;
      end
      else begin
        if (enb) begin
          twiddleReg_re <= twiddleS_re;
        end
      end
    end



  // Twiddle ROM2
  assign Twiddle_im_table_data[0] = 12'sb000000000000;
  assign Twiddle_im_table_data[1] = 12'sb111100111000;
  assign Twiddle_im_table_data[2] = 12'sb111001111000;
  assign Twiddle_im_table_data[3] = 12'sb110111000111;
  assign twiddleS_im = Twiddle_im_table_data[twdlAddr];



  always @(posedge clk or posedge reset)
    begin : TWIDDLEROM_IM_process
      if (reset == 1'b1) begin
        twiddleReg_im <= 12'sb000000000000;
      end
      else begin
        if (enb) begin
          twiddleReg_im <= twiddleS_im;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_process
      if (reset == 1'b1) begin
        twdlOctantReg <= 3'b000;
      end
      else begin
        if (enb) begin
          twdlOctantReg <= twdlOctant;
        end
      end
    end



  always @(posedge clk or posedge reset)
    begin : intdelay_1_process
      if (reset == 1'b1) begin
        twdl45Reg <= 1'b0;
      end
      else begin
        if (enb) begin
          twdl45Reg <= twdl45;
        end
      end
    end



  // Radix22TwdlOctCorr
  always @(twdl45Reg, twdlOctantReg, twiddleReg_im, twiddleReg_re) begin
    Radix22TwdlOctCorr_cast_0 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_2 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_4 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_6 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_8 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_10 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_3 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_9 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_1 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_7 = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast = 13'sb0000000000000;
    Radix22TwdlOctCorr_cast_5 = 13'sb0000000000000;
    Radix22TwdlOctCorr_twdlIn_re = twiddleReg_re;
    Radix22TwdlOctCorr_twdlIn_im = twiddleReg_im;
    if (twdl45Reg) begin
      case ( twdlOctantReg)
        3'b000 :
          begin
            Radix22TwdlOctCorr_twdlIn_re = 12'sb001011010100;
            Radix22TwdlOctCorr_twdlIn_im = 12'sb110100101100;
          end
        3'b010 :
          begin
            Radix22TwdlOctCorr_twdlIn_re = 12'sb110100101100;
            Radix22TwdlOctCorr_twdlIn_im = 12'sb110100101100;
          end
        3'b100 :
          begin
            Radix22TwdlOctCorr_twdlIn_re = 12'sb110100101100;
            Radix22TwdlOctCorr_twdlIn_im = 12'sb001011010100;
          end
        default :
          begin
            Radix22TwdlOctCorr_twdlIn_re = 12'sb001011010100;
            Radix22TwdlOctCorr_twdlIn_im = 12'sb110100101100;
          end
      endcase
    end
    else begin
      case ( twdlOctantReg)
        3'b000 :
          begin
          end
        3'b001 :
          begin
            Radix22TwdlOctCorr_cast = {twiddleReg_im[11], twiddleReg_im};
            Radix22TwdlOctCorr_cast_0 =  - (Radix22TwdlOctCorr_cast);
            Radix22TwdlOctCorr_twdlIn_re = Radix22TwdlOctCorr_cast_0[11:0];
            Radix22TwdlOctCorr_cast_5 = {twiddleReg_re[11], twiddleReg_re};
            Radix22TwdlOctCorr_cast_6 =  - (Radix22TwdlOctCorr_cast_5);
            Radix22TwdlOctCorr_twdlIn_im = Radix22TwdlOctCorr_cast_6[11:0];
          end
        3'b010 :
          begin
            Radix22TwdlOctCorr_twdlIn_re = twiddleReg_im;
            Radix22TwdlOctCorr_cast_7 = {twiddleReg_re[11], twiddleReg_re};
            Radix22TwdlOctCorr_cast_8 =  - (Radix22TwdlOctCorr_cast_7);
            Radix22TwdlOctCorr_twdlIn_im = Radix22TwdlOctCorr_cast_8[11:0];
          end
        3'b011 :
          begin
            Radix22TwdlOctCorr_cast_1 = {twiddleReg_re[11], twiddleReg_re};
            Radix22TwdlOctCorr_cast_2 =  - (Radix22TwdlOctCorr_cast_1);
            Radix22TwdlOctCorr_twdlIn_re = Radix22TwdlOctCorr_cast_2[11:0];
            Radix22TwdlOctCorr_twdlIn_im = twiddleReg_im;
          end
        3'b100 :
          begin
            Radix22TwdlOctCorr_cast_3 = {twiddleReg_re[11], twiddleReg_re};
            Radix22TwdlOctCorr_cast_4 =  - (Radix22TwdlOctCorr_cast_3);
            Radix22TwdlOctCorr_twdlIn_re = Radix22TwdlOctCorr_cast_4[11:0];
            Radix22TwdlOctCorr_cast_9 = {twiddleReg_im[11], twiddleReg_im};
            Radix22TwdlOctCorr_cast_10 =  - (Radix22TwdlOctCorr_cast_9);
            Radix22TwdlOctCorr_twdlIn_im = Radix22TwdlOctCorr_cast_10[11:0];
          end
        default :
          begin
            Radix22TwdlOctCorr_twdlIn_re = twiddleReg_im;
            Radix22TwdlOctCorr_twdlIn_im = twiddleReg_re;
          end
      endcase
    end
    twdl_3_19_re_1 = Radix22TwdlOctCorr_twdlIn_re;
    twdl_3_19_im_1 = Radix22TwdlOctCorr_twdlIn_im;
  end



  assign twdl_3_19_re = twdl_3_19_re_1;

  assign twdl_3_19_im = twdl_3_19_im_1;

endmodule  // TWDLROM_3_19

