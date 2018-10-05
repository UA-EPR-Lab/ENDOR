#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	6.500000
XWID	12.000000
*
*	Documentational Text:
*
TITL	'mimsENDOR CYP3A4 nl 2960 G 001'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'RF'
IRUNI	''
IIUNI	''
XUNI	'MHz'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    mmlockart
DATE    07/10/18
TIME    16:44:20
CMNT    10K
SAMP    CYP3A4 no ligand
SFOR    4mm tube
STAG    P
EXPT    PLS
OXS1    TADC
AXS1    B0VL
AXS2    NONE
AXS3    
MWPW    1.99e-07
A1CT    0.296
B0VL    0.296
A1SW    0.12
MWFQ    9.758864e+09
AVGS    850
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        100.0
Power              0.0001990 mW
PowerAtten         60 dB

.DVC     endor, 1.0

EIEENDORFreq       14.9021800 MHz/3.5 kG
EIEIsotope         H1
EIERFSweepDir      Same
EIEStaticField     3480.000 G
EIEStaticRF        1.000 MHz
ENDORType          EIF
RF1Atten           1.000 dB
RF1FreqPos         10.000 MHz
RF1StartFreq       1.000 MHz
RF1SweepWidth      29.000 MHz
RF2Atten           60.000 dB
RF2FreqPos         1.000 MHz
RF2StartFreq       1.000 MHz
RF2SweepWidth      30.000 MHz
RFSrcMixing        Add
SumAtten           12.000 dB
SumAttenStart      10.000 dB
SumAttenWidth      70.000 dB

.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        2960.00 G
Delay              0.0 s
FieldFlyback       On
FieldPosition      2960.000 G
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     2930.000 G
SweepDirection     Up
SweepWidth         1200.0 G
WidthTM            200.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.758864 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        4.00 dB
ELDORAtt           30 dB
FrequencyA         9.500000 GHz
VideoBW            200 MHz
VideoGain          33 dB

.DVC     ftEpr, 1.0

AutoTimeOut        True
AveragesPerScan    1
ELDORFreqStart     9.300000 GHz
ELDORFreqWidth     0.700000 GHz
FTAcqModeSlct      Run from PulseSPEL
FieldIsStatic      False
GradIntPulse       False
GrdEnable          False
LastXAxis          Time
LastYAxis          Magnetic Field
MMWaveLOFreq       26.880000 GHz
PCycleAllowed      False
PCycleOn           False
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     Mims ENDOR ESE
PlsSPELGlbTxt      \
;\n\
;       ENDOR package\n\
;     PulseSPEL general variables, definitions & conventions\n\
\n\
begin defs\n\
\n\
; Variables     ; Comments:\n\
\n\
p0   = 24      ; 90 pulse length\n\
p1   = 32      ; 180 pulse length\n\
p2   = 400     ; inversion pulse\n\
p3   = 16      ; non selective pulse\n\
\n\
pg   = 22      ; Integrator Gate Width (RESERVED)\n\
d0   = 430     ; Detection Trigger offset time\n\
\n\
d1   = 240     ; tau\n\
d2   = 400     ;\n\
d3   = 1000    ; T\n\
d4   = 124     ;\n\
d5   = 400     ;\n\
\n\
d9   = 0       ;\n\
d30  = 4       ; 1st time scale increment\n\
d31  = 4       ; 2nd time scale increment\n\
\n\
a    = 50      ; (RESERVED) number of transient averages (TRANS)\n\
               ; if greater than 1024, a must be a multiple of 1024\n\
h    = 1       ; (CONVENTION) number of shots per point (INTG and SPT)\n\
n    = 850     ; (CONVENTION) number of sweeps to accumulate\n\
\n\
srt  = 1000 * srtu ; (RESERVED) SRT - Shot Repetition Time (srtu=1.02 us)\n\
\n\
w    = 800     ; (RESERVED) Split point of sweeps (ns)\n\
               ; i.e. pulse separation where TWT gate pulse can be split.\n\
dx   = 0       ; (RESERVED) 1st time scale starting value\n\
dy   = 0       ; (RESERVED) 2nd time scale starting value\n\
;\n\
;\n\
; PULSED ENDOR PARAMETERS\n\
;\n\
;\n\
rf1  = 6500    ; start frequency of RF-Source 1 in kHz\n\
rf2  = 1000    ; start frequency of RF-Source 2 in kHz (pump frequency in ENQOR and TRIPLE)\n\
\n\
df21 = 12000   ; RF sweep width in kHz of rf1\n\
df22 = 30000   ; RF sweep width in kHz of rf2\n\
\n\
df1  = 12 MHz  ; on resonance RF frequency (for nutationsMims.exp)\n\
df2  = 10 MHz  ; off resonance RF frequency (for nutationsMims.exp)\n\
\n\
p10  = 9000    ; RF pulse width\n\
p11  = 8000    ; RF pulse width\n\
d10  = 16000   ; mw delay during RF pulse\n\
d11  = 8000    ; mw delay during RF pulse\n\
\n\
d12  = 1000    ; delay before RF pulse\n\
d13  = 1000    ; delay after RF pulse\n\
d14  = 200     ; delay between RF pulses\n\
\n\
d15  = 1000    ; Nuclear Spin Echo tau\n\
\n\
df11 = 100     ; frequency increment in kHz of rf1 (will be recalculated)\n\
df12 = 100     ; frequency increment in kHz of rf2 (will be recalculated)\n\
	       ; or shift for TPPI	\n\
df9  = 100000  ;\n\
end defs\n\
\n\
; Note:\n\
;       Variables Marked 'RESERVED' = variable may be used for other times,\n\
;                                       however only this variable can change\n\
;                                        the indicated parameter.\n\
;                                       eg. only 'a' will change the number of\n\
;                                           transients averaged in SpecJet\n\
\n\
\n\
\n
PlsSPELLISTSlct    4-step
PlsSPELPhPrgEx     Normal
PlsSPELPrg         
PlsSPELPrgTxt      \
;\n\
;\n\
; Mims ENDOR with Echo Detection\n\
;  12/2004  MKB  removed d1 from T period, set d10=p10\n\
;   6/2008  MKB  changed exp2-4,8 for stochastic ENDOR\n\
;\n\
;\n\
;\n\
begin defs      ; dimension of data arrays [sx,sy]\n\
 dim s[512]     ;       for 3-pulse echo set-up\n\
 dim1 s[128]    ;       for 3-pulse echo field sweep\n\
 dim2 s[1024]    ;       for Mims ENDOR\n\
 dim3 s[256,16] ;       for 2D Mims ENDOR vs tau\n\
 dim4 s[512,16] ;       for 2D Mims ENDOR vs B0\n\
 dim5 s[128]    ;       for Rabi Oscillations\n\
 dim6 s[512]    ;       for Mims TRIPLE\n\
 dim7 s[512,16];       for 2D-Mims TRIPLE\n\
 dim8 s[451]    ;       for reflected echo detection\n\
 dim9 s[512]    ;       for reflected echo setup\n\
end defs\n\
;\n\
; Phase Cycle Programs\n\
;\n\
begin lists "none"\n\
 ph1  +x\n\
 ph2  +x\n\
 ph3  +x\n\
 ph4  +x\n\
 asg1 +a\n\
 bsg1 +b\n\
 asg2 +a\n\
 bsg2 +b\n\
end lists\n\
\n\
begin lists1 "4-step"\n\
 ph1 +x -x +x -x\n\
 ph2 +x +x -x -x\n\
 ph3 +x +x +x +x\n\
 ph4 +x +x +x +x\n\
 asg1 +a -a -a +a\n\
 bsg1 +b -b -b +b\n\
 asg2 -a +a +a -a\n\
 bsg2 -b +b +b -b\n\
end lists1\n\
\n\
begin lists2 "4-step both Triple"\n\
 ph1 +x -x +x -x\n\
 ph2 +x +x -x -x\n\
 ph3 +x +x +x +x\n\
 ph4 +x +x +x +x\n\
 asg1 +a -a -a +a\n\
 bsg1 +a -a -a +a\n\
 asg2 -a +a +a -a\n\
 bsg2 +a -a -a +a\n\
end lists2\n\
\n\
\n\
;\n\
; 3-pulse echo set-up\n\
;\n\
begin exp "3P ESE Setup" [QUAD TRANS]\n\
        ; QUAD detection of Transient\n\
\n\
     d10=p10\n\
     shot i=1 to a              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       dig [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
end exp\n\
\n\
;\n\
; 3-pulse echo field sweep\n\
;\n\
begin exp1 "3P ESE Field Sweep" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
  d10=p10\n\
  for k=1 to n                  ; field sweep averaging loop\n\
   totscans (n)                 ; output total number of scans\n\
\n\
   bsweep x=1 to sx             ; field sweep loop\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
   next x                       ; end of field sweep loop\n\
\n\
   scansdone (k)                ; output number of scans done\n\
  next k                        ; end of sweep averaging loop\n\
\n\
end exp1\n\
\n\
;\n\
; Mims ENDOR\n\
;\n\
begin exp2 "Mims ENDOR ESE" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
  df1=rf1                       ; start frequency of RF-Source 1\n\
  df11=sx - 1                   ; calculate RF increment\n\
  df11=df21/df11                ; RF sweep width / number of point\n\
  d10=p10\n\
\n\
  for k=1 to n                  ; RF sweep averaging loop\n\
   totscans (n)                 ; output total number of scans\n\
\n\
;   rf1=df1                      ; reset RF frequency for averaging\n\
   df1=rf1\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
;    rfsweep x=1 to sx           ; RF sweep loop\n\
    sweep x=1 to sx\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
;       p10 [RF1]                ; RF-pulse\n\
       p10 ["RF 1", df1]\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
;     rf1=rf1+df11               ; increment RF frequency (rf1) by df11\n\
     df1=df1+df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
   scansdone (k)                ; output number of scans done\n\
  next k                        ; end of averaging loop\n\
\n\
end exp2\n\
\n\
;\n\
; Mims ENDOR vs tau\n\
;\n\
begin exp3 "Mims ENDOR vs. tau" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
 df1=rf1                        ; start frequency of RF 1\n\
 df11=sx - 1                    ; calculate RF increment\n\
 df11=df21/df11                 ; RF sweep width / number of point\n\
 dy=d1                          ; define y-axis\n\
 d10=p10\n\
\n\
 for k=1 to n                  ; RF sweep averaging loop\n\
  totscans (n)                  ; output total number of scans\n\
\n\
  for y=1 to sy                 ; tau step loop\n\
\n\
   df1=rf1                      ; reset RF frequency for averaging\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
    sweep x=1 to sx           ; RF sweep loop\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       dy                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1", df1]                ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       dy                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
     df1=df1+df11               ; increment RF frequency (rf1) by df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
   dy=dy+d31                    ; increment tau step\n\
  next y                                ; end of tau step loop\n\
  dy=d1\n\
   \n\
  scansdone (k)         ; output number of scans done\n\
 next k                 ; end of averaging loop\n\
\n\
\n\
end exp3\n\
\n\
;\n\
; Mims ENDOR vs B0\n\
;\n\
begin exp4 "Mims ENDOR vs. B0" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
 df1=rf1                        ; start frequency of RF 1\n\
 df11=sx - 1                    ; calculate RF increment\n\
 df11=df21/df11                 ; RF sweep width / number of point\n\
 d10=p10\n\
\n\
 for k=1 to n                  ; RF sweep averaging loop\n\
  totscans (n)                  ; output total number of scans\n\
\n\
  bsweep y=1 to sy              ; B0 scan loop\n\
\n\
   df1=rf1                      ; reset RF frequency for averaging\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
    sweep x=1 to sx           ; RF sweep loop\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1", df1]                ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
     df1=df1+df11               ; increment RF frequency (rf1) by df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
  next y                        ; end of B0 scan loop\n\
   \n\
  scansdone (k)         ; output number of scans done\n\
 next k                 ; end of averaging loop\n\
\n\
\n\
end exp4\n\
\n\
;\n\
; Rabi Oscillations in Mims ENDOR sequence\n\
;\n\
begin exp5 "Trans-Nutation (Rabi Osc)" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
  df1=rf1                       ; frequency of RF-Source 1\n\
  d8=p10                        ; store starting RF pulse width\n\
  d10=sx-1                      ; calculate maximum mw-delay during RF pulse\n\
  d10=d10*d30                   ; maximum step added to p10\n\
  d10=p10+d10                   ; the overall maximum value of p10\n\
\n\
  for k=1 to n                  ; averaging loop\n\
    totscans (n)                ; output total number of scans\n\
\n\
    p10=d8                      ; reset RF pulse for averaging loop\n\
    dx=p10                      ; assignment of x-axis\n\
\n\
    sweep x=1 to sx             ; RF pulse Width sweep loop\n\
      shot i=1 to h             ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1",df1]            ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
      next i                    ; end of accumulation loop\n\
\n\
     p10=p10+d30                ; increment RF pulse width\n\
     dx=dx+d30                  ; determine x-axis\n\
    next x                      ; end RF pulse Width sweep loop\n\
\n\
  scansdone (k)                 ; output number of scans done\n\
 next k                 ; end of averaging loop\n\
\n\
end exp5\n\
\n\
;\n\
; Mims TRIPLE\n\
;\n\
begin exp6 "Mims TRIPLE" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
  df1=rf1                       ; start frequency of RF-Source 1\n\
  df11=sx - 1                   ; calculate RF increment\n\
  df11=df21/df11                ; RF sweep width / number of point\n\
  d10=p10\n\
  d11=p11                       ; pump pulse\n\
  df2=rf2                       ; second rf frequency\n\
\n\
  for k=1 to n                  ; RF sweep averaging loop\n\
   totscans (n)                 ; output total number of scans\n\
\n\
   df1=rf1                      ; reset RF frequency for averaging\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
    sweep x=1 to sx           ; RF sweep loop\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1", df1]        ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d12                      ; delay between\n\
       p11 ["RF 1", df2]        ; Pump RF-pulse\n\
       d11                      ; mw-Delay during Pump RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
     df1=df1+df11               ; increment RF frequency (rf1) by df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
   df1=rf1                      ; reset RF frequency for averaging\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
    sweep x=1 to sx           ; RF sweep loop for difference\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1", df1]                ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d12                      ; delay between\n\
       p11 [*]                  ; Pump RF-pulse\n\
       d11                      ; mw-Delay during Pump RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg2]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
     df1=df1+df11               ; increment RF frequency (rf1) by df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
\n\
   scansdone (k)                ; output number of scans done\n\
  next k                        ; end of averaging loop\n\
\n\
end exp6\n\
\n\
;\n\
; Mims 2D-TRIPLE equal length pulses\n\
;\n\
begin exp7 "Mims 2D-TRIPLE" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
  df1=rf1                       ; start frequency of RF-Source 1\n\
  df11=sx - 1                   ; calculate RF increment\n\
  df11=df21/df11                ; RF sweep width / number of point\n\
  d10=p10\n\
\n\
  df2=rf2                       ; start frequency of RF-Source 1\n\
  df12=sy - 1                   ; calculate RF increment\n\
  df12=df22/df12                ; RF sweep width / number of point\n\
  d11=d10\n\
  p11=p10\n\
 \n\
  for k=1 to n                  ; RF sweep averaging loop\n\
   totscans (n)                 ; output total number of scans\n\
\n\
   df2=rf2                      ; reset RF frequency for averaging\n\
   dy=rf2                       ; assignment of x-axis\n\
   for y=1 to sy\n\
\n\
   df1=rf1                      ; reset RF frequency for averaging\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
    sweep x=1 to sx           ; RF sweep loop\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1", df1]        ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d12                      ; delay between\n\
       p11 ["RF 1", df2]        ; Pump RF-pulse\n\
       d11                      ; mw-Delay during Pump RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
     df1=df1+df11               ; increment RF frequency (rf1) by df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
     df2=df2+df12               ; increment RF frequency (rf1) by df11\n\
     dy=dy+df12                 ; determine x-axis\n\
    next y                      ; end RF sweep loop\n\
\n\
   scansdone (k)                ; output number of scans done\n\
  next k                        ; end of averaging loop\n\
\n\
end exp7\n\
\n\
\n\
;\n\
; Mims ENDOR with refocused echo\n\
;\n\
begin exp8 "Mims ENDOR refocused ESE" [QUAD INTG]\n\
        ; QUAD detection with Integration\n\
\n\
  df1=rf1                       ; start frequency of RF-Source 1\n\
  df11=sx - 1                   ; calculate RF increment\n\
  df11=df21/df11                ; RF sweep width / number of point\n\
  d10=p10\n\
\n\
  for k=1 to n                  ; RF sweep averaging loop\n\
   totscans (n)                 ; output total number of scans\n\
\n\
   df1=rf1                      ; reset RF frequency for averaging\n\
   dx=rf1                       ; assignment of x-axis\n\
\n\
    sweep x=1 to sx           ; RF sweep loop\n\
     shot i=1 to h              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       p10 ["RF 1", df1]                ; RF-pulse\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d4                       ; delay for refocus\n\
       p1 [ph3]                 ; refocusing pi pulse\n\
       d4                       ; delay\n\
       d0                       ; constant acquisition delay\n\
       acq [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
     df1=df1+df11               ; increment RF frequency (rf1) by df11\n\
     dx=dx+df11                 ; determine x-axis\n\
    next x                      ; end RF sweep loop\n\
\n\
   scansdone (k)                ; output number of scans done\n\
  next k                        ; end of averaging loop\n\
\n\
end exp8\n\
\n\
;\n\
; 3-pulse echo set-up\n\
;\n\
begin exp9 "refocused echo Setup" [QUAD TRANS]\n\
        ; QUAD detection of Transient\n\
\n\
     d10=p10\n\
     shot i=1 to a              ; accumulation loop\n\
       d9                       ; DAF (Delay After Flash)\n\
       p0 [ph1]                 ; mw- 1st pulse (90 degrees)\n\
       d1                       ; tau\n\
       p0 [ph2]                 ; mw- 2nd pulse (90 degrees)\n\
       d12                      ; Delay after Preparation (d12 > p0)\n\
       d10                      ; mw-Delay during RF pulse \n\
       d13                      ; Delay after RF pulse\n\
       p0 [ph3]                 ; mw- 3rd pulse (90 degrees)\n\
       d1                       ; tau\n\
       d4                       ; delay for refocus\n\
       p1 [ph3]                 ; refocusing pi pulse\n\
       d4                       ; delay\n\
       d0                       ; constant acquisition delay\n\
       dig [sg1]                ; acquisition\n\
     next i                     ; end of accumulation loop\n\
\n\
end exp9\n\
\n
Psd1               {2;33,4;0} 0,0,256,496,11496,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,44,44,44,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd10              {2;33,4;0} 0,260,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd11              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd12              {2;33,4;0} 0,260,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd13              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd14              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd15              {2;33,4;0} 0,0,1240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd16              {2;33,4;0} 0,0,1240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd17              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd18              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd19              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd2               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd20              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd21              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd22              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd23              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd24              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd25              {2;33,4;0} 0,262,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd26              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd27              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd28              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd29              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd3               {2;33,4;0} 0,0,150,11390,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,482,242,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd4               {2;33,4;0} 0,0,0,11240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,452,212,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd5               {2;33,4;0} 0,0,11910,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd6               {2;33,4;0} 0,260,260,500,11500,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,24,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd7               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd8               {2;33,4;0} 0,260,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd9               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadDetect         True
RF1Prg             {2;33,3;0.000[MHz]} 0.000,2.000,6.500,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.012,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
RF2Prg             {2;33,3;0.000[MHz]} 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
ReplaceMode        Off
ShotRepTime        1020.00 us
ShotsPLoop         1
SmoothAllowed      False
SmoothPoints       1
SptProgress        100 %
StochMode          True
SweepsPExp         850
TriggerTimeOut     11 s
XAxisQuant         Time
XSpecRes           1024
YAxisQuant         Magnetic Field
YSpecRes           1

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         850
NbScansDone        850
NbScansToDo        850
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
