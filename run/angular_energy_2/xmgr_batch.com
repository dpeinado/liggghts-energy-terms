FOCUS G0
READ NXY "INFILE"
view 0.140000, 0.070000, 1.250000, 0.980000
    world 0, 0, 1500, 600000
xaxis  on
    xaxis  type zero false
    xaxis  offset 0.000000 , 0.000000
    xaxis  bar on
    xaxis  bar color 1
    xaxis  bar linestyle 1
    xaxis  bar linewidth 1.0
    xaxis  label "Time Step"
    xaxis  label layout para
    xaxis  label place auto
    xaxis  label char size 1.000000
    xaxis  label font 0
    xaxis  label color 1
    xaxis  label place normal
    xaxis  tick on
    xaxis  tick major 500
    xaxis  tick minor ticks 1
    xaxis  tick default 4
    xaxis  tick place rounded true
    xaxis  tick in
    xaxis  tick major size 1.000000
    xaxis  tick major color 1
    xaxis  tick major linewidth 1.0
    xaxis  tick major linestyle 1
    xaxis  tick major grid on
    xaxis  tick minor color 1
    xaxis  tick minor linewidth 1.0
    xaxis  tick minor linestyle 2
    xaxis  tick minor grid on
    xaxis  tick minor size 0.500000
    xaxis  ticklabel on
    xaxis  ticklabel format general
    xaxis  ticklabel prec 5
    xaxis  ticklabel formula ""
    xaxis  ticklabel append ""
    xaxis  ticklabel prepend ""
    xaxis  ticklabel angle 0
    xaxis  ticklabel skip 0
    xaxis  ticklabel stagger 0
    xaxis  ticklabel place normal
    xaxis  ticklabel offset auto
    xaxis  ticklabel offset 0.000000 , 0.010000
    xaxis  ticklabel start type auto
    xaxis  ticklabel start 0.000000
    xaxis  ticklabel stop type auto
    xaxis  ticklabel stop 0.000000
    xaxis  ticklabel char size 1.000000
    xaxis  ticklabel font 0
    xaxis  ticklabel color 1
    xaxis  tick place both
    xaxis  tick spec type none
    yaxis  on
    yaxis  type zero false
    yaxis  offset 0.000000 , 0.000000
    yaxis  bar on
    yaxis  bar color 1
    yaxis  bar linestyle 1
    yaxis  bar linewidth 1.0
    yaxis  label "Energy"
    yaxis  label layout para
    yaxis  label place auto
    yaxis  label char size 1.000000
    yaxis  label font 0
    yaxis  label color 1
    yaxis  label place normal
    yaxis  tick on
    yaxis  tick major 100000
    yaxis  tick minor ticks 1
    yaxis  tick default 6
    yaxis  tick place rounded true
    yaxis  tick in
    yaxis  tick major size 1.000000
    yaxis  tick major color 1
    yaxis  tick major linewidth 1.0
    yaxis  tick major linestyle 1
    yaxis  tick major grid on
    yaxis  tick minor color 1
    yaxis  tick minor linewidth 1.0
    yaxis  tick minor linestyle 2
    yaxis  tick minor grid on
    yaxis  tick minor size 0.500000
    yaxis  ticklabel on
    yaxis  ticklabel format general
    yaxis  ticklabel prec 5
    yaxis  ticklabel formula ""
    yaxis  ticklabel append ""
    yaxis  ticklabel prepend ""
    yaxis  ticklabel angle 0
    yaxis  ticklabel skip 0
    yaxis  ticklabel stagger 0
    yaxis  ticklabel place normal
    yaxis  ticklabel offset auto
    yaxis  ticklabel offset 0.000000 , 0.010000
    yaxis  ticklabel start type auto
    yaxis  ticklabel start 0.000000
    yaxis  ticklabel stop type auto
    yaxis  ticklabel stop 0.000000
    yaxis  ticklabel char size 1.000000
    yaxis  ticklabel font 0
    yaxis  ticklabel color 1
    yaxis  tick place both
    yaxis  tick spec type none
    legend on
    legend loctype view
    legend 1.094, 0.98
    legend box color 1
    legend box pattern 1
    legend box linewidth 1.0
    legend box linestyle 1
    legend box fill color 0
    legend box fill pattern 1
    legend font 0
    legend char size 1.000000
    legend color 1
    legend length 4
    legend vgap 1
    legend hgap 1
    legend invert false
AUTOSCALE
LEGEND ON

SWAP G0.S3 AND G0.S4
SWAP G0.S5 AND G0.S6
SWAP G0.S6 AND G0.S7
SWAP G0.S6 AND G0.S8
SWAP G0.S3 AND G0.S6


S0 LINE COLOR 1
S1 LINE COLOR 1
S2 LINE COLOR 1
S3 LINE COLOR 1
S4 LINE COLOR 1
S5 LINE COLOR 1
S6 LINE COLOR 1
S7 LINE COLOR 1
S8 LINE COLOR 1
S9 LINE COLOR 1

S0 LEGEND "Ekin"
S1 LEGEND "Ekrot"
S2 LEGEND "Un"
S3 LEGEND "Econs"
S4 LEGEND "Dnv"
S5 LEGEND "Wt"
S6 LEGEND "Dtv"
S7 LEGEND "Dtf"
S8 LEGEND "HDE"
S9 LEGEND "Etot"

S0 LINESTYLE 1
S1 LINESTYLE 1
S2 LINESTYLE 1
S3 LINESTYLE 1
S4 LINESTYLE 4
S5 LINESTYLE 1
S6 LINESTYLE 3
S7 LINESTYLE 3
S8 LINESTYLE 2
S9 LINESTYLE 1

S0 LINEWIDTH 2
S1 LINEWIDTH 2
S2 LINEWIDTH 2
S3 LINEWIDTH 2
S4 LINEWIDTH 2
S5 LINEWIDTH 2
S6 LINEWIDTH 2
S7 LINEWIDTH 2
S8 LINEWIDTH 4
S9 LINEWIDTH 2


S0 SYMBOL 1
S1 SYMBOL 2
S2 SYMBOL 3
S3 SYMBOL 4
S4 SYMBOL 0
S5 SYMBOL 0
S6 SYMBOL 8
S7 SYMBOL 9
S8 SYMBOL 0
S9 SYMBOL 10

S0 SYMBOL LINEWIDTH 2
S1 SYMBOL LINEWIDTH 2
S2 SYMBOL LINEWIDTH 2
S3 SYMBOL LINEWIDTH 2
S4 SYMBOL LINEWIDTH 2
S5 SYMBOL LINEWIDTH 2
S6 SYMBOL LINEWIDTH 2
S7 SYMBOL LINEWIDTH 2
S8 SYMBOL LINEWIDTH 2
S9 SYMBOL LINEWIDTH 2

S0 SYMBOL color 1
S1 SYMBOL color 1
S2 SYMBOL color 1
S3 SYMBOL color 1
S4 SYMBOL color 1
S5 SYMBOL color 1
S6 SYMBOL color 1
S7 SYMBOL color 1
S8 SYMBOL color 1
S9 SYMBOL color 1

S0 SYMBOL size 0.8
S1 SYMBOL size 0.8
S2 SYMBOL size 0.8
S3 SYMBOL size 0.8
S4 SYMBOL size 0.8
S5 SYMBOL size 0.8
S6 SYMBOL size 0.8
S7 SYMBOL size 0.8
S8 SYMBOL size 0.8
S9 SYMBOL size 0.8


S0 SYMBOL skip 10
S1 SYMBOL skip 10
S2 SYMBOL skip 10
S3 SYMBOL skip 10
S4 SYMBOL skip 10
S5 SYMBOL skip 10
S6 SYMBOL skip 10
S7 SYMBOL skip 10
S8 SYMBOL skip 10
S9 SYMBOL skip 10
world xmax 1500
REDRAW
PRINT TO "OUTFILE"
HARDCOPY DEVICE "PDF"
PAGE SIZE 2660, 2048
PRINT
