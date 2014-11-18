#!/bin/sh

message() {
    bold=`tput bold`
    normal=`tput sgr0`
    echo
    echo "This script runs the octave file to obtain the values of the variables"
    echo "to be replaced in the declarations section of the UPPAAL model."
    echo
    echo "${bold}Usage:${normal}"
    echo "$0 {-g} {model} {n} {WLANStandard} {channelWidthMHz}"
    echo "        {dataRateMbps} {controlRateMbps} {appLayerPayLoadBytes}"
    echo
    echo "${bold}flags${normal}"
    echo "    -g  ::=  update global definitions."
    echo
    echo "${bold}options${normal}"
    echo "    model  ::=  {path to UPPAAL .xml file}"
    echo
    echo "    n   ::=  {number of nodes}"
    echo 
    echo "    WLANStandard  ::=  {0|1|2}"
    echo "        0 = 802.11a"
    echo "        1 = 802.11g-only (fast slot duration)"
    echo "        2 = 802.11g-mixed (slow slot duration, 802.11b legacy)"
    echo 
#    echo "    m  ::=  {number of backoff stages}"
#    echo 
#    echo "    W  ::=  {initial value for the maximum contention window}"
#    echo 
#    echo "    Tbeta ::=  {sensing time the sender node waits before"
#    echo "                   transmitting in a secondary channel. Zero"
#    echo "                   for upper-bound studies and for standard DCF.}"
#    echo 
    echo "    channelWidthMHz  ::=  {5|10|20}"
    echo 
#    echo "    numberOfChannels  ::=  {integer}"
#    echo 
    echo "    dataRateMbps  ::=  {6|9|12|18|24|36|48|54}"
    echo "        This parameter refers to the modulation scheme "
    echo "    employed for data packets. Each modulation scheme is"
    echo "    represented by the data rate it achieves in the standard"
    echo "    20 MHz wide channel. E.g. BPSK 1/2 => 6Mbps (in 20 MHz)."
    echo "    The final data rate also account the channelWidthMHz."
    echo "    E.g. 12 Mbps in 10 MHz leads to 6 Mbps."
    echo 
    echo "    controlRateMbps  ::=  {6|9|12|18|24|36|48|54}"
    echo "        See dataRateMbps"
    echo 
#    echo "    rtsCtsEnabled  ::=  {0|1|2}"
#    echo "        0 means two-way handshake (DATA-ACK)"
#    echo "        1 RTS-CTS-DATA-ACK in all channels"
#    echo "        2 RTS-CTS-DATA-ACK only in primary channels"
#    echo 
    echo "    appLayerPayLoadBytes ::= {size of payload in bytes}"
    echo 
#    echo "    transportProtocol  ::=  {0|1}"
#    echo "        0 = UDP"
#    echo "        1 = TCP"
#    echo 
#    echo "    ipProtocol  ::=  {0|1}"
#    echo "        0 = ipv4"
#    echo "        1 = ipv6"
    exit 1
}

DECLFILE=declarations
HEADER="// Values extracted from octave:"
TRAILER="// End extracted values"

parseArguments() {
    MODELFILE=$2
    N=$3
    WLANStandard=$4
    channelWidthMHz=$5
    dataRateMbps=$6
    controlRateMbps=$7
    rtsCtsEnabled='0'
    appLayerPayLoadBytes=$8
}

prepareDescriptionFile() {
    echo -e "$HEADER" > "$DECLFILE"
    echo -e "// Do not modify this or the End commentary" >> "$DECLFILE"
    octave -q <<EOF >> $DECLFILE
source("octave.m")
N = $N

[TsifsMicroSecs, TdifsMicroSecs, TslotMicroSecs, TsymbolMicroSecs] = compute80211TimingParameters($channelWidthMHz,$WLANStandard)
                 
msduDataBits = $appLayerPayLoadBytes * 8;

[macPayloadTxTimeInMicroSecs, macAndPhyHeadersTxTimeMicroSecs, ackTxTimeInMicroSecs, rtsTxTimeInMicroSecs, ctsTxTimeInMicroSecs] = computeTimeToTx80211FramesWithPhyHeaders (msduDataBits, $dataRateMbps, $controlRateMbps, TsymbolMicroSecs, $channelWidthMHz, $WLANStandard, $rtsCtsEnabled);

H = macAndPhyHeadersTxTimeMicroSecs;
TL = macPayloadTxTimeInMicroSecs;
delta = 1;                              %from octave.m
Tack = ackTxTimeInMicroSecs

Tdata = H + TL + delta
Ts = H + TL + delta + TsifsMicroSecs + Tack + delta + TdifsMicroSecs

if ($channelWidthMHz == 20)
    aPHY_RX_START_DELAY = 25;               %from IEEE Std 802.11-2012, p. 1623
elseif ($channelWidthMHz == 10)
    aPHY_RX_START_DELAY = 49;
elseif ($channelWidthMHz == 5)
    aPHY_RX_START_DELAY = 97;
endif

ACKTimeout = TsifsMicroSecs + TslotMicroSecs + aPHY_RX_START_DELAY

%% optimal (minimum) time to transmit
% T1 = Ts;
% T2 = T1 + (1sifs + Ts);
% T3 = T2 + (2sifs + Ts);
%...
Tn = sum([0:N-1])*TsifsMicroSecs + N*Ts
EOF
    awk -i inplace 'NR<=2 {print $0}; NR>2 {print "const int "$0";"}' $DECLFILE
    echo "$TRAILER" >> $DECLFILE
}

replaceDescription() {
     awk -i inplace '{
        if(match($0, h)){
            skip=1;
            next;
        }else if(match($0, t)){
            system("cat "f);
            skip=0;
            next;
        }else if(skip==1){
            next;
        }
    }1' h="$HEADER" t="$TRAILER" f="$DECLFILE" "$MODELFILE"
}

[ "$1" = "" ] && message

while getopts :hg opt;
do
    case "$opt" in
	g)
	    parseArguments ${@}
	    prepareDescriptionFile
	    replaceDescription
	    rm $DECLFILE
            exit 0
            ;;
	h) message; exit 0;;
	?) message; exit 0;;
    esac
done
message
