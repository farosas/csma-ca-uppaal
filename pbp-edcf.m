%    Performance model for the PbP-EDCF MAC Protocol    
%    Copyright (c) 2014 Saulo Queiroz (ssaulojorge@gmail.com).
%    All rights reserved.
%    Please, if you use this code as base for your research
%    include the following citation:
%    Saulo Queiroz and Roberto Hexsel, 2014. Translating Full Duplex into Capacity Gains for the High-Priority Traffic Classes of IEEE 
%    802.11. In Proceedings of the 30th Annual ACM SIGAPP Symposium on Applied Computing (SAC '15). ACM, New York, NY, USA. 

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%
%========================================>Compute the saturation throughput following the Markovian Bianchi Model
%---------> n = number of nodes
%
%--------->WLANStandard = {0,1,2}
%---> 0 = 802.11a
%---> 1 = 802.11g-only (fast slot duration)
%---> 2 = 802.11g-mixed (slow slot duration, 802.11b legacy)
%
%---------> m = number of backoff stages
%
%---------> W = initial value for the maximum contention window
%
%---------> Tbeta = sensing time the sender node waits before
%                   transmitting in a secondary channel. Zero
%                   for upper-bound studies and for standard DCF.
%
%--------->channelWidthMHz = {5,10,20} MHz
%
%--------->numberOfChannels integer
%
%--------->dataRateMbps = {6,9,12,18,24,36,48,54}
%    Actually this parameter refers to the modulation scheme 
%    employed for data packets. Each modulation scheme is
%    represented by the data rate it achieves in the standard
%    20 MHz wide channel. E.g. BPSK 1/2 => 6Mbps (in 20 MHz).
%    The final data rate also account the channelWidthMHz.
%    E.g. 12 Mbps in 10 MHz leads to 6 Mbps.
%
%--------->controlRateMbps = {6,9,12,18,24,36,48,54}. See dataRateMbps
%
%--------->rtsCtsEnabled = {0, 1, 2}
%---> 0 means two-way handshake (DATA-ACK)
%---> 1 RTS-CTS-DATA-ACK in all channels
%---> 2 RTS-CTS-DATA-ACK only in primary channels
%
%--------->appLayerPayLoadBytes: size of payload in bytes
%
%--------->transportProtocol = {0,1}
%---> 0 = UDP; 
%---> 1 = TCP
%
%--------->ipProtocol = {0,1}
%---> 0 = ipv4; 
%---> 1 = ipv6;
%
%====>References
%---> http://www.oreillynet.com/pub/a/wireless/2003/08/08/wireless_throughput.html?page=1
%---> A Case for Adapting Channel Width in Wireless Networks
%---> IEEE 802.11-2012 Standard
function [saturationThroughputMbps, tau, p]= \
         BianchiPbPDCFMarkovModel(n, WLANStandard, m, W, Tbeta, channelWidthMHz, numberOfChannels, \
         					dataRateMbps, controlRateMbps, rtsCtsEnabled, \
         					appLayerPayLoadBytes, transportProtocol, ipProtocol)

  %====================================================>  Preparing MSDU, the MAC payload. It is: app layer + UDP/TCP + IP + LLC
  msduDataBits =  appLayerPayLoadBytes * 8;
  msduDataBits += retTransportProtocolHeaderBits(transportProtocol);
  msduDataBits += retIPHeaderBits(ipProtocol);
  msduDataBits += retLLCHeaderBits();

  %====================================================>  Computing 802.11 Timing Parameters 
  % To calculate the time taken by data and headers we need first to
  % calculate the time duration of a single symbol. It depends on the
  % channel width and the WLAN standard.
  [TsifsMicroSecs, TdifsMicroSecs, TslotMicroSecs, TsymbolMicroSecs] = compute80211TimingParameters(channelWidthMHz, WLANStandard);

  %====================================================>  Computing the Time required to send:
  %  - macpayload (msduDataBits), 
  %  - MAC header/fcs and  PHY header of the MSDU frame
  %  - Time taken by ACK, RTS and CTS. Note RTS/CTS returns zero
  %    if they are not enabled.
  [macPayloadTxTimeInMicroSecs, macAndPhyHeadersTxTimeMicroSecs, ackTxTimeInMicroSecs, rtsTxTimeInMicroSecs, ctsTxTimeInMicroSecs] = \
                    computeTimeToTx80211FramesWithPhyHeaders (msduDataBits, dataRateMbps, controlRateMbps, TsymbolMicroSecs, \
                                                              channelWidthMHz, WLANStandard, rtsCtsEnabled);

  %====================================================>  Compute the basic node probabilities: p, tau0, tayC
  % These probabilities refer to the behavio of a node, as capture by the Markov model.
  % p = collision probability (for PbP-DCF refers to the  primary channel only)
  % tau0 = transmission probability (for PbP-DCF refers to the  primary channel only)
  % tauC = transmission probability in the C-th secondary channel (zero for standard DCF)
  % overallTxProbability = system transmission probability
  [p, tau0, tauC, overallTxProbability] = computeMarkovProbabilities (n, m, W, numberOfChannels);  
  tau = overallTxProbability;

  %====================================================>  Compute the channel probabilities
  %---------> Basic (preliminary) channel probabilities:
  [Ptr0, Ps0] = basicPreliminaryChannelProbabilities (tau0, n);
  [PtrC, PsC] = basicPreliminaryChannelProbabilities (tauC, n); % zero for standard DCF
  %---------> Channel probabilities
  [Psuc0, Pcol0, Pidl0] = basicChannelProbabilities (Ptr0, Ps0);
  [PsucC, PcolC, PidlC] = basicChannelProbabilities (PtrC, PsC);
 

  %====================================================>  Computing Saturation Throughput S
  %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-> The denominator of S
  %------------------------>  Preliminaries variables
  %  H      = duration of phy and mac overheads in microsecs
  %  TL     = time in micro secs. to send the MSDU (payload offered by the MAC layer) 
  %  delta  = prop. delay in microseconds
  %  Tack   = duration of ack including its phy overheads
  %  Trts   = duration of RTS including its phy overheads
  %  Tcts   = duration of CTS including its phy overheads
  %  Ts     = Time taken by a successful primary channel transmission
  %  Tcol     = Time spent by a collision in the primary channel
  %  Tslinha      = Time taken by a successful secondary channel transmission
  %  NsigmaTL     = Number of empty time slots that fits into the TL time interval
  %  NsigmaSlinha = Number of empty time slots that fits into the Tslinha time interval
  H = macAndPhyHeadersTxTimeMicroSecs;
  TL = macPayloadTxTimeInMicroSecs;
  delta = 1; %prop. delay in 
  Tack = ackTxTimeInMicroSecs;
  Trts = rtsTxTimeInMicroSecs;
  Tcts = ctsTxTimeInMicroSecs;

  if (rtsCtsEnabled == 1 || rtsCtsEnabled == 2)
      %-----> In both settings primary channel 
      %       transmissions employ RTS/CTS     
      Ts = Trts + delta + TsifsMicroSecs + Tcts + delta + TsifsMicroSecs + \
           H + TL + delta + TsifsMicroSecs + Tack + delta + TdifsMicroSecs;
      Tcol = Trts + delta + TdifsMicroSecs;
      %-----> Secondary channel transmissions 
      %       have not DIFS but have Tbeta
      Tslinha = Ts - (TdifsMicroSecs + Tbeta);
      if (rtsCtsEnabled == 2)
          %-----> In this setting  ONLY the primary
          %       transmissions employ RTS/CTS
          Tslinha -= (Trts + delta + TsifsMicroSecs + Tcts + delta + TsifsMicroSecs);
      endif
  else
      %-----> In this setting the two-way handshake
      %      is employed for all channels
      Ts = H + TL + delta + TsifsMicroSecs + Tack + delta + TdifsMicroSecs;
      Tcol = H + TL + delta + TdifsMicroSecs;
      Tslinha = Ts - (TdifsMicroSecs + Tbeta);
  endif
  NsigmaTL = ceil (TL/TslotMicroSecs);
  NsigmaSlinha = ceil (Tslinha/TslotMicroSecs);

  saturationThroughputSDenominator = Psuc0 * Ts + Pcol0 * Tcol +  Pidl0 * TslotMicroSecs;

  %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-> The numerator of S
  %------------------------>  Preliminaries
  %---------> Expected payload for primary channels (S numerator part A)
  % expectedpayLoad = E[L] = payload in primary channel
  L =  msduDataBits;

  %---------> Expected payload for secondary channels (S numerator part B)
  %   Lsuc and Lcol are the expected secondary payLoad when primary channel
  %  either contains a successful transmission or a collision, respectively.
  %   OBS.: This part becomes zero in the throughput for the standard DCF
  Lcol = Lsuc = L / 2;

  %---------> Expected payload Lsigma for secondary channels when primay is idle
  Lsigma = msduDataBits*(TslotMicroSecs/TL)*(NsigmaTL/NsigmaSlinha);


  %---------> Expected (average) payload in the primary channel.
  %      We assume a system where all packets have the same size.
  primaryChannelPayload = Psuc0 * L;

  %---------> Expected (average) secondary channel payload when 
  %  there is a collision in the primary channel and a successful 
  %  transmission in the secondary channel c
  secondaryPayloadWhenPrimaryCollides = (Pcol0 * PsucC) * Lcol;

  %---------> Expected (average) secondary channel payload when 
  %  there is a successful transmission in the primary and 
  %  secondary channels, respectively secondary channel c
  secondaryPayloadWhenPrimarySucceeds = (Psuc0 * PsucC) * Lcol;

  %---------> Expected (average) secondary channel payload when 
  %  time slot is empty in the primary channel
  secondaryPayloadWhenPrimayIsEmpty = (Pidl0 * PsucC) * Lcol;

  %---------> Total expected (average) added by all secondary channels
  % Note that Nc = 1 for the standard IEEE 802.11 DCF, then Ldelta = 0.
  Ldelta = (numberOfChannels - 1)*(secondaryPayloadWhenPrimaryCollides + secondaryPayloadWhenPrimarySucceeds + secondaryPayloadWhenPrimayIsEmpty);

  %---------> Numerator of the DCF / PbP-DCF Saturation Throughput
  saturationThroughputSNumerator = primaryChannelPayload + Ldelta;


  %=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-> The DCF / PbP-DCF Saturation Throughput S
  %---------> S in bits per micro seconds
   S = saturationThroughputSNumerator / saturationThroughputSDenominator;
  
  %---------> S in megabits per second
  %           Online calculator: http://www.endmemo.com/convert/data%20transfer.php
  saturationThroughputMbps = S * 0.953674;
endfunction



function  [Ptr, Ps] = basicPreliminaryChannelProbabilities (tau, n)
  Ptr = 1 - (1 - tau)^n;
  %----> tau = 0 means the probability to access the PbP-DCF secondary channel
  %      in the standard DCF, for instance
  if (Ptr != 0)
    Ps  = (n * tau * (1 - tau)^(n-1)) / Ptr;
  else
    Ps = 0;
  endif
endfunction

function  [Psuc, Pcol, Pidl] = basicChannelProbabilities (Ptr, Ps)
  Psuc = Ps * Ptr;
  Pcol = (1 - Ps) * Ptr;
  Pidl = (1 - Ptr);
endfunction

%
function [p, tau0, tauC, overallTxProbability] = computeMarkovProbabilities (n, m, W, numberOfChannels)
  %---------> Solving the non-linear system of equations between tau0 and p
  %           value of p and tau0 will be assigned to x(1) and x(2), respectively.
  [x,info]= fsolve(@(x) tau0pSystemOfEquantions(x, n, m, W, numberOfChannels), [0.00001;0.99999]);
  p = x(1);
  tau0 = x(2);
  if (numberOfChannels == 1)
    tauC = 0;
    overallTxProbability = tau0;
  else
    tauC = (1-p)*tau0;
    overallTxProbability = tau0 + (numberOfChannels-1)*tauC;
  endif
endfunction

function y=tau0pSystemOfEquantions(x, n, m, W, Nc)
  %p==x(1) and tau0 == x(2) and 
  %y(1) = b000/(1-x(1)) - x(2);
  % Note that when Nc=1 (PbP-DCF becomes the standard DCF) tau0 becomes
  % Eq. 6 of 'Performance Analysis of the IEEE 802.11 Distributed Coordination Function'
  % G. Bianchi, JSAC 2000.
  y(1) = ((2*(1-2*x(1))) / ((1-2*x(1))*(W+1) + x(1)*W*(1-(2*x(1))^m) + (Nc-1)*(1-x(1))*(2*W)*( 1 - x(1)*(1+ (2*x(1))^m)))) - x(2);

  % Bianchi's transmission probability tau (tau0 for PbP-DCF)
  % as function of the collision probability:
  % p = 1 - (1 - tau0)^(n-1) => 0 = - p + 1 - (1 - tau0)^(n-1)
  y(2)=-x(1)+1-(1-x(2))^(n-1); 
endfunction 


% Given bytes from the layer above MAC we compute:
% -> time to send the mac payload (msdu)
% -> time to send the mac overheads (header and trailer-fcs)
% -> time to send the phy overheads (plcp header and preambles
% -> time to send 802.11 overhead frames with their respective
%    PHY overheads (ACK + PHY headers and  RTS/CTS + PHY headers)
function [macPayloadTxTimeInMicroSecs, macAndPhyHeadersTxTimeMicroSecs, ackTxTimeInMicroSecs,\
         rtsTxTimeInMicroSecs, ctsTxTimeInMicroSecs] = computeTimeToTx80211FramesWithPhyHeaders (msduDataBits, dataRateMbps,\
                                                                                                 controlRateMbps, TsymbolMicroSecs,\
                                                                                                 channelWidthMHz, WLANStandard,\
                                                                                                 rtsCtsEnabled)
  %==========================> Preparing the MAC overheads (header and trailer)
  % We keep the MAC headers/overheads aside from the msdu (mac payload)
  % The reason is that in the saturation throughput equation we need
  % to clearly separate the time taken with MAC header from the time
  % taken with the MAC payload.
  % Note that: mpduDataBits = msduDataBits + macHeaderOfDataBits;
  [macHeaderBits, macFCSTrailerBits] = ret80211MACOverheadsBits();
  macHeaderOfDataBits = macHeaderBits + macFCSTrailerBits;


  %==========================> Remark: Adding the PHY bits that account as payload
  % The PHY layer has two type of overheads: 1) header + preamble, which
  % are sent at a constant rate (BPSK 1/2 which yields 6 Mbps in 20 MHz)
  % 2) modulation dependent bits (service, tail and pad bits) that are
  % added to the PHY payload (to make it multiple of the total data bits 
  % carried in a symbol generated at a given data rate). These bits are 
  % sent at the same rate of the phy payload (PSDU). In practice, the PSDU
  % is the MPDU i.e. MAC Payload (MSDU) + MAC overheads (header and fcs-trailer).
  % However we the PSDU = MSDU + PHY bits of type 2 (service and tail).
  % We calculate PHY preamble and header later, since they are constants.
  [servicePhyBits, tailPhyBits] = retRateDependentPhyOverHeadBits();
  msduDataBits += (servicePhyBits + tailPhyBits);


  %==========================> Computing 802.11 Timing Parameters 
  % To calculate the time taken by data and headers we need first to
  % calculate the time duration of a single symbol. It depends on the
  % channel width and the WLAN standard.
  %[TsifsMicroSecs, TdifsMicroSecs, TslotMicroSecs, TsymbolMicroSecs] = Compute80211TimingParameters(channelWidthMHz, WLANStandard);

  %==========================> Computing the number of symbols
  %---------> Number of symbols of the MAC payload
  msduDataInNumberOfOFDMSymbols = computeNumberOfOFDMSymbols(msduDataBits, dataRateMbps);
  %---------> Number of symbols of MAC overheads
  macHeaderOfDataBitsInNumberOfOFDMSymbols = computeNumberOfOFDMSymbols(macHeaderOfDataBits, dataRateMbps);


  %==========================> Computing the time to send the symbols
  %---------> Time to send only the MAC payload (denoted as TL in the paper).
  macPayloadTxTimeInMicroSecs = TsymbolMicroSecs*msduDataInNumberOfOFDMSymbols;

  %---------> Time to send MAC and PHY overheads (H in the paper)
  %--> Time to send only the MAC overheads (Hmac in Bianchi paper)
  macHeaderTxTimeInMicroSecs = TsymbolMicroSecs*macHeaderOfDataBitsInNumberOfOFDMSymbols;
  %--> Time to send the PHY overheads (Hphy in the paper)
  [plcpPreambleMicroSecs, plcpHeaderMicroSecs, signalExtensionMicroSecs] = retFixedPhyOverheads(channelWidthMHz, WLANStandard);
  phyHeaderTxTimeInMicroSecs = plcpPreambleMicroSecs + plcpHeaderMicroSecs + signalExtensionMicroSecs;
  %--> Time to send all MAC and PHY overheads (H in the paper)
  macAndPhyHeadersTxTimeMicroSecs = macHeaderTxTimeInMicroSecs + phyHeaderTxTimeInMicroSecs;


  %+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=> Preparing ACK ...
  %--> The 802.11 ACK is only 14 bytes and doesn't have TCP/UDP, IP and LLC layers.
  ackBits = 112;
  %--> Passing ACK to PHY also adds service and tail bits 
  ackBits += (servicePhyBits + tailPhyBits);

  %==========================> Computing the number of symbols  
  %---------> Number of symbols of the ACK. Note it is sent at rate controlRateMbps
  ackInNumberOfOFDMSymbols = computeNumberOfOFDMSymbols(ackBits, controlRateMbps);
  %---------> Number of symbols of MAC overheads
  macHeaderOfDataBitsInNumberOfOFDMSymbols = computeNumberOfOFDMSymbols(macHeaderOfDataBits, dataRateMbps);

  %==========================> Computing the time to send the symbols
  %---------> Time to send only the ACK without and PHY overheads
  ackTxTimeInMicroSecs = TsymbolMicroSecs * ackInNumberOfOFDMSymbols;
 
  %--> Adding the PHY overheads (PLCP and Preamble). Note it is the
  %    same as for data frames since it always use same modulation
  %    irrespective of the type of 802.11 frame (control,data,mngt.)
  ackTxTimeInMicroSecs += phyHeaderTxTimeInMicroSecs;

  %====================================================>  Computing the Time of one RTS and one CTS if needed.
  if (rtsCtsEnabled != 1 || rtsCtsEnabled != 2)
    %--> The 802.11 RTS and CTS are 20 and 14 respectively. They don't have TCP/UDP, IP and LLC layers.
    rtsBits = 160;
    ctsBits = 112;
    %--> Passing RTS and CTS to PHY also adds service and tail bits 
    rtsBits += (servicePhyBits + tailPhyBits);
    ctsBits += (servicePhyBits + tailPhyBits);

    %==========================> Computing the number of symbols  
    %---------> Number of symbols of RTS and CTS. Note they are sent at rate controlRateMbps
    rtsInNumberOfOFDMSymbols = computeNumberOfOFDMSymbols(rtsBits, controlRateMbps);
    ctsInNumberOfOFDMSymbols = computeNumberOfOFDMSymbols(ctsBits, controlRateMbps);

    %==========================> Computing the time to send the symbols
    %---------> Time to send only the RTS and CTS without and PHY overheads
    rtsTxTimeInMicroSecs = TsymbolMicroSecs * rtsInNumberOfOFDMSymbols;
    ctsTxTimeInMicroSecs = TsymbolMicroSecs * ctsInNumberOfOFDMSymbols;

    %--> Adding the PHY overheads (PLCP and Preamble). Note it is the
    %    same as for data frames since it always use same modulation
    %    irrespective of the type of 802.11 frame (control,data,mngt.)
    rtsTxTimeInMicroSecs += phyHeaderTxTimeInMicroSecs;
    ctsTxTimeInMicroSecs += phyHeaderTxTimeInMicroSecs;
  else
    rtsTxTimeInMicroSecs = ctsTxTimeInMicroSecs = 0;
 endif
endfunction


% The time taken to send the PLCP header and preamble does not depend
% on the set data rate. The standard mandates that the header (without
% the service bits) must be sent at BPSK 1/2 (which yields 6 Mbps in 
% 20 MHz). Thus the PLCP header takes only one symbol (4 micro secs
% in 20 MHz channel).  The PLCP preamble consists in the first set of 
% signals sent and are used only for synchronization purposes. It 
% takes 16 micro seconds. See sections 18.3.2 'PLCP Frame Format' and 
% 18.3.3 'PLCP preamble (SYNC)' in the IEEE 802.11-2012 standard.
function [plcpPreambleMicroSecs, plcpHeaderMicroSecs, signalExtensionMicroSecs] = retFixedPhyOverheads(channelWidth, WLANStandard)
  scalingRatio = 20 / channelWidth;

  plcpPreambleMicroSecs = 16 * scalingRatio;
  plcpHeaderMicroSecs = 4 * scalingRatio;
  switch (WLANStandard)
    case 0 
       signalExtensionMicroSecs = 0;
    case 1
     signalExtensionMicroSecs = 6;
  endswitch
endfunction

function [TsifsMicroSecs, TdifsMicroSecs, TslotMicroSecs, TsymbolMicroSecs] = compute80211TimingParameters(channelWidth, WLANStandard)
  % similar to Chandra et al 2008, we use a scalling ratio
  % to recalculate all 20 MHz-based timing parameters to
  % the value corresponding to channelWidth
  scalingRatio = 20 / channelWidth;
  
  %OBS.: To improve readibility we will set each parameter
  %      aside from each other even if they share some
  %      test conditions.

  %======> Calculating the standard 20MHz SIFS 
  %        in microseconds (us).
  %        It depends on the WLAN standard.
  switch (WLANStandard)
    case 0 % 802.11a
      T_20MHz_sifs = 16;
    case 1 % 802.11g-only
      T_20MHz_sifs = 10;
    case 2 % 802.11g-mixed
      T_20MHz_sifs = 10;
  endswitch

  %======> SLOT duration depends standard and channel width
  %         but cannot be computed from the scaling ratio.
  if (WLANStandard == 2)
    % if there is no legacy, 11g BSS can 
    % employ the slow (long) slot duration 
    % IEEE 802.11-2012 standard, section 19.4.5 
    TslotMicroSecs = 20;
  else 
    % In case of 11a or 11g-only, short preamble
    % is enabled for 20 MHz i.e. 9 us. See sections 18.3.8.7
    % and 19.4.5 of IEEE 802.11-2012 standard, respectivelly.
    switch (channelWidth)
      case 20
       TslotMicroSecs = 9;          
      case 10
       TslotMicroSecs = 13;
      case 5
       TslotMicroSecs = 21;      
    endswitch
  endif

  %======> Symbol duration time for the standard 20 MHz chanel
  T_20MHz_symbol = 4;

  %======> Resizing timing parameters according to channel width (scale)
  TsymbolMicroSecs = T_20MHz_symbol * scalingRatio;
  TsifsMicroSecs = T_20MHz_sifs * scalingRatio;
  TdifsMicroSecs = 2*TslotMicroSecs + TsifsMicroSecs;  
endfunction


function [numberOfOFDMSymbols] = computeNumberOfOFDMSymbols(psduDataBits, dataRateMbps)
  % 1. Calculate the number of bits per OFDM symbol at rate dataRateMbps
  %  See Section 18.3.2.3 "Modulation-dependent parameters" Table 18-4 
  %  IEEE Std 802.11-2012. Corresponds to N_{DBPS} in the table
  numberOfDataBitsPerSymbol = 4 * dataRateMbps;

  % 2. Calculating the number of OFDM symbols required to transmit psduDataBits bits
  %  (Section 18.3.5.4 "Pad bits (PAD)" Equation 18-11; IEEE Std 802.11-2012)
  % Note that the pad bits are implicit in the ceil function.
  numberOfOFDMSymbols = ceil (psduDataBits / numberOfDataBitsPerSymbol);
endfunction

% see section 'Pad bits (PAD)' (number 18.3.5.4 or 17.3.5.3 in 802.11-2012, 
% 802.11-2007, respec.)
function [servicePhyBits, tailPhyBits, padPhyBits] = retRateDependentPhyOverHeadBits(psduBits, dataRateMbps)
  servicePhyBits =  16;
  tailPhyBits = 6; 
endfunction

function [tpOverheadBits]=retTransportProtocolHeaderBits(transportProtocol)
  switch (transportProtocol)
    case 0
      tpOverheadBits = 64; % UDP Header in bits: 8 bytes * 8
    case 1 
      tpOverheadBits = 160; % TCP Header (without options) in bits: 20 bytes * 8
  endswitch
endfunction

function [ipOverheadBits]=retIPHeaderBits(ipProtocol)
  switch (ipProtocol)
    case 0
      ipOverheadBits = 160; % IPv4 Header (without options) in bits: 20 bytes * 8
    case 1 
      ipOverheadBits = 320; % IPv6 Header (without options) in bits: 40 bytes * 8
  endswitch
endfunction

function [llcOverheadBits]=retLLCHeaderBits()
  llcOverheadBits = 64; % 8 bytes * 8
endfunction

function [macHeaderBits, macFCSTrailerBits]=ret80211MACOverheadsBits()
  macHeaderBits = 192; %24 bytes * 8;
  macFCSTrailerBits = 32; % 4 bytes * 8
endfunction





function totalSymbols = numberOfSymbols(totalBits, dataRateMbps)
  % Computing the Ndbps (as established by IEEE 802.11 standards)
  % At least in 802.11g/a, a single symbol encodes 4xR bits
  % given a modulation R = {6,9,12,18,24,48,54} 
  % (IEEE Std 802.11-2007, section 17.3.2.2, table 17-3)
  dataBitsPerOFDMSymbol = nominalRate * 4;

  % ceil means: if we need n<1 symbol, we need 1 symbol
  totalSymbols = ceil(totalBits / dataBitsPerOFDMSymbol);
endfunction




%  /* Explanação pacotes
%   * ======================================================> 1.PPDU = PSDU + PHY Headers/Overheads
%   * | |=========================================> 1.1 PSDU = the payload of the PHY Layer i.e. the MPDU
%   * | |==========================> 1.1.1 MPDU = MAC Headers/overheads + MDSU
%   * | |  |===================> 1.1.1.1 MAC Headers/overheads
%   * | |  |  |==========> MAC header  = 24 bytes
%   * | |  |  |==========> MAC Frame Check Sequence (trailer) = 4 bytes
%   * | |  |===================> 1.1.1.2 MSDU i.e. the payload of the MAC frame
%   * | |  |  |==========> Application layer Payload = 1460
%   * | |  |  |==========> UDP/TCP header = 8 bytes
%   * | |  |  |==========> IP header = 20 bytes
%   * | |  |  |==========> LLC header = 8 bytes
%   * |=========================================> 1.2 PHY Headers/Overheads
%   * | |==========> 1.2.1 PHY PLCP Preamble = (20/BW)*16 microseconds
%   * | |==========> 1.2.2 PHY PLCP Header = (20/BW)*4  microseconds
%   * OBS: vide WifiPhy::GetPlcpPreambleDurationMicroSeconds e WifiPhy::GetPlcpHeaderDurationMicroSeconds   
%   * Total = 1524 corresponde a MACTotalDataFrameSizeInBits no script octave (tudo.m)
%   */

function [phyHeaderDurationMicroSecs, phyPayloadDurationMicroSecs] = \
         calculatePPDUDurationMicroSecs(a,b)
phyHeaderDurationMicroSecs = a;
endfunction

% Compute the Physical Protocol Data Unit in bytes, i.e. PSDU (PHY Payload) + PCLP Header + PLCP Preamble
function [psduBytes]=calculatePSDUSize(appLayerPayLoadBytes,transportProtocol)
  % the MPDU (MAC payload + headers) is the PHY payload, i.e. the PSDU
  psduBytes = calculateMPDUSize(appLayerPayLoadBytes,transportProtocol);
endfunction

% Compute the time overhead (in micro-seconds) 
% introduced by the PHY Layer, i.e. plcpPreamble, plcpHeader 
function [plcpHeaderMicroSecs, plcpPreambleMicroSecs, signalExtensionMicroSeconds]=calculatePhyOverhead(WLANStandard)
  % The 802.11a/g PHY headers are in BPSK 1/2 which 
  % yields 6 Mbps in 20 MHz channels
  plcpHeaderMicroSecs =  4;
  plcpPreambleMicroSecs = 16;
  switch (WLANStandard)
  case 0
    % 802.11a
    signalExtensionMicroSeconds = 0;
  case 1
    % 802.11g
    signalExtensionMicroSeconds = 6;
  case 2 
    %not working
  endswitch
endfunction


% Compute the MAC Protocol Data Unit in bytes, i.e. MSDU + MAC Header + Trailler
function [mpduBits]=calculateMPDUSize(appLayerPayLoadBytes,transportProtocol)
  % calculate the MAC payload
  msduBits = calculateMSDUSize(appLayerPayLoadBytes,transportProtocol);
  macHeaderBits = 192; %24 bytes * 8;
  macTrailerFCSBits = 32; % 4 bytes * 8
  mpduBits = msduBits + macHeaderBits + macTrailerFCSBits;
endfunction

% Compute the MAC payload (i.e. the MSDU) in Bytes
function [msduBits]=calculateMSDUSize(appLayerPayLoadBytes,transportProtocol)
  msduBits = appLayerPayLoadBytes*8;
  switch (transportProtocol)
    case 0
      % UDP
      msduBits += 64; % 8 bytes * 8
    case 1 
      % TCP (standard without optional fields)
      msduBits += 160; % 20 bytes * 8
  endswitch
  msduBits += 64; % LLC Header: 8 bytes * 8
endfunction

