Formal model of the 802.11 DCF using UPPAAL
---

This is a formal model of the IEEE 802.11 Distributed Coordination Function. It was created with [UPPAAL](http://uppaal.org).

Some assumptions were made:

- Basic Access only (i.e. no RTS/CTS);
- Saturated medium;
- Start time of the STAs is the same;
- Same Tx time for all STAs;
- No capture effect;
- No EIFS;
- Infinite retries when the transmission fails.

The values of the timing parameters were computed according to the [standard](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6178212) with the help of an octave script from [1].

The complete model is shown [here](https://raw.githubusercontent.com/fabianorosas/csma-ca-uppaal/master/csmaca.png).

--

[1] - Saulo Queiroz and Roberto Hexsel, 2014. Translating Full Duplex into Capacity Gains for the High-Priority Traffic Classes of IEEE 802.11. In Proceedings of the 30th Annual ACM SIGAPP Symposium on Applied Computing (SAC '15). ACM, New York, NY, USA.