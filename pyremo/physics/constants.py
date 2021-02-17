# flake8: noqa
"""collection of common Remo physics constants.
"""

# descritions from REMO
#!C**
#!C**   BESCHREIBUNG DER VARIABLEN:
#!C**   T0       :   NULLPUNKT DER TEMPERATUR               (K)
#!C**   R        :   GASKONSTANTE FUER TROCKENE LUFT        (J/(KG*K))
#!C**   RD       :   GASKONSTANTE FUER WASSERDAMPF          (J/(KG*K))
#!C**   WCP      :   SPEZ.WAERME  FUER TROCKENE LUFT        (J/(KG*K))
#!C**   WLK      :   VERDAMPFUNGSWAERME                     (J/KG)
#!C**   WLF      :   GEFRIERWAERME                          (J/KG)
#!C**   WLS      :   SUBLIMATIONSWAERME                     (J/KG)
#!C**   G        :   ERDBESCHLEUNIGUNG                      (M/S**2)
#!C**   RERD     :   MITTLERER ERDRADIUS                    (M)
#!C**   STAG     :   MITTLERER STERNTAG                     (S)
#!C**   RHF      :   DICHTE DES WASSERS                     (KG/M**3)
#!C**   SIGMA    :   BOLTZMANN-KONSTANTE                    (W/(M**2*K**4)
#!C**   SOKO     :   SOLARKONSTANTE                         (W/M**2)
#!C**
#!C**
#!C**   BESCHREIBUNG DER VARIABLEN:
#!C**   B1       :   PARAMETER ZUR BERECHNUNG DES SAETTIGUNGS-  (PA)
#!C**   B2W      :   DAMPFDRUCKES UEBER WASSER (W) UND EIS (E)  (--)
#!C**   B2E      :                 "                            (--)
#!C**   B3       :                 "                            (K)
#!C**   B4W      :                 "                            (K)
#!C**   B4E      :                 "                            (K)
#!C**   UC1      :   PARAMETER ZUR BERECHNUNG DES BEDECKUNGS-   (--)
#!C**   UC2      :   GRADES MIT WOLKEN IM UNGESAETTIGTEN FALL   (--)
#!C**   UCL      :                 "
#!C**   RHDE     :   PARAMETER ZUR BERECHNUNG DES SCHNEEBEDECK- (--)
#!C**                TEN BODENS
#!C**   AKS2     :   PARAMETER BEI DER HORIZONTALDIFFUSION 2.O. (--)
#!C**   AKS4     :   PARAMETER BEI DER HORIZONTALDIFFUSION 4.O. (--)
#!C**   AKT      :   VON KARMAN-KONSTANTE                       (--)
#!C**
#!C**   VERFASSER:   D.MAJEWSKI


# these are from druint
T0 = 273.15
R = 287.05
RD = 461.51
WCP = 1005.0
WLK = 2.501e6
WLF = 0.334e6
WLS = 2.835e6
G = 9.80665
RERd = 6371.229e3
STAg = 86164.09054
B1 = 610.78
B2W = 17.2693882
B2E = 21.8745584
B3 = 273.16
B4W = 35.86
B4E = 7.66
RHDe = 1.0 / 2.5


# from preprocessor RemapToRemo

# T0 = 273.16
# R = 2.8705E2 differs
# RD = 4.6151E2
# WCP = 1.005E3
# WLK = 2.501E6
# WLF = 0.334E6
# WLS = 2.835E6
# G = 9.80665
# RERd = 6.371229E6
# STAg = 86164.09054


PI = 3.141592654
RADdeg = 57.29577951
DEGrad = 0.0174532925
RDRd = R / RD
RDDrm1 = RD / R - 1.0
EMRdrd = 1.0 - RDRd
#!
#!     DEFAULT-WERTE FUER /COMPAR/
#!
# B1 = 610.78
# B2W = 17.2693882
# B2E = 21.8745584
# B3 = 273.16
# B4E = 7.66
# B4W = 35.86
