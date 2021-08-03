# Growth-Distribution-and-Dynamic-Inefficiency-in-Turkey

GDETR_main.m generates the figures and tables in the paper.

GDETR_matlab.xls is the main data file (see below).

The original data file organized and used by Altug et al. (2008, EREH) is
located at 

   http://myweb.sabanciuniv.edu/alpayf/files/2010/04/erehdata.xls

 and has 13 data columns. The GDETR_matlab.xls file is organized in the 
 following way:
 
 (*) From Altug et al. (2008, EREH):

       A1:A83 year              (A3:A85 in the original erehdata.xls file)
       B1:B83 Real GDP          (B3:B85 in the original erehdata.xls file)
       C1:C83 Capital stock     (E3:E85 in the original erehdata.xls file)
       D1:D83 Employment        (H3:H85 in the original erehdata.xls file)
       E1:E83 Human capital     (K3:K85 in the original erehdata.xls file)
 
 (*) From the Penn World Tables v.9:

       F1:F83 Human capital (PWT 9)
              F1:F27  1923-1949  no data ( = NaN )
              F28:F83 1950-2005  Turkey's "hc" data from the PWT 9
