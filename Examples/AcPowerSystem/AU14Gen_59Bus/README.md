# AU14Gen_59Bus

This example is based on reports of IEEE-PES-TR18 entitled "Benchmark Systems for Small-Signal Stability Analysis and Control" from IEEE PES Task Force on Benchmark Systems for Stability Controls. 

A summary of this report has been published in the paper  entitled "Benchmark Models for the Analysis and Control of Small-Signal Oscillatory Dynamics in Power Systems".

On the link (https://cmte.ieee.org/pes-psdp/benchmark-systems-2/), reports and data files for six benchmark systems have been provided. 
	
The current example is based on the Australian 14 generator system report entitled "Simplified 14-Generator Model of the South East Australian Power System" and its accompanying data files. 

In the AU14GenModelData_Ver04.zip, load flow files are available for different load cases.

We prepared this example based on load case 1 data, and it is possible to prepare files for other load cases.

The example has been prepared using LF_Case01_R4_S.raw, which is available in AU14GenModelData_Ver04.zip.

# Notes

## Bus data (LF_Case01_R4_S.raw file: BUS DATA and LOAD DATA)

Bus types are based on the data file. 
	In the data file, nineteen buses 1,3,4,5,6,7,20,21,32,35,36,37,38,46,51,52,53,57, and 59 are generator buses.
	Bus No. 1 has been considered as ref = 1. 
	Other generator buses as PV = 2. 
	Remained buses (niether ref, nor PV) as PQ = 3.
	Main 14 generator buses are: 1,3,4,5,6,20,21,35,36,37,38,51,52, and 53.
	In buses 7,32,46,57, and 59, Pg is 0.

Parameters (for power flow analysis) in Bus data (Vsp (pu), theta (rad), PLi (pu), and QLi (pu)) are based on the data file.

For area, original bus numbers' first digit, refers to the corresponding area. 
	Australian 14 generator system covers 5 areas.
	In this example, we converted three digit bus numbers to 1-59 buses. 
	Also, we consider whole system as a single area, so we set all areas as 1 (1 ref bus; 1 area).
  
## Bus data for generators (LF_Case01_R4_S.raw file: GENERATOR DATA)
	
For PGi (pu) and QGi (pu), the data file values were set.

For Qmin (pu) and Qmax (pu), two typical large numbers were set (-99,99).

## Apparatus data

For generator parameters (J (pu), D (pu), wL (pu), and R (pu)), values from the report "Simplified 14-Generator Model of the South East Australian Power System" were set.

In the report, values for D and R are equal 0. 

We set them as 0.25*J and 0.0001, respectively.

It is okay to consider original values.

## NetworkLine_IEEE data (LF_Case01_R4_S.raw file: BRANCH DATA and TRANSFORMER DATA)
	
For transmission lines, r, x, and b values from BRANCH DATA and ratio = 1 were entered.

For transmission line (29-30), x value was -0.03370, which was set positive.

For transformers, x and ratio values from TRANSFORMER DATA were used.

Bus numbers from LF_Case01_R4_S.raw file were changed as follows:

(101,1)
(102,2)
(201,3)
(202,4)
(203,5)
(204,6)
(205,7)
(206,8)
(207,9)
(208,10)
(209,11)
(210,12)
(211,13)
(212,14)
(213,15)
(214,16)
(215,17)
(216,18)
(217,19)
(301,20)
(302,21)
(303,22)
(304,23)
(305,24)
(306,25)
(307,26)
(308,27)
(309,28)
(310,29)
(311,30)
(312,31)
(313,32)
(314,33)
(315,34)
(401,35)
(402,36)
(403,37)
(404,38)
(405,39)
(406,40)
(407,41)
(408,42)
(409,43)
(410,44)
(411,45)
(412,46)
(413,47)
(414,48)
(415,49)
(416,50)
(501,51)
(502,52)
(503,53)
(504,54)
(505,55)
(506,56)
(507,57)
(508,58)
(509,59)

## Basic data

System frequency is 50 Hz.

# References

[1] C. Canizares et al., "Benchmark Models for the Analysis and Control of Small-Signal Oscillatory Dynamics in Power Systems," in IEEE Transactions on Power Systems, vol. 32, no. 1, pp. 715-722, Jan. 2017, doi: 10.1109/TPWRS.2016.2561263.

[2] Benchmark Systems for Small-Signal Stability Analysis and Control. Available online at https://cmte.ieee.org/pes-psdp/benchmark-systems-2/

[3] M.J. Gibbard & D.J. Vowles, "Simplified 14-Generator Model of the South East Australian Power System", Revision 4, School of Electrical & Electronic Engineering, The University of Adelaide, South Australia, 14 June 2014. References [3] is source of data.