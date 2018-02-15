# TrackletEmulation
Updated 11/10/2017: now separates cluster-finding into two layers to more accurately emulate the firmware.
L1_cluster does phi-clustering in each etabin, and L2_cluster stitches everything together.

11/21/2017: updated to use an algorithm closer to the firmware algorithm. Also some bug fixes.

02/15/2018: Setup a rough version for Code Sharing with RAL. Includes:
	Trigger Turn On Code
	Truth Matching for Efficiency SUmmary
	L1 Rates in Min Bias Sample
	Copy of L1TrackNtupler with Fast Jets
	(To compile the necessary C++ code and ROOT Macros)
	make all 	
	
