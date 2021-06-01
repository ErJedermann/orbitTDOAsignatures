# orbitTDOAsignatures
This repository about 'Orbit-based Authentication Using TDOA Signatures in Satellite Networks' belongs to the paper with the same name. Available at <TODO>.

ABSTRACT
Given the nature of satellites orbiting the Earth on a fixed trajectory, in principle, it is interesting to investigate how this invariant can be exploited for security purposes. In particular, satellite orbit information can be retrieved from public databases. Using time difference of arrival (TDOA) measurements from multiple receivers, we can check this orbit information against a corresponding TDOA-based signature of the satellite. In that sense, we propose an orbit-based authentication scheme for down-link satellite communications in this paper. To investigate the properties and fundamentals of our novel TDOA signature scheme we study two satellite systems at different altitudes: Iridium and Starlink.

Clearly, many challenging questions with respect to the feasibility and effectiveness of this authentication scheme arise; to name some: how many receivers are necessary, how should they be distributed, and how many consecutive measurements do we need for the TDOA signatures. We address these questions by a full factorial experimental design using a simulation framework, we developed for that purpose. Besides a deep understanding about the effects of the major factors on the authentication performance, we find that in adequate configurations, even under a versatile attacker, the orbit-based authentication scheme is able to achieve low false authentication rates well below 1% at false rejection rates of about 2%, for both, Iridium and Starlink satellites.


This is the code that was used to simulate the satellites and their TDOA signatures. 
Required software:
- Python: numpy, sgp4 (https://pypi.org/project/sgp4/), astropy (https://www.astropy.org/), plotly (https://plotly.com/python/)
- R
Usage:
- The first experiment:
	- This experiment generates the data, used for the ANOVA and figure 2.
	- Execute the 'execute_experiment1.py'. In this file the variable 'ir_sats' (line 14) can be set to True/False if Iridium (True) or Starlink (False) satellites should be used.
	- The raw data is stored in "./data/experiment1/run_iridium1/" (line 16) or "./data/experiment1/run_starlink1/" (line 27). For multiple parallel executions it is recommended to change the numbers in the lines, so that results of other executions are not overwritten.
	- For the paper 6 executions with Iridium satellites and 6 executions with Starlink satellites were made. The generated raw data can be found in the folders "./data/experiment1/run_{iridium/starlink}{1-6}/".
- Make figure 2:
	- (Figure 1 in the paper is a static image, that explains the concept of orbit-based authentication using TDOA signatures. Figure 2 is the first graph one generated from data.)
	- Execute the 'make_figure2.py'. 
	- It takes the data stored in "./data/experiment1/run_{iridium/starlink}{1-6}/" and creates the contour plots, that shows the impact of distribution diameter, numer of messages and number of receivers on the F-beta score. The fourth figure (distribution diameter and numer of messages of starlink satellites) is used in the paper.
- ANOVA:
	- To make the ANOVA-table, first execute 'prepare_table_1.py'. It combies the data from the first experiment.
	- Than execute 'ANOCA_combined.R' to get the ANOVA-table.
- The second experiment:
	- This experiment generates the data, used for figure 3.
	- Execute the 'execute_experiment2.py'. Again, this file has a variable 'ir_sats' (line 14) that can be set to True/False if Iridium (True) or Starlink (False) satellites should be used.
	- The raw data is stored in "./data/experiment2/iridium/" (line 16) or "./data/experiment2/starlink/" (line 27). (In the current setup no parallel executions are planned since it is much faster than experiment 1 (still takes a few hours).
- Make figure 3: 
	- Execute the 'make_figure3.py'. 
	- It takes the data stored in "./data/experiment2/{iridium/starlink}/" and creates the line graphs, that shows the false-rejection-rate and false-acceptance-rate of the algorithms. Both generated graphs are used in the paper.


