# RF-sheet-resistance

This code can be used to extract the sheet resistance of conductive films (e.g., conductive fabric, MXene, etc.).

This code needs two-port microstrip transmission line S-parameters in the form of .s2p files (dB/angle degrees). 

# The following parameters need to be edited into the main.mlx file: 
1. filename: Name of the .s2p file from the network analyzer,
2. CL: connector loss in dB (optional, default 0), 
3. dl: connector length (meters), from the point of calibration to the device under 
test (DUT)/transmission line top, 
4. ub, lb: upper and lower bound of optimization, doesn't need to change, 
5. m: number of the frequency of interest in the freq vector. The sheet resistance value is calculated
at this frequency, 
6. len: length of the transmission line / DUT in meters,
7. w: width of top layer in meters, 
8. rad_eff: radiation efficiency (%) of the transmission line, typically less than 5% for electrically small
transmission lines. Can be derived from simulation.

If you find this code useful and use it in your research, please cite the following paper:
```
@INPROCEEDINGS{9329989,
  author={Md Abu Saleh Tajin and Ariana S. Levitt and Yuqiao Liu and Chelsea E. Amanatides and Caroline L. Schauer and Genevieve Dion and Kapil R. Dandekar}, 
  booktitle={2020 IEEE International Symposium on Antennas and Propagation and North American Radio Science Meeting}, 
  title={Extraction of Knitted RFID Antenna Design Parameter from Transmission Line Measurements}, 
  year={2020},
  volume={},
  number={},
  pages={1551-1552},
  doi={10.1109/IEEECONF35879.2020.9329989}}
```  
