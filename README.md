# DNSM-special-distribution
1.the allfuc.py is required since all the formula is inside

2. edit the ppar.py file to set parameters(including M1, M2, maxwellian dispersion, ...)
  
3.run ppar.py it will output a par.csv file
  
4.all the codes bellow require allfun.py and par.csv to run

5.then ip.py to get initial position (ip.py are required to read par.py) 
  it will also make a new direction named as the time sending the ip.py, all the file will put in this file
  ip.py will output a IRZ.csv as the initial R and Z as the initial position

6.then run isv.py to get initial system velocity

7.then run to.py inorder to get track orbital(this takes most of the simulation time)
