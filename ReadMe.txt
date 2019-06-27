*********************************
Use linux cmd to run it !
*********************************

1, Put the all the files into the same path.

2, put the .hbf files ( Px,Py you want to solve) into the same path

3, Configuration file modification,

    change the "rx.hbf" into the name of file for Pentex.hbf
    change the "ry.hbf" into the name of file for Pentey.hbf

    The parametre for "solution" is  the name of the file of estimated solution
    or we put the name of one of the hbf files Pentex.hbf/Pentey.hbf

4, cd into the path, type in cmd, 
    ./feelpp_h3_surface_rec --config-file rec.cfg
     to run the application.

5, The result is stored as "z.hbf" in the path indicated in cmd window.

6, if "bash" problem occurs, use  the command:
	chmod 777 feelpp_h3_surface_rec 