Welcome to a Python simulation of ideal gas law. The work contains three files: "Main", "Classes" and "tests". "Main" contains a code which will output data based on a number of parameters that can be set. "Classes" contains the module with classes needed to run the simulation. Please do not alter anything inside it! "tests" contains snippets of code that were used to design the bulk of the code. These may not be functional anymore, and are there just to demonstrate the development process.


To run the simulation, please first run the "Classes" file.

Note: Before you run the "Main" file, be advised that this program creates four small external files (in binary) in the directory where the code is run. They are used to record the data and are overwritten on every runthrough. When you are done with the simulation, remember to delete them.

Next run the "Main" file. You can change the default parameters: the number of particles (labelled "Number_of_Balls", default 20) and number of frames (labelled "Number_of_Frames", default 1000). The code will output the volume, pressure and temperature in the console. To obtain graphs, please toggle the desired options from "False" to "True" and run again. There are three: A distribution of speeds (maxwell-boltzmann), instantaneous pressure and total energy throughout the simulation. You can also switch the visualisation off by toggling "animate" to "False".

Exceeding the allowed number of particles (80) will raise an exception. To overrun this, please increase the "volume" value in the "Simulation" class inside the "Classes" and get rid of the exception. Please note that for large numbers of particles a sufficiently large container is needed. It is perhaps better to just increase the number of frames to obtain better statistical results.

