Guide lines on using the virus competition model:

Content:
	1. running simulations
	2. data storage
	3. data analysis


1. running simulations 

	There's 2 types of models: comp and meta

	1) comp model

		a) explanation
		comp model is a single population model where you can simulate 
		segmented and non-segmented viruses to compete with each other. They don't have to
		compete if you set the parameter 'N1r' to 1 or 0.
		comp model has 3 different versions.
		ver 1 is a prototype, and practically useless. I first made the model with python so it has .py and .c (comp1.1)
		ver 2 has flexible population size with carrying capacity K (comp1.2). The latest one 1.2.3 is only logically sound and work fine with commandx.py 
		ver 3 is a WF model with constant population (comp1.3). 1.3 will only work perfectly. 1.3.2 is garbage now and 1.3.3 is iffy.
		comp2.3 and comp3.3 adds 3 segments and 8 segments respectively but needs a lot more work to run fine.

		b) how to run
		same as running meta model.
		To run simply type: 
		```python command.py comp(version).c (code)``` 

	2) meta model

		A) explanation
		More detailed explanation of the model is described in the Yoon2019.docx
		meta model is a multi metapopulation model where you can simulate the competition between segmented and non-segmented viruses in multiple demes.
		Each step, the viruses can migrate to the migration pool and transmit to other hosts with some transmission rate. 

		The model consists of 4 or 5 steps and a recording step:

		a) mutation:
			In this step, viruses in each class (class denote a set of virus with equal amount of deleterious mutations in each segments)
			goes through mutation according to the set mutation rate.

		b) reassortment:
			In this step, segmented virus goes through reassortment process, where their two segments get mixed up.

		c) repoduction:
			Viruses reproduce reflecting their fitness determined by their mutational load (amount of mutation they have) and the carrying capacity of a deme.

		d) migration:
			Viruses migrate among the hosts dependent on the set migration and transmission rate. Some proportion of viruses from each deme
			migrate to the pool, and they get equally redistributed to every deme in proportion of transmission rate.

		e) evolution:
			This step is optional depending on whether the user wants to run the simulation where segmented population is present from the start or run the simulation with the segmented population occuring from non-segmented viruses after some number of generation. When the parameter value of 'e' is above zero, the latter mode is activated and the evolution step is run in the simulation. In this step, some number of non-segmented viruses is turned into the segmented viruses in proportion to e.


		B) how to run
		meta model has more easy and convenient method of running the simulation.
		Instead of making multiple commands for different parameter sets in command*.py, here we use a csv file called 'paraminfo.csv' to store all the 
		different parameter combinations. Each row contains a combination of parameters that is used for a simulation that is described on the first column, 'description'. 

		To run the set of parameters on a row. You specify the code on a 'code' column for the row on shell command.

		When making some parameters, theres are some rules to follow:

		-set pop2initlen and pop1initlen to 0. Its just used for parsing in the c file.

		-pop2init and pop1init has to have a number between 0 and 1 followed by '~'sign. Each '(number)~' is whatever segment's initial frequency at a host.
		For example with hostnum=3, pop2init= '0.5~0.4~0.3~', pop1init= '0.5~0.6~0.7~', the three demes are initialized with 0.5, 0.4, and 0.3 of their population with segmented viruses and 0.5, 0.6, and 0.7 with non-segmented viruses. Ideally, the proportion of each type of viruses in a deme should not exceed 1.
		As the number of hosts increases, there's a way to abbreviate the pop#init parameter value. Putting a number after a colon followed by'(proportion)~ ' is equivalent to writing '(proportion)~' the stated number of times. For example '3:1~2:0.2~' is equivalent to '1~1~1~0.2~0.2~'
		Another thing to note is that is a user writes fewer '(proportion)~' than the amount of hosts, the rest of the host that didn't get specified gets '0~'. For example, when hostnum=3 and pop2init='2:1~', then the last host is automatically assigned with '0~'

		-if you want random seed, type 'random'. If you need a consistent random number, put any negative integer below -9223372036854775808.
		To run, simply type:

		```python command.py meta1.1.3.c (code)```

2. data storage
	All data is stored in Reassortment/data/(destination parameter value).
	Each simulation data is stored in a .csv file.
	There are several ways one can record the simulation data.
	-krecord=0 records mean mutational load, whereas krecord=1 records minimum mutational load in the population. krecord=2 records distributional data over mutational load.
	-with timestep=1, every population information is recorded every generation and only every repetition with timestep=0.
	-with SHORTCSV 1, only total population size and mutational load is recorded, and with SHORTCSV 0, every deme's population size and mutational load is also recorded.
	-with EACHSTEPREC 1, data is recorded every step of the simulation.

	In data, 'rep' and 'gen' denote repetition and generation number of the simulation, 'pop1.0' and 'pop2.0' the total population size of the non-segmented and segmented viruses, and 'k1.0' and 'k2.0' the mean mutational load of non-segmented and segmented viruses.
	When the data contain information for every host, the host number is identified by the number after the dot. For example 'pop1.4' denote population size of the non-segmented virus in host 4.
	In histogram data, number after p denote virus type (1 for non-segmented and 0 for segmented) and k denote mutational load. For example 'p1k20' denote number of non-segmented viruses with 20 mutations. 
	

3. data analysis.
	For plotting total mutational load and population size for both types of viruses and plotting histograms of mutation load in certain reps and generations, use data_viewer2.ipynb using Jupyter notebook.






















