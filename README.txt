README

==========================================================================
Setting up a new cruise to process
==========================================================================
The code looks for three folders to put data in:
- Processed Data & Plots for each YEAR_DAY(for now set to: PROCESSED/ )
- A Log file for each YEAR_DAY (for now set to: LOGS/ )
- A directory for saving intermittant files during processing that will get 
overwritten by the next YEAR_DAY processed (for now set to: TEMP/ )
	
1.  Create a directory for each of these, if it doesn't already exist
and update readIngestParameters.m to reference the full path for each of these files.
	
2.  Make sure prepacs.exe is executable by your system.  You might need to 
put it on your path.
	
3.  Update readIngestParameters.m to hold all the variables you want for 
processing
		
	Change the paths for your DATA/ LOGS/ PROCESSED/ and TEMP/ directories. 
	Change other variables based on how you want to process the data

4.  Update the import methods for each data type you are using (ac, tsg, flow, etc.)
(See next section)
	
5.  Run NAAMES_IngestManager.m

6.  Run PreProcessingManger.m

7.  Run ProcessingManager.m

8.  Run OutputManager.m

9.  Once you have the code working for the new cruise, you can continue to run the code manually or all at once with a script:

1) block by block, one after the other:
    1.  Set all processing parameters in "readIngestParameters.m"
    2.  Set the YEAR_DAY in "readIngestParameters.m"
    3.  Open IngestManager, set RUN_FROM_SCRIPT = false, and run it.
    4.  Open PreProcesssingManager, set LOAD_FROM_DISK = false, and run it.
    5.  Open ProcessingManager, set LOAD_FROM_DISK = false, and run it.
    6.  Open Output Manager and run it.

2) automatically with the script "RUNLeg.m"
    1.  Set all processing parameters in "readIngestParameters.m"
    2.  Change the YEAR_DAY to be the one you want in RUNLeg
    3.  Change RUNFROMSCRIPT in IngestManager to TRUE
    4.  Make sure LOADFROMDISK is set to FALSE in PreProcessingManager, 
        ProcessingManager, and OutputManager



==========================================================================
Setting up data files
==========================================================================
In this code, the data sources are decoupled from the data types. This means
that different cruises can have the data stored in different formats and 
very little code needs to be changed.  

Before any code can be run, the import functions for importing your raw data 
need to be created.  You can do this manually, or automatically in Matlab.
Four examples are given:
    importfileACS     - the script used for TaraMed to import data output 
                        by prepacs.exe for .bin AC files
    importfileACS_dat - the script used for TaraMed to import data output
                        by prepacs.exe for .dat AC files
    importfileFlow    - the script used for TaraMed to import flow files
                        containing flow and valve data
    importfileTSG     - the script used for TaraMed to import TSG files
                        containing 
It is the job of the "importfiles" to know the file layout of the file and 
identify the data correctly.  The easiest thing might be to edit one of the existing
import files to fit your data.

When you run IngestManager, the objects that call the import methods and assign the data appropriately are:
    ACFileLoader
    GPSFileLoader
    TSGFileLoader

It is the job of the FileLoaders to call the correct import method and assign 
the raw data appropriately to the correct data object

After running IngestManager, the assigned data is now in the data objects:
    ACData
    TemperatureData
    SalinityData
    GPSData
    FlowData
    ValveData

It now no longer matters how the raw data was initally stored, it is now in 
a uniform data structure.

==========================================================================
General code organization
==========================================================================
The AC processing code consists of the following types of code:
A)  Scripts to hold the processing logic
B)  Data objects to hold the data, and the functions to process the data
C)  Output data files, plots, and a processing log

Data is processed by YEAR_DAY - a 24 hour period of data which is part of a 
CRUISE_LEG which is part of a CRUISE.

==========================================================================
What are the main processing files?
==========================================================================
The code specifically contains the following scripts:

1)  RUN_LEG.m              - a script to run one/more days of code all at once
2)  readIngestParameters.m - a matlab function, editable to set all the 
                             parameters needed for processing
3)  IngestManager          - a script to call the objects and methods needed  
                             for ingesting raw data files
4)  PreProcessingManager   - a script to call the objects and methods needed
                             for finding the TSW/FSW transitions in a/c data,
                             flow, and valve data; checking transitions 
                             (manually if necessary); synchronizing a/c data
                             with flow/valve data; separating TSW/FSW
5)  ProcessingManager      - a script to call the objects and methods needed
                             to process the separated a/c TSW/FSW data
6)  OutputManager          - a script to do the final processing needed to 
                             output various file formats (SeaBASS) and final
                             plots


==========================================================================
Editing TSW/FSW transitions manually
==========================================================================
1.  Open PreProcessingManager
2.  Set MANUAL_MODE = true;
3.  Run the code IN SECTIONS.  Follow the directions when you get to the manual 
    sections.
4.  If you need to re-run the code with the data you set manually, set 
    REPROCESS_MANUAL = true; (keep MANUAL_MODE = true;)