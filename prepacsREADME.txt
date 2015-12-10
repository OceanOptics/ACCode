 Command line formats: 
 prepacs devfilename binfilename outfilename 
 prepacs devfilename binfilename outfilename switch 
 
 This program preps an acs binary file to output an 
 ascii tab-delimited data file with one header row. 
 
 If the switch '+wet' is used, devfile info is  
 prepended to the outfile so that it can be 
 replayed in wetview. Note that, because wetview 
 skips the 1st two data packets (and prepacs does 
 not), the acs time as determined by prepacs will 
 be about 500 msec (2 packets' worth) later for 
 a given data line than that from wetview. 
 
 default action (no switch): 
 the first data column is time in msec, followed by 
 first the 'c' acs data channels, then the 'a'  
 channels, each in order of increasing wavelength. 
 
 auxiliary and dark acs values from the binary data 
 can be included in the outputfile if a valid switch 
 is used, although the auxiliary values may be 
 artifactual, depending on the acs configuration. 
 
 valid switches and their actions: 
     +aux  : also output 4 auxiliary channels; 
     +dark : also output 4      dark channels; 
     +all  : also output all 8 extra channels; 
     +wet  : also output all 8 extra channels, and, 
             outfile can be replayed in wetview. 
 
 
 Last updated by Russell Desiderio, 18-Dec-2009. 