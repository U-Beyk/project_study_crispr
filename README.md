# Project study
Project study of my Master of science.

To run this tool please execute the main.py from the folder it is located in.   
Also you need to  download the sql dump of the cripsr cas database:  
https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=ccpp_db.zip  
Save the sql dump in "crispr_cas_db/db_parser/Crispr_Cas_Database_SQL_Dump.sql".

You also need to download RNAmotiFold from here:  
https://github.com/RNABioInfo/RNAmotiFold  
Save it in the same folder the main.py is located in. 

Be aware, that it is possible that the RNAmotiFold algorith might not work due to changes in the tool.

Many of the classes in this project are generic, abstract, or serve as parent classes for others. This design choice was made to promote code reuse. Much of the code is expected to be refactored over time to be more generic, rather than hardcoded. The only part that will remain implementation-specific is the database processing logic, which needs to be rewritten for each distinct database, and a small part of the analysis.  
Additionally, this project currently lacks testing and exception handling. Since it was primarily developed for internal data analysis and not intended as a widely used tool, these aspects were initially deprioritized. However, they may be introduced in a future version of the project. The documentation is also kept to a minimum.